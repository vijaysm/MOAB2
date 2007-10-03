/* compile-time settings:

   FORTRAN naming convention
     default      cpgs_setup, etc.
     -DUPCASE     CPGS_SETUP, etc.
     -DUNDERSCORE cpgs_setup_, etc.

   -DMPI             parallel version (sequential otherwise)
   -DCRYSTAL_STATIC  avoid some message exchange at the risk of
                     crashing b/c of insufficient buffer size
   
   -DINITIAL_BUFFER_SIZE=expression
      ignored unless CRYSTAL_STATIC is defined.
      arithmetic expression controlling the initial buffer size for the crystal
      router; this needs to be large enough to hold the intermediate messages
      during all stages of the crystal router
      
      variables that can be used in expression include
         num   - the number of processors
         n     - the length of the global index array

*/

/* default for INITIAL_BUFFER_SIZE */
#ifdef CRYSTAL_STATIC
#  ifndef INITIAL_BUFFER_SIZE
#    define INITIAL_BUFFER_SIZE 2*(3*num+n*9)
#  endif
#endif

/* FORTRAN usage:

   integer hc, np
   call crystal_new(hc,comm,np)  ! get a crystal router handle (see fcrystal.c)

   integer hgs
   integer n, max_vec_dim
   integer*? global_index_array(1:n) ! type corresponding to slong in "types.h"

   call cpgs_setup(hgs,hc,global_index_array,n,max_vec_dim)
     sets hgs to new handle

   !ok to call crystal_done(hc) here, or any later time

   call cpgs_op(hgs, u, op)
     integer handle, op : 1-add, 2-multiply, 3-min, 4-max
     real    u(1:n) - same layout as global_index_array provided to cpgs_setup

   call cpgs_op_vec(hgs, u, d, op)
     integer op : 1-add, 2-multiply, 3-min, 4-max
     integer d    <= max_vec_dim
     real    u(1:d, 1:n) - vector components for each node stored together

   call cpgs_op_many(hgs, u1, u2, u3, u4, u5, u6, d, op)
     integer op : 1-add, 2-multiply, 3-min, 4-max
     integer d : in {1,2,3,4,5,6}, <= max_vec_dim
     real    u1(1:n), u2(1:n), u3(1:n), etc.
     
     same effect as: call cpgs_op(hgs, u1, op)
                     if(d.gt.1) call cpgs_op(hgs, u2, op)
                     if(d.gt.2) call cpgs_op(hgs, u3, op)
                     etc.
     with possibly some savings as fewer messages are exchanged
   
   call cpgs_free(hgs)
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#ifdef USE_MPI
#  include <mpi.h>
#endif

#include "errmem.h"     
#include "types.h"
#include "minmax.h"
#include "tuple_list.h"
#ifdef USE_MPI
#  include "crystal.h"  
#  include "transfer.h"

typedef struct {
  uint np;           /* number of processors to communicate with          */
  uint *target;      /* int target[np]: array of processor ids to comm w/ */
  uint *nshared;     /* nshared[i] = number of points shared w/ target[i] */
  uint *sh_ind;      /* list of shared point indices                      */
  MPI_Request *reqs; /* pre-allocated for MPI calls                       */
  real *buf;         /* pre-allocated buffer to receive data              */
  uint maxv;         /* maximum vector size                               */
} nonlocal_info;

#else
   typedef void crystal_data;
#endif

typedef struct {
  sint *local_cm; /* local condense map */
#ifdef USE_MPI
  nonlocal_info *nlinfo;
  MPI_Comm comm;
#endif
} gs_data;

#define OP_ADD 1
#define OP_MUL 2
#define OP_MIN 3
#define OP_MAX 4
#define OP_BPR 5

/*--------------------------------------------------------------------------
   Local Execution Phases
  --------------------------------------------------------------------------*/

#define DO_SET(a,b) b=a
#define DO_ADD(a,b) a+=b
#define DO_MUL(a,b) a*=b
#define DO_MIN(a,b) if(b<a) a=b
#define DO_MAX(a,b) if(b>a) a=b
#define DO_BPR(a,b) \
  do { uint a_ = a; uint b_ = b; \
       for(;;) { if(a_<b_) b_>>=1; else if(b_<a_) a_>>=1; else break; } \
       a = a_; \
     } while(0)


#define LOOP(op) do { \
  sint i,j; \
  while((i=*cm++) != -1) \
    while((j=*cm++) != -1) \
      op(u[i],u[j]); \
} while(0)
  
static void local_condense(real *u, int op, const sint *cm)
{
  switch(op) {
    case OP_ADD: LOOP(DO_ADD); break;
    case OP_MUL: LOOP(DO_MUL); break;
    case OP_MIN: LOOP(DO_MIN); break;
    case OP_MAX: LOOP(DO_MAX); break;
    case OP_BPR: LOOP(DO_BPR); break;
  }
}

static void local_uncondense(real *u, const sint *cm)
{
  LOOP(DO_SET);
}

#undef LOOP

#define LOOP(op) do { \
  sint i,j,k; \
  while((i=*cm++) != -1) { \
    real *pi=u+n*i; \
    while((j=*cm++) != -1) { \
      real *pj=u+n*j; \
      for(k=n;k;--k) { op(*pi,*pj); ++pi, ++pj; } \
    } \
  } \
} while(0)

static void local_condense_vec(real *u, uint n, int op, const sint *cm)
{
  switch(op) {
    case OP_ADD: LOOP(DO_ADD); break;
    case OP_MUL: LOOP(DO_MUL); break;
    case OP_MIN: LOOP(DO_MIN); break;
    case OP_MAX: LOOP(DO_MAX); break;
    case OP_BPR: LOOP(DO_BPR); break;
  }
}

static void local_uncondense_vec(real *u, uint n, const sint *cm)
{
  LOOP(DO_SET);
}

#undef LOOP

/*--------------------------------------------------------------------------
   Non-local Execution Phases
  --------------------------------------------------------------------------*/

#ifdef USE_MPI

static nonlocal_info *nlinfo_alloc(uint np, uint count, uint maxv)
{
  nonlocal_info *info = tmalloc(nonlocal_info,1);
  info->np = np;
  info->target = tmalloc(uint,2*np+count);
  info->nshared = info->target + np;
  info->sh_ind = info->nshared + np;
  info->reqs = tmalloc(MPI_Request,2*np);
  info->buf = tmalloc(real,2*count*maxv);
  info->maxv = maxv;
  return info;
}

static void nlinfo_free(nonlocal_info *info)
{
  free(info->buf);
  free(info->reqs);
  free(info->target);
  free(info);
}

static void nonlocal(real *u, int op, const nonlocal_info *info, MPI_Comm comm)
{
  MPI_Status status;
  uint np = info->np, i;
  MPI_Request *reqs = info->reqs;
  uint *targ = info->target;
  uint *nshared = info->nshared;
  uint *sh_ind = info->sh_ind;
  uint id;
  real *buf = info->buf, *start;
  { int i; MPI_Comm_rank(comm,&i); id=i; }
  for(i=0;i<np;++i) {
    uint c = nshared[i];
    start = buf;
    for(;c;--c) *buf++ = u[*sh_ind++];
    MPI_Isend(start,nshared[i]*sizeof(real),MPI_UNSIGNED_CHAR,
              targ[i],id,comm,reqs++);
  }
  start = buf;
  for(i=0;i<np;++i) {
    MPI_Irecv(start,nshared[i]*sizeof(real),MPI_UNSIGNED_CHAR,
              targ[i],targ[i],comm,reqs++);
    start+=nshared[i];
  }
  for(reqs=info->reqs,i=np*2;i;--i) MPI_Wait(reqs++,&status);
  sh_ind = info->sh_ind;
# define LOOP(OP) do { \
    for(i=0;i<np;++i) { \
      uint c; \
      for(c=nshared[i];c;--c) { OP(u[*sh_ind],*buf); ++sh_ind, ++buf; } \
    } \
  } while(0)
  switch(op) {
    case OP_ADD: LOOP(DO_ADD); break;
    case OP_MUL: LOOP(DO_MUL); break;
    case OP_MIN: LOOP(DO_MIN); break;
    case OP_MAX: LOOP(DO_MAX); break;
    case OP_BPR: LOOP(DO_BPR); break;
  }
# undef LOOP
}

static void nonlocal_vec(real *u, uint n, int op,
                         const nonlocal_info *info, MPI_Comm comm)
{
  MPI_Status status;
  uint np = info->np, i;
  MPI_Request *reqs = info->reqs;
  uint *targ = info->target;
  uint *nshared = info->nshared;
  uint *sh_ind = info->sh_ind;
  uint id;
  real *buf = info->buf, *start;
  uint size = n*sizeof(real);
  { int i; MPI_Comm_rank(comm,&i); id=i; }
  for(i=0;i<np;++i) {
    uint ns=nshared[i], c=ns;
    start = buf;
    for(;c;--c) memcpy(buf,u+n*(*sh_ind++),size), buf+=n;
    MPI_Isend(start,ns*size,MPI_UNSIGNED_CHAR,targ[i],id,comm,reqs++);
  }
  start = buf;
  for(i=0;i<np;++i) {
    int nsn=n*nshared[i];
    MPI_Irecv(start,nsn*size,MPI_UNSIGNED_CHAR,targ[i],targ[i],comm,reqs++);
    start+=nsn;
  }
  for(reqs=info->reqs,i=np*2;i;--i) MPI_Wait(reqs++,&status);
  sh_ind = info->sh_ind;
# define LOOP(OP) do { \
    for(i=0;i<np;++i) { \
      uint c,j; \
      for(c=nshared[i];c;--c) { \
        real *uu=u+n*(*sh_ind++); \
        for(j=n;j;--j) { OP(*uu,*buf); ++uu, ++buf; } \
      } \
    } \
  } while(0)
  switch(op) {
    case OP_ADD: LOOP(DO_ADD); break;
    case OP_MUL: LOOP(DO_MUL); break;
    case OP_MIN: LOOP(DO_MIN); break;
    case OP_MAX: LOOP(DO_MAX); break;
    case OP_BPR: LOOP(DO_BPR); break;
  }
# undef LOOP
}

static void nonlocal_many(real **u, uint n, int op,
                          const nonlocal_info *info, MPI_Comm comm)
{
  MPI_Status status;
  uint np = info->np, i;
  MPI_Request *reqs = info->reqs;
  uint *targ = info->target;
  uint *nshared = info->nshared;
  uint *sh_ind = info->sh_ind;
  uint id;
  real *buf = info->buf, *start;
  { int i; MPI_Comm_rank(comm,&i); id=i; }
  for(i=0;i<np;++i) {
    uint c, j, ns = nshared[i];
    start = buf;
    for(j=0;j<n;++j) {real*uu=u[j]; for(c=0;c<ns;++c) *buf++=uu[sh_ind[c]];}
    sh_ind+=ns;
    MPI_Isend(start,n*ns*sizeof(real),MPI_UNSIGNED_CHAR,targ[i],id,comm,reqs++);
  }
  start = buf;
  for(i=0;i<np;++i) {
    int nsn = n*nshared[i];
    MPI_Irecv(start,nsn*sizeof(real),MPI_UNSIGNED_CHAR,
              targ[i],targ[i],comm,reqs++);
    start+=nsn;
  }
  for(reqs=info->reqs,i=np*2;i;--i) MPI_Wait(reqs++,&status);
  sh_ind = info->sh_ind;
# define LOOP(OP) do { \
    for(i=0;i<np;++i) { \
      uint c,j,ns=nshared[i]; \
      for(j=0;j<n;++j) { \
        real *uu=u[j]; \
        for(c=0;c<ns;++c) { OP(uu[sh_ind[c]],*buf); ++buf; } \
      } \
      sh_ind+=ns; \
    } \
  } while(0)
  switch(op) {
    case OP_ADD: LOOP(DO_ADD); break;
    case OP_MUL: LOOP(DO_MUL); break;
    case OP_MIN: LOOP(DO_MIN); break;
    case OP_MAX: LOOP(DO_MAX); break;
    case OP_BPR: LOOP(DO_BPR); break;
  }
# undef LOOP
}
#endif

/*--------------------------------------------------------------------------
   Combined Execution
  --------------------------------------------------------------------------*/

void gs_op(real *u, int op, const gs_data *data)
{
  local_condense(u,op,data->local_cm);
#ifdef USE_MPI
  nonlocal(u,op,data->nlinfo,data->comm);
#endif
  local_uncondense(u,data->local_cm);
}

void gs_op_vec(real *u, uint n, int op, const gs_data *data)
{
#ifdef USE_MPI
  if(n>data->nlinfo->maxv)
    fail("%s: initialized with max vec size = %d,"
         " but called with vec size = %d\n",__FILE__,data->nlinfo->maxv,n);
#endif
  local_condense_vec(u,n,op,data->local_cm);
#ifdef USE_MPI
  nonlocal_vec(u,n,op,data->nlinfo,data->comm);
#endif
  local_uncondense_vec(u,n,data->local_cm);
}

void gs_op_many(real **u, uint n, int op, const gs_data *data)
{
  uint i;
#ifdef USE_MPI
  if(n>data->nlinfo->maxv)
    fail("%s: initialized with max vec size = %d,"
         " but called with vec size = %d\n",__FILE__,data->nlinfo->maxv,n);
#endif
  for(i=0;i<n;++i) local_condense(u[i],op,data->local_cm);
#ifdef USE_MPI
  nonlocal_many(u,n,op,data->nlinfo,data->comm);
#endif
  for(i=0;i<n;++i) local_uncondense(u[i],data->local_cm);
}

/*--------------------------------------------------------------------------
   Setup
  --------------------------------------------------------------------------*/

gs_data *gs_data_setup(uint n, const ulong *label,
                       uint maxv, crystal_data *crystal)
{
  gs_data *data=tmalloc(gs_data,1);
  tuple_list nonzero, primary;
  const int nz_index=0, nz_size=1, nz_label=0;
  const int pr_nzindex=0, pr_index=1, pr_count=2, pr_size=3, pr_label=0;
#ifdef USE_MPI
  tuple_list shared;
  const int pr_proc=0;
  const int sh_dproc=0, sh_proc2=1, sh_index=2, sh_size=3, sh_label=0;
#else
  buffer buf;
#endif
#ifdef USE_MPI
  MPI_Comm_dup(crystal->comm,&data->comm);
#else
  buffer_init(&buf,1024);
#endif

  /* construct list of nonzeros: (index ^, label) */
  tuple_list_init_max(&nonzero,nz_size,1,0,n);
  {
    uint i; sint *nzi = nonzero.vi; slong *nzl = nonzero.vl;
    for(i=0;i<n;++i)
      if(label[i]!=0) 
        nzi[nz_index]=i,
        nzl[nz_label]=label[i],
        nzi+=nz_size, ++nzl, nonzero.n++;
  }

  /* sort nonzeros by label: (index ^2, label ^1) */
#ifndef USE_MPI
  tuple_list_sort(&nonzero,nz_size+nz_label,&buf);
#else
  tuple_list_sort(&nonzero,nz_size+nz_label,&crystal->all->buf);
#endif

  /* build list of unique labels w/ lowest associated index:
     (index in nonzero ^, primary (lowest) index in label, count, label) */
  tuple_list_init_max(&primary,pr_size,1,0,nonzero.n);
  {
    uint i;
    sint  *nzi=nonzero.vi, *pi=primary.vi;
    slong *nzl=nonzero.vl, *pl=primary.vl;
    sint last=-1;
    for(i=0;i<nonzero.n;++i,nzi+=nz_size,++nzl) {
      if(nzl[nz_label]==last) {
        ++pi[-pr_size+pr_count];
        continue;
      }
      last=nzl[nz_label];
      pi[pr_nzindex]=i;
      pi[pr_index]=nzi[nz_index];
      pl[pr_label]=nzl[nz_label];
      pi[pr_count]=1;
      pi+=pr_size, ++pl; primary.n++;
    }
  }

  /* calculate size of local condense map */
  {
    uint i, count=1; sint *pi=primary.vi;
    for(i=primary.n;i;--i,pi+=pr_size)
      if(pi[pr_count]>1) count+=pi[pr_count]+1;
    data->local_cm = tmalloc(sint,count);
  }

  /* sort unique labels by primary index:
     (nonzero index ^2, primary index ^1, count, label ^2) */
#ifndef USE_MPI
  tuple_list_sort(&primary,pr_index,&buf);
  buffer_free(&buf);
#else
  tuple_list_sort(&primary,pr_index,&crystal->all->buf);
#endif
  
  /* construct local condense map */
  {
    uint i, n; sint *pi=primary.vi;
    sint *cm = data->local_cm;
    for(i=primary.n;i;--i,pi+=pr_size) if((n=pi[pr_count])>1) {
      uint j; sint *nzi=nonzero.vi+nz_size*pi[pr_nzindex];
      for(j=n;j;--j,nzi+=nz_size) *cm++ = nzi[nz_index];
      *cm++ = -1;
    }
    *cm++ = -1;
  }
  tuple_list_free(&nonzero);
  
#ifndef USE_MPI
  tuple_list_free(&primary);
#else
  /* assign work proc by label modulo np */
  {
    uint i; sint *pi=primary.vi; slong *pl=primary.vl;
    for(i=primary.n;i;--i,pi+=pr_size,++pl)
      pi[pr_proc]=pl[pr_label]%crystal->num;
  }
  gs_transfer(1,&primary,pr_proc,crystal); /* transfer to work procs */
  /* primary: (source proc, index on src, useless, label) */
  /* sort by label */
  tuple_list_sort(&primary,pr_size+pr_label,&crystal->all->buf);
  /* add sentinel to primary list */
  if(primary.n==primary.max) tuple_list_grow(&primary);
  primary.vl[primary.n] = -1;
  /* construct shared list: (proc1, proc2, index1, label) */
  tuple_list_init_max(&shared,sh_size,1,0,primary.n);
  {
    sint *pi1=primary.vi, *si=shared.vi;
    slong lbl, *pl1=primary.vl, *sl=shared.vl;
    for(;(lbl=pl1[pr_label])!=-1;pi1+=pr_size,++pl1) {
      sint *pi2=pi1+pr_size; slong *pl2=pl1+1;
      for(;pl2[pr_label]==lbl;pi2+=pr_size,++pl2) {
        if(shared.n+2>shared.max)
          tuple_list_grow(&shared),
          si=shared.vi+shared.n*sh_size, sl=shared.vl+shared.n;
        si[sh_dproc] = pi1[pr_proc];
        si[sh_proc2] = pi2[pr_proc];
        si[sh_index] = pi1[pr_index];
        sl[sh_label] = lbl;
        si+=sh_size, ++sl, shared.n++;
        si[sh_dproc] = pi2[pr_proc];
        si[sh_proc2] = pi1[pr_proc];
        si[sh_index] = pi2[pr_index];
        sl[sh_label] = lbl;
        si+=sh_size, ++sl, shared.n++;
      }
    }
  }
  tuple_list_free(&primary);
  gs_transfer(1,&shared,sh_dproc,crystal); /* transfer to dest procs */
  /* shared list: (useless, proc2, index, label) */
  /* sort by label */
  tuple_list_sort(&shared,sh_size+sh_label,&crystal->all->buf);
  /* sort by partner proc */
  tuple_list_sort(&shared,sh_proc2,&crystal->all->buf);
  /* count partner procs */
  {
    uint i, count=0; sint proc=-1,*si=shared.vi;
    for(i=shared.n;i;--i,si+=sh_size)
      if(si[sh_proc2]!=proc) ++count, proc=si[sh_proc2];
    data->nlinfo = nlinfo_alloc(count,shared.n,maxv);
  }
  /* construct non-local info */
  {
    uint i; sint proc=-1,*si=shared.vi;
    uint *target  = data->nlinfo->target;
    uint *nshared = data->nlinfo->nshared;
    uint *sh_ind  = data->nlinfo->sh_ind;
    for(i=shared.n;i;--i,si+=sh_size) {
      if(si[sh_proc2]!=proc)
        proc=si[sh_proc2], *target++ = proc, *nshared++ = 0;
      ++nshared[-1], *sh_ind++=si[sh_index];
    }
  }
  tuple_list_free(&shared);
#endif
  return data;
}

void gs_data_free(gs_data *data)
{
  free(data->local_cm);
#ifdef USE_MPI
  nlinfo_free(data->nlinfo);
  MPI_Comm_free(&data->comm);
#endif
  free(data);
}

