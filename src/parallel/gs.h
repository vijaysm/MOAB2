#ifndef GS_H
#define GS_H

/* requires "types.h", and, when MPI is defined, "crystal.h" */
#if !defined(TYPES_H) || ( defined(MPI) && !defined(CRYSTAL_H) )
#warning "gs.h" requires "types.h" and "crystal.h"
#endif

/*  typedef struct gs_data_ gs_data; */

#ifndef MPI
#  define crystal_data void
#endif

#ifdef USE_MPI
typedef struct {
  uint np;           /* number of processors to communicate with          */
  uint *target;      /* int target[np]: array of processor ids to comm w/ */
  uint *nshared;     /* nshared[i] = number of points shared w/ target[i] */
  uint *sh_ind;      /* list of shared point indices                      */
  slong *slabels;    /* list of signed long labels (not including gid)    */
  ulong *ulabels;    /* list of unsigned long labels                      */
  MPI_Request *reqs; /* pre-allocated for MPI calls                       */
  real *buf;         /* pre-allocated buffer to receive data              */
  uint maxv;         /* maximum vector size                               */
} nonlocal_info;
#endif

typedef struct {
  sint *local_cm; /* local condense map */
#ifdef USE_MPI
  nonlocal_info *nlinfo;
  MPI_Comm comm;
#endif
} moab_gs_data;

moab_gs_data *moab_gs_data_setup(uint n, const long *label, const ulong *ulabel,
                       uint maxv, const unsigned int nlabels,
                       const unsigned int nulabels,
                       crystal_data *crystal);

#ifndef MPI
#  undef crystal_data
#endif

void moab_gs_data_free(moab_gs_data *data);
void moab_gs_op(real *u, int op, const moab_gs_data *data);
void moab_gs_op_vec(real *u, uint n, int op, const moab_gs_data *data);
void moab_gs_op_many(real **u, uint n, int op, const moab_gs_data *data);

#define GS_OP_ADD 1
#define GS_OP_MUL 2
#define GS_OP_MIN 3
#define GS_OP_MAX 4
#define GS_OP_BPR 5

#endif

