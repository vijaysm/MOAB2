/*------------------------------------------------------------------------------
  
  Tuple list definition and utilities
  
  Conceptually, a tuple list is a list of n records or tuples,
  each with mi integers, ml longs, and mr reals
  (these types are defined in "types.h" as sint, slong, real;
   it may be that sint==slong)
  
  There are three arrays, one for each type (vi,vl,vr),
  with records layed out contiguously within each array

  ----------------------------------------------------------------------------*/

#ifndef TUPLE_LIST_H
#define TUPLE_LIST_H

/* requires "errmem.h" and "types.h" */
#if !defined(ERRMEM_H) || !defined(TYPES_H) || !defined(MINMAX_H) || !defined(SORT_H)
#warning "tuple_list.h" requires "errmem.h" and "types.h" and  "minmax.h" and "sort.h"
#endif

typedef struct tuple_list {
  unsigned mi,ml,mul,mr;
  uint n, max;
  sint *vi; slong *vl; ulong *vul; real *vr;
} tuple_list;

/* storage layed out as: vi[max][mi], vl[max][ml], vr[max][mr]
   where "tuple" i is given by (vi[i][0:mi-1],vl[i][0:ml-1],vr[i][0:mr-1]).
   only the first n tuples are in use */

static void tuple_list_init_max(tuple_list *tl,
                                unsigned mi, unsigned ml, unsigned mul, 
                                unsigned mr, uint max)
{
  tl->n=0; tl->max=max;
  tl->mi=mi,tl->ml=ml,tl->mul=mul,tl->mr=mr;
  tl->vi=(max*mi ? tmalloc(sint, max*mi) : 0);
  tl->vl=(max*ml ? tmalloc(slong,max*ml) : 0);
  tl->vul=(max*mul ? tmalloc(ulong,max*mul) : 0);
  tl->vr=(max*mr ? tmalloc(real, max*mr) : 0);
}

static void tuple_list_free(tuple_list *tl) {
  free(tl->vi), free(tl->vl), free(tl->vul), free(tl->vr);
}

static void tuple_list_resize(tuple_list *tl, uint max)
{
  tl->max = max;
  if (tl->vi || (tl->max*tl->mi))
    tl->vi=trealloc(sint, tl->vi,tl->max*tl->mi);
  if (tl->vl || (tl->max*tl->ml))
    tl->vl=trealloc(slong,tl->vl,tl->max*tl->ml);
  if (tl->vul || (tl->max*tl->mul))
    tl->vul=trealloc(ulong,tl->vul,tl->max*tl->mul);
  if (tl->vr || (tl->max*tl->mr))
    tl->vr=trealloc(real, tl->vr,tl->max*tl->mr);
}

static void tuple_list_grow(tuple_list *tl)
{
  tuple_list_resize(tl,(tl->max ? tl->max+tl->max/2+1 : 2));
}

void tuple_list_permute(tuple_list *tl, uint *perm, void *work);
void tuple_list_sort(tuple_list *tl, unsigned key, buffer *buf);

#endif

