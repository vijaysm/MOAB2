#ifndef TRANSFER_H
#define TRANSFER_H

#ifdef USE_MPI

#if !defined(TUPLE_LIST_H) || !defined(CRYSTAL_H)
#warning "transfer.h" requires "tuple_list.h" and "crystal.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*------------------------------------------------------------------------------
  
  Transfer
 
  Treats one integer (not long) member of the tuple list as a target proc;
  Sends out tuples accordingly, using the crystal router.
  Target proc member overwritten with source proc.
  
  dynamic: non-zero if the tuple list should grow to accomodate arrivals
  tl:      the tuple list
  pf:      which tuple member specifies target proc
  crystal: an initialized crystal router structure (cf. crystal.h)

  ----------------------------------------------------------------------------*/
void gs_transfer(int dynamic, tuple_list* tl,
                 unsigned pf, crystal_data *crystal);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif

#endif

