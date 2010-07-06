/**
 * Copyright 2006 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */
#include "iRel.h"

#include "Lasso.hpp"
#include "AssocPairC.hpp"

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <math.h>
#include <stdio.h>

const bool debug = false;

// TAG_CHECK_SIZE is like CHECK_SIZE except it checks for and makes the allocated memory
// size a multiple of sizeof(void*), and the pointer is assumed to be type char*
#define TAG_CHECK_SIZE(array, allocated, size)  \
  if (NULL != array && 0 != array ## _allocated && array ## _allocated < (size)) {\
    iRel_processError(iBase_MEMORY_ALLOCATION_FAILED, \
          "Allocated array not large enough to hold returned contents.");\
    RETURN(iBase_MEMORY_ALLOCATION_FAILED);\
  }\
  if (NULL == array) {\
    allocated=(size); \
    if (allocated%sizeof(void*) != 0) allocated=((size)/sizeof(void*)+1)*sizeof(void*);\
    array = (char*)malloc(allocated); \
    if (NULL == array) {iRel_processError(iBase_MEMORY_ALLOCATION_FAILED, \
          "Couldn't allocate array.");return iBase_MEMORY_ALLOCATION_FAILED; }\
  }

#define CHECK_SIZE(type, array, allocated_size, size)  \
          if (NULL == *array || *allocated_size == 0) {\
            *array = (type *) malloc(sizeof(type) * size); \
            *allocated_size = size;} \
          else if (*allocated_size < size) { \
             std::cerr << "   Array passed in is non-zero but too short." << std::cerr; }

iBase_Error iRel_LAST_ERROR;

#define RETURN(a) {iRel_LAST_ERROR.error_type = *ierr = a; return;}
#define iRel_processError(a, b) {sprintf(iRel_LAST_ERROR.description, "%s", b); iRel_LAST_ERROR.error_type = *ierr = a;}

void
iRel_dtor(iRel_Instance instance, int *ierr) 
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  delete lasso;

  RETURN(iBase_SUCCESS);
}

void iRel_createRelation (
  iRel_Instance instance,
  iBase_Instance iface1,
  const int ent_or_set1,
  const int iface_type1,
  iBase_Instance iface2,
  const int ent_or_set2,
  const int iface_type2,
  iRel_RelationHandle *rel,
  int *ierr)
{
  // assume we want an AssocPairC
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  *rel = reinterpret_cast<iRel_RelationHandle>(
    new AssocPairC(
      iface1, static_cast<RelationType>(ent_or_set1),
      static_cast<IfaceType>(iface_type1),
      iface2, static_cast<RelationType>(ent_or_set2),
      static_cast<IfaceType>(iface_type2), lasso)
    );
  
  RETURN(iBase_SUCCESS);
}

void iRel_getRelationInfo (
  iRel_Instance instance,
  iRel_RelationHandle rel,
  iBase_Instance *iface1,
  int *ent_or_set1,
  int *iface_type1,
  iBase_Instance *iface2,
  int *ent_or_set2,
  int *iface_type2,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *assoc_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == assoc_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find relation pair.");
    RETURN(iBase_FAILURE);
  }

  *iface1      = assoc_pair->iface_instance(0);
  *ent_or_set1 = assoc_pair->ent_or_set(0);
  *iface_type1 = assoc_pair->iface_type(0);
  *iface2      = assoc_pair->iface_instance(1);
  *iface_type2 = assoc_pair->iface_type(1);
  *ent_or_set2 = assoc_pair->ent_or_set(1);

  RETURN(iBase_SUCCESS);
}

void iRel_destroyRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *assoc_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == assoc_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find relation pair.");
    RETURN(iBase_FAILURE);
  }

  int result = lasso->delete_pair(assoc_pair);

  RETURN(result);
}

void iRel_findRelations(
  iRel_Instance instance,
  iBase_Instance iface,
  iRel_RelationHandle **relations,
  int *relations_allocated,
  int *relations_size,
  int *ierr
  )
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  std::vector<AssocPair*> tmp_pairs;
  lasso->find_pairs(iface, tmp_pairs);
  if (tmp_pairs.empty()) {
    *relations_size = 0;
  }
  else {
    CHECK_SIZE(iRel_RelationHandle, relations, relations_allocated, 
               (int)tmp_pairs.size());
    *relations_size = tmp_pairs.size();
	for (size_t i=0; i<tmp_pairs.size(); ++i) {
      (*relations)[i] = reinterpret_cast<iRel_RelationHandle>(tmp_pairs[i]);
    }
  }

  RETURN(iBase_SUCCESS);
}

void iRel_setEntEntRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle ent1,
  iBase_EntityHandle ent2,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find relation pair.");
    RETURN(iBase_FAILURE);
  }

  int result = this_pair->set_assoc_tags(ent1, ent2);
  
    // xxx - need to check whether either is a set, and act accordingly!
  RETURN(result);
}

void iRel_setEntSetRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle ent1,
  iBase_EntitySetHandle set2,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find relation pair.");
    RETURN(iBase_FAILURE);
  }

  int result = this_pair->set_assoc_tags(ent1, set2);
  
  // xxx - need to check whether either is a set, and act accordingly!
  RETURN(result);
}

void iRel_setSetEntRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle set1,
  iBase_EntityHandle ent2,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find relation pair.");
    RETURN(iBase_FAILURE);
  }

  int result = this_pair->set_assoc_tags(set1, ent2);
  
  // xxx - need to check whether either is a set, and act accordingly!
  RETURN(result);
}

void iRel_setSetSetRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle set1,
  iBase_EntitySetHandle set2,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find relation pair.");
    RETURN(iBase_FAILURE);
  }

  int result = this_pair->set_assoc_tags(set1, set2);
  
  // xxx - need to check whether either is a set, and act accordingly!
  RETURN(result);
}

void iRel_setEntEntArrRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle ent1,
  int switch_order,
  iBase_EntityHandle *ent_array_2,
  int num_entities,
  int *ierr)
{
  int result = iBase_SUCCESS;
  for (int i = 0; i < num_entities; i++) {
    int tmp_result;
    if (switch_order) {
      iRel_setEntEntRelation(instance, rel,
                             ent_array_2[i], ent1, &tmp_result);
    }
    else {
      iRel_setEntEntRelation(instance, rel,
                             ent1, ent_array_2[i], &tmp_result);
    }
      
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}

void iRel_setEntSetArrRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle ent1,
  int switch_order,
  iBase_EntitySetHandle *set_array_2,
  int num_sets,
  int *ierr)
{
  int result = iBase_SUCCESS;
  for (int i = 0; i < num_sets; i++) {
    int tmp_result;
    if (switch_order) {
      iRel_setSetEntRelation(instance, rel,
                             set_array_2[i], ent1, &tmp_result);
    }
    else {
      iRel_setEntSetRelation(instance, rel,
                             ent1, set_array_2[i], &tmp_result);
    }
      
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}

void iRel_setSetEntArrRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle set1,
  int switch_order,
  iBase_EntityHandle *ent_array_2,
  int num_entities,
  int *ierr)
{
  int result = iBase_SUCCESS;
  for (int i = 0; i < num_entities; i++) {
    int tmp_result;
    if (switch_order) {
      iRel_setEntSetRelation(instance, rel,
                             ent_array_2[i], set1, &tmp_result);
    }
    else {
      iRel_setSetEntRelation(instance, rel,
                             set1, ent_array_2[i], &tmp_result);
    }
      
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}

void iRel_setSetSetArrRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle set1,
  int switch_order,
  iBase_EntitySetHandle *set_array_2,
  int num_sets,
  int *ierr)
{
  int result = iBase_SUCCESS;
  for (int i = 0; i < num_sets; i++) {
    int tmp_result;
    if (switch_order) {
      iRel_setSetSetRelation(instance, rel,
                             set_array_2[i], set1, &tmp_result);
    }
    else {
      iRel_setSetSetRelation(instance, rel,
                             set1, set_array_2[i], &tmp_result);
    }
      
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}

void iRel_setEntArrEntArrRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle *ent_array_1,
  int num_entities1,
  iBase_EntityHandle *ent_array_2,
  int num_entities2,
  int *ierr)
{
  if (num_entities1 != num_entities2) {
    iRel_processError(iBase_NOT_SUPPORTED, "setEntArrEntArrRelation doesn't "
                      "support different #'s of entities.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  int result = iBase_SUCCESS;
  for (int i = 0; i < num_entities1; i++) {
    int tmp_result;
    iRel_setEntEntRelation(instance, rel,
                           ent_array_1[i], ent_array_2[i], &tmp_result);
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}

void iRel_setEntArrSetArrRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle *ent_array_1,
  int num_entities1,
  iBase_EntitySetHandle *set_array_2,
  int num_sets2,
  int *ierr)
{
  if (num_entities1 != num_sets2) {
    iRel_processError(iBase_NOT_SUPPORTED, "setEntArrSetArrRelation doesn't "
                      "support different #'s of entities.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  int result = iBase_SUCCESS;
  for (int i = 0; i < num_entities1; i++) {
    int tmp_result;
    iRel_setEntSetRelation(instance, rel,
                           ent_array_1[i], set_array_2[i], &tmp_result);
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}

void iRel_setSetArrEntArrRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle *set_array_1,
  int num_sets1,
  iBase_EntityHandle *ent_array_2,
  int num_entities2,
  int *ierr)
{
  if (num_sets1 != num_entities2) {
    iRel_processError(iBase_NOT_SUPPORTED, "setSetArrEntArrRelation doesn't "
                      "support different #'s of entities.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  int result = iBase_SUCCESS;
  for (int i = 0; i < num_sets1; i++) {
    int tmp_result;
    iRel_setSetEntRelation(instance, rel,
                           set_array_1[i], ent_array_2[i], &tmp_result);
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}

void iRel_setSetArrSetArrRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle *set_array_1,
  int num_sets1,
  iBase_EntitySetHandle *set_array_2,
  int num_sets2,
  int *ierr)
{
  if (num_sets1 != num_sets2) {
    iRel_processError(iBase_NOT_SUPPORTED, "setSetArrSetArrRelation doesn't "
                      "support different #'s of entities.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  int result = iBase_SUCCESS;
  for (int i = 0; i < num_sets1; i++) {
    int tmp_result;
    iRel_setSetSetRelation(instance, rel,
                           set_array_1[i], set_array_2[i], &tmp_result);
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}

void iRel_getEntEntRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle ent1,
  int switch_order,
  iBase_EntityHandle *ent2,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find relation pair.");
    RETURN(iBase_FAILURE);
  }
  
  int iface_no = (switch_order ? 1 : 0);
  
  if (this_pair->ent_or_set(!iface_no) > 0) { // iface2 is sets
    iRel_processError(iBase_INVALID_ENTITY_HANDLE,
                      "Expected EntitySet, got Entity");
    RETURN(iBase_INVALID_ENTITY_HANDLE);
  }
  
  int result = this_pair->get_assoc_tags(iface_no, &ent1, 1, ent2);
  RETURN(result);
}

void iRel_getEntSetRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle ent1,
  int switch_order,
  iBase_EntitySetHandle *ent2,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find relation pair.");
    RETURN(iBase_FAILURE);
  }
  
  int iface_no = (switch_order ? 1 : 0);
  
  if (this_pair->ent_or_set(!iface_no) == 0) { // iface2 is not sets
    iRel_processError(iBase_INVALID_ENTITY_HANDLE,
                      "Expected Entity, got EntitySet");
    RETURN(iBase_INVALID_ENTITY_HANDLE);
  }
  
  int result = this_pair->get_assoc_tags(iface_no, &ent1, 1, ent2);
  
  RETURN(result);
}

void iRel_getSetEntRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle ent1,
  int switch_order,
  iBase_EntityHandle *ent2,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find relation pair.");
    RETURN(iBase_FAILURE);
  }
  
  int iface_no = (switch_order ? 1 : 0);
  
  if (this_pair->ent_or_set(!iface_no) > 0) { // iface2 is sets
    iRel_processError(iBase_INVALID_ENTITY_HANDLE,
                      "Expected EntitySet, got Entity");
    RETURN(iBase_INVALID_ENTITY_HANDLE);
  }
  
  int result = this_pair->get_assoc_tags(iface_no, &ent1, 1, ent2);
  
  RETURN(result);
}

void iRel_getSetSetRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle ent1,
  int switch_order,
  iBase_EntitySetHandle *ent2,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find relation pair.");
    RETURN(iBase_FAILURE);
  }
  
  int iface_no = (switch_order ? 1 : 0);
  
  if (this_pair->ent_or_set(!iface_no) == 0) { // iface2 is not sets
    iRel_processError(iBase_INVALID_ENTITY_HANDLE,
                      "Expected Entity, got EntitySet");
    RETURN(iBase_INVALID_ENTITY_HANDLE);
  }
  
  int result = this_pair->get_assoc_tags(iface_no, &ent1, 1, ent2);
  
  RETURN(result);
}

void iRel_getEntEntArrRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle ent1,
  int switch_order,
  iBase_EntityHandle **ent_array_2,
  int *ent_array_2_allocated,
  int *ent_array_2_size,
  int *ierr)
{
  int result;

  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find relation pair.");
    RETURN(iBase_FAILURE);
  }
  int iface_no = (switch_order ? 1 : 0);
  
  if (this_pair->ent_or_set(!iface_no) > 0) { // iface2 is sets
    iBase_EntitySetHandle ent_set2 = 0;

    iRel_getEntSetRelation(instance, rel,
                           ent1, switch_order,
                           &ent_set2, &result);
    if (iBase_SUCCESS != result) RETURN(result);

    result = this_pair->get_entities(!iface_no, -1, ent_set2, ent_array_2,
                                     ent_array_2_allocated, ent_array_2_size);
  }
  else {
    // otherwise put the set into the ent list
    CHECK_SIZE(iBase_EntityHandle, ent_array_2, ent_array_2_allocated, 1);
    *ent_array_2_size = 1;
    iRel_getEntEntRelation(instance, rel, ent1, switch_order, 
                           *ent_array_2, &result);
  }
  
  RETURN(result);
}

void iRel_getSetEntArrRelation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle set1,
  int switch_order,
  iBase_EntityHandle **ent_array_2,
  int *ent_array_2_allocated,
  int *ent_array_2_size,
  int *ierr)
{
  int result;

  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find relation pair.");
    RETURN(iBase_FAILURE);
  }
  int iface_no = (switch_order ? 1 : 0);
  
  if (this_pair->ent_or_set(!iface_no) > 0) { // iface2 is sets
    iBase_EntitySetHandle set2 = 0;

    iRel_getSetSetRelation(instance, rel,
                           set1, switch_order, &set2, &result);
    if (iBase_SUCCESS != result) RETURN(result);

    result = this_pair->get_entities(!iface_no, -1, set2, ent_array_2,
                                     ent_array_2_allocated, ent_array_2_size);
  }
  else {
    // otherwise put the set into the ent list
    CHECK_SIZE(iBase_EntityHandle, ent_array_2, ent_array_2_allocated, 1);
    *ent_array_2_size = 1;
    iRel_getSetEntRelation(instance, rel, set1, switch_order, 
                           *ent_array_2, &result);
  }
  
  RETURN(result);
}

void iRel_getEntArrEntArrRelation(
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle *ent_array_1,
  int ent_array_1_size,
  int switch_order,
  iBase_EntityHandle **ent_array_2,
  int *ent_array_2_allocated,
  int *ent_array_2_size,
  int **offset,
  int *offset_allocated,
  int *offset_size,
  int *ierr)
{
  std::vector<iBase_EntityHandle> tmp_ents;
  int result = iBase_SUCCESS;
  CHECK_SIZE(int, offset, offset_allocated, (ent_array_1_size+1));
  *offset_size = ent_array_1_size+1;

  // get assocpair so we can check whether ents are related to
  // entities or sets or both
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  bool switched = false;
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find relation pair.");
    RETURN(iBase_FAILURE);
  }

  for (int i = 0; i < ent_array_1_size; i++) {
    (*offset)[i] = tmp_ents.size();
    int tmp_result;
    iBase_EntityHandle *tmp_array = NULL;
    int tmp_array_size, tmp_array_allocated;
    iRel_getEntEntArrRelation(instance, rel,
                              ent_array_1[i], 
                              switch_order,
                              &tmp_array, &tmp_array_allocated,
                              &tmp_array_size, &tmp_result);
    if (iBase_SUCCESS == result && NULL != tmp_array) {
      std::copy(tmp_array, tmp_array+tmp_array_size,
                std::back_inserter(tmp_ents));
      free(tmp_array);
    }
    
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }

  (*offset)[ent_array_1_size] = tmp_ents.size();
  CHECK_SIZE(iBase_EntityHandle, ent_array_2, ent_array_2_allocated,
             (int)tmp_ents.size());
  *ent_array_2_size = tmp_ents.size();
  std::copy(tmp_ents.begin(), tmp_ents.end(), *ent_array_2);
  
  RETURN(result);
}

void iRel_getEntArrSetArrRelation(
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle *ent_array_1,
  int ent_array_1_size,
  int switch_order,
  iBase_EntitySetHandle **set_array_2,
  int *set_array_2_allocated,
  int *set_array_2_size,
  int *ierr)
{
  std::vector<iBase_EntityHandle> tmp_ents;
  int tmp_result, result = iBase_SUCCESS;
  CHECK_SIZE(iBase_EntitySetHandle, set_array_2, set_array_2_allocated,
             ent_array_1_size);
  *set_array_2_size = ent_array_1_size;

  // get assocpair so we can check whether ents are related to
  // entities or sets or both
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  bool switched = false;
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find relation pair.");
    RETURN(iBase_FAILURE);
  }

  for (int i = 0; i < ent_array_1_size; i++) {
    iRel_getEntSetRelation( instance, rel, ent_array_1[i],
                            switch_order, (*set_array_2)+i,
                            &tmp_result );
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}

void iRel_getSetArrEntArrRelation(
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle *set_array_1,
  int set_array_1_size,
  int switch_order,
  iBase_EntityHandle **ent_array_2,
  int *ent_array_2_allocated,
  int *ent_array_2_size,
  int **offset,
  int *offset_allocated,
  int *offset_size,
  int *ierr)
{
  std::vector<iBase_EntityHandle> tmp_ents;
  int result = iBase_SUCCESS;
  CHECK_SIZE(int, offset, offset_allocated, (set_array_1_size+1));
  *offset_size = set_array_1_size+1;

  // get assocpair so we can check whether ents are related to
  // entities or sets or both
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  bool switched = false;
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find relation pair.");
    RETURN(iBase_FAILURE);
  }

  for (int i = 0; i < set_array_1_size; i++) {
    (*offset)[i] = tmp_ents.size();
    int tmp_result;
    iBase_EntityHandle *tmp_array = NULL;
    int tmp_array_size, tmp_array_allocated;
    iRel_getSetEntArrRelation(instance, rel,
                              set_array_1[i], 
                              switch_order,
                              &tmp_array, &tmp_array_allocated,
                              &tmp_array_size, &tmp_result);
    if (iBase_SUCCESS == result && NULL != tmp_array) {
      std::copy(tmp_array, tmp_array+tmp_array_size,
                std::back_inserter(tmp_ents));
      free(tmp_array);
    }
    
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }

  (*offset)[set_array_1_size] = tmp_ents.size();
  CHECK_SIZE(iBase_EntityHandle, ent_array_2, ent_array_2_allocated,
             (int)tmp_ents.size());
  *ent_array_2_size = tmp_ents.size();
  std::copy(tmp_ents.begin(), tmp_ents.end(), *ent_array_2);
  
  RETURN(result);
}

void iRel_getSetArrSetArrRelation(
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle *set_array_1,
  int set_array_1_size,
  int switch_order,
  iBase_EntitySetHandle **set_array_2,
  int *set_array_2_allocated,
  int *set_array_2_size,
  int *ierr)
{
  std::vector<iBase_EntityHandle> tmp_ents;
  int tmp_result, result = iBase_SUCCESS;
  CHECK_SIZE(iBase_EntitySetHandle, set_array_2, set_array_2_allocated,
             set_array_1_size );
  *set_array_2_size = set_array_1_size;

  // get assocpair so we can check whether ents are related to
  // entities or sets or both
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find relation pair.");
    RETURN(iBase_FAILURE);
  }

  for (int i = 0; i < set_array_1_size; i++) {
    iRel_getSetSetRelation( instance, rel, set_array_1[i],
                            switch_order, (*set_array_2)+i,
                            &tmp_result );
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}

static int
get_gids_and_dims(AssocPair *assoc_pair,
                  int iface_no,
                  iBase_EntityHandle *ents,
                  int ents_size,
                  int is_set,
                  std::vector<int> &ents_gids,
                  std::vector<int> &ents_dims)
{
  int result;
  iBase_EntitySetHandle *sets = reinterpret_cast<iBase_EntitySetHandle*>(ents);

  ents_gids.resize(ents_size);
  if (is_set) {
    result = assoc_pair->get_gid_tags(iface_no, sets, ents_size, &ents_gids[0]);
  }
  else {
    result = assoc_pair->get_gid_tags(iface_no, ents, ents_size, &ents_gids[0]);
  }

  if (iBase_SUCCESS != result && iBase_TAG_NOT_FOUND != result) {
    return result;
  }

  ents_dims.resize(ents_size, -1);
  if (is_set) {
    result = assoc_pair->get_dim_tags(iface_no, sets, ents_size, &ents_dims[0]);
  }
  else {
    int *ents_dims_ptr = &ents_dims[0];
    int ents_dims_alloc = ents_dims.size(), ents_dims_size;
    result = assoc_pair->get_ents_dims(iface_no, ents, ents_size,
                                       &ents_dims_ptr, &ents_dims_alloc,
                                       &ents_dims_size);
  }

  if (iBase_SUCCESS != result && iBase_TAG_NOT_FOUND != result) {
    return result;
  }
}

static void
iRel_inferArrArrRelations(Lasso *lasso,
                          AssocPair *assoc_pair,
                          iBase_EntityHandle *ents1,
                          const int ents1_size,
                          int is_set1,
                          iBase_EntityHandle *ents2,
                          const int ents2_size,
                          int is_set2,
                          int *ierr)
{
  // get dimension and global id tags for sets
  int tmp_index = 0;
  int ents_index[2] = {tmp_index, 1-tmp_index};
  iBase_EntityHandle *ents[2] = {ents1, ents2};
  iBase_EntitySetHandle *sets[2] = {
    reinterpret_cast<iBase_EntitySetHandle*>(ents1),
    reinterpret_cast<iBase_EntitySetHandle*>(ents2) };
  int ents_size[2] = {ents1_size, ents2_size};
  
  int result;

  std::vector<int> ents_gids, ents_dims;
  std::map<const int, iBase_EntityHandle> ents_gid_map[4];

  get_gids_and_dims(assoc_pair, 0, ents1, ents1_size, is_set1, ents_gids,
                    ents_dims);
  for (int i = 0; i < ents1_size; i++) {
    int dim = ents_dims[i];
    if (0 <= dim && 3 >= dim)
      ents_gid_map[dim][ents_gids[i]] = ents1[i];
  }

  get_gids_and_dims(assoc_pair, 1, ents2, ents2_size, is_set2, ents_gids,
                    ents_dims);
  for (int i = 0; i < ents2_size; i++) {
    int dim = ents_dims[i];

    // only check entities for which the dimension entry is in a reasonable
    // range
    if (0 > dim || 3 < dim) continue;

    // there's a match if there's an entity with that dimension with matching id
    std::map<const int, iBase_EntityHandle>::iterator iter =
      ents_gid_map[dim].find(ents_gids[i]);

      // if it matches, set the relation tags for those entities
    if (iter != ents_gid_map[dim].end()) {
      if (!is_set1 && !is_set2) {
        result = assoc_pair->set_assoc_tags(
          (*iter).second,
          ents2[i]);
      }
      else if (is_set1 && !is_set2) {
        result = assoc_pair->set_assoc_tags(
          (iBase_EntitySetHandle)(*iter).second,
          ents2[i]);
      }
      else if (!is_set1 && is_set2) {
        result = assoc_pair->set_assoc_tags(
          (*iter).second,
          (iBase_EntitySetHandle)ents2[i]);
      }
      else { // is_set1 && is_set2
        result = assoc_pair->set_assoc_tags(
          (iBase_EntitySetHandle)(*iter).second,
          (iBase_EntitySetHandle)ents2[i]);
      }

      if (iBase_SUCCESS != result) {
        RETURN(result);
      }
    }
  }

  RETURN(iBase_SUCCESS);
}

void
iRel_inferEntArrEntArrRelations(Lasso *lasso,
                                AssocPair *assoc_pair,
                                iBase_EntityHandle *ents1,
                                const int ents1_size,
                                iBase_EntityHandle *ents2,
                                const int ents2_size,
                                int *ierr)
{
  iRel_inferArrArrRelations(lasso, assoc_pair, ents1, ents1_size, 0,
                            ents2, ents2_size, 0, ierr);
}

void
iRel_inferEntArrSetArrRelations(Lasso *lasso,
                                AssocPair *assoc_pair,
                                iBase_EntityHandle *ents1,
                                const int ents1_size,
                                iBase_EntitySetHandle *ents2,
                                const int ents2_size,
                                int *ierr)
{
  iRel_inferArrArrRelations(lasso, assoc_pair, ents1, ents1_size, 0,
                            (iBase_EntityHandle*)ents2, ents2_size, 1, ierr);
}

void
iRel_inferSetArrEntArrRelations(Lasso *lasso,
                                AssocPair *assoc_pair,
                                iBase_EntitySetHandle *ents1,
                                const int ents1_size,
                                iBase_EntityHandle *ents2,
                                const int ents2_size,
                                int *ierr)
{
  iRel_inferArrArrRelations(lasso, assoc_pair, (iBase_EntityHandle*)ents1,
                            ents1_size, 1, ents2, ents2_size, 0, ierr);
}

void
iRel_inferSetArrSetArrRelations(Lasso *lasso,
                                AssocPair *assoc_pair,
                                iBase_EntitySetHandle *ents1,
                                const int ents1_size,
                                int is_set1,
                                iBase_EntitySetHandle *ents2,
                                const int ents2_size,
                                int is_set2,
                                int *ierr)

{
  iRel_inferArrArrRelations(lasso, assoc_pair, (iBase_EntityHandle*)ents1,
                            ents1_size, 1, (iBase_EntityHandle*)ents2, 
                            ents2_size, 1, ierr);
}

void iRel_inferAllRelations (
  iRel_Instance instance,
  iRel_RelationHandle rel,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find relation pair.");
    RETURN(iBase_FAILURE);
  }

  // get all entities in those interfaces
  int result;

  iBase_EntityHandle *ents1 = NULL;
  int ents1_alloc = 0, ents1_size;
  if (this_pair->ent_or_set(0) > 0)
    result = this_pair->get_all_sets(0, (iBase_EntitySetHandle**)&ents1,
                                     &ents1_alloc, &ents1_size);
  else
    result = this_pair->get_all_entities(0, -1, &ents1, &ents1_alloc,
                                         &ents1_size);
  if (iBase_SUCCESS != result) {
    RETURN(result);
  }
  
  iBase_EntityHandle *ents2 = NULL;
  int ents2_alloc = 0, ents2_size;
  if (this_pair->ent_or_set(1) > 0)
    result = this_pair->get_all_sets(1, (iBase_EntitySetHandle**)&ents2,
                                     &ents2_alloc, &ents2_size);
  else
  result = this_pair->get_all_entities(1, -1, &ents2, &ents2_alloc,
                                       &ents2_size);
  if (iBase_SUCCESS != result) {
    RETURN(result);
  }
  
  iRel_inferArrArrRelations(lasso, this_pair, 
                            ents1, ents1_size, this_pair->ent_or_set(0),
                            ents2, ents2_size, this_pair->ent_or_set(1),
                            &result);

  free(ents1); free(ents2);
  RETURN(result);
}

void iRel_inferEntRelations (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle entity,
  int iface_no,
  int *ierr)
{
  iRel_inferEntArrRelations(instance, rel, &entity, 1, iface_no, ierr);
}

void iRel_inferSetRelations (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle entity,
  int iface_no,
  int *ierr)
{
  iRel_inferSetArrRelations(instance, rel, &entity, 1, iface_no, ierr);
}

static void iRel_inferArrRelations (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle *entities,
  int entities_size,
  int is_set,
  int iface_no,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find relation pair.");
    RETURN(iBase_FAILURE);
  }

  if (0 > iface_no || 1 < iface_no) {
    iRel_processError(iBase_INVALID_ARGUMENT,
                      "Interface number must be 0 or 1");
    RETURN(iBase_INVALID_ARGUMENT);
  }
  else if (( is_set && this_pair->ent_or_set(iface_no) == 0) ||
           (!is_set && this_pair->ent_or_set(iface_no)  > 0)) {
    iRel_processError(iBase_INVALID_ARGUMENT, "is_set must match entOrSet in "
                      "call to inferArrRelations");
    RETURN(iBase_INVALID_ARGUMENT);
  }
  
  // get all entities in iface2
  int result;

  iBase_EntityHandle *ents1 = entities;
  int ents1_size = entities_size;
  iBase_EntityHandle *ents2 = NULL;
  int ents2_alloc = 0, ents2_size;
  if (this_pair->ent_or_set(1-iface_no) > 0)
    result = this_pair->get_all_sets(!iface_no, 
                                     (iBase_EntitySetHandle**)&ents2, 
                                     &ents2_alloc, &ents2_size);
  else
    result = this_pair->get_all_entities(!iface_no, -1, 
                                         &ents2, &ents2_alloc, &ents2_size);
  if (iBase_SUCCESS != result) {
    RETURN(result);
  }

  // switch so that entity lists always go into inferArrArrRelations in
  // forward order wrt assoc_pair
  if (1 == iface_no) {
    std::swap(ents1, ents2);
    std::swap(ents1_size, ents2_size);
  }
  
  iRel_inferArrArrRelations(lasso, this_pair,
                            ents1, ents1_size, 
                            this_pair->ent_or_set(0),
                            ents2, ents2_size, 
                            this_pair->ent_or_set(1), &result);

  free(1 == iface_no ? ents1 : ents2);
  
  RETURN(result);
}

void iRel_inferEntArrRelations (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle *entities,
  int entities_size,
  int iface_no,
  int *ierr)
{
  iRel_inferArrRelations(instance, rel, entities, entities_size, 0, iface_no,
                         ierr);
}

void iRel_inferSetArrRelations (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle *entities,
  int entities_size,
  int iface_no,
  int *ierr)
{
  iRel_inferArrRelations(instance, rel, (iBase_EntityHandle*)entities,
                         entities_size, 1, iface_no, ierr);
}

void iRel_newRel(/* in */ const char * /* options */,
                 iRel_Instance *instance,
                 int *ierr,
                 const int options_len) 
{
  if (0 != options_len) {
    iRel_processError(iBase_NOT_SUPPORTED, "No options for iRel factory have "
                      "been implemented.");
    *instance = NULL;
  }
  
  *instance = new Lasso();
  RETURN(iBase_SUCCESS);
}

Lasso::~Lasso() 
{
  for (std::vector<AssocPair*>::iterator vit = assocPairs.begin();
       vit != assocPairs.end(); vit++) 
    delete *vit;
}

//! find a pair equivalent to these ifaces, passed as pointer to
//! SIDL interface or interface instance
AssocPair *Lasso::find_pair(void *iface0, void *iface1, bool *switched) 
{
  for (std::vector<AssocPair*>::iterator vit = assocPairs.begin();
       vit != assocPairs.end(); vit++) {
    if ((*vit)->equivalent(iface0, iface1, switched)) return *vit;
  }
  
  return NULL;
}

//! find a pair with the right types
AssocPair *Lasso::find_pair(IfaceType type1, IfaceType type2, bool *switched) 
{
  for (std::vector<AssocPair*>::iterator vit = assocPairs.begin();
       vit != assocPairs.end(); vit++) {
    if ((*vit)->equivalent(type1, type2, switched)) return *vit;
  }
  
  return NULL;
}

void Lasso::find_pairs(void *iface, std::vector<AssocPair*> &iface_pairs)
{
  for (std::vector<AssocPair*>::iterator vit = assocPairs.begin();
       vit != assocPairs.end(); vit++) {
    if ((*vit)->contains(iface)) iface_pairs.push_back(*vit);
  }
}

int Lasso::delete_pair(AssocPair *this_pair) 
{
  if (std::find(assocPairs.begin(), assocPairs.end(),
                this_pair) == assocPairs.end())
    return iBase_FAILURE;
  
  // remove the pair from the list
  assocPairs.erase(std::remove(assocPairs.begin(),
                               assocPairs.end(), this_pair), 
                   assocPairs.end());
  
  // then delete it
  delete this_pair;

  return iBase_SUCCESS;
}
