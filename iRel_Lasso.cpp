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
#define iRel_processError(a, b) {sprintf(iRel_LAST_ERROR.description, b); iRel_LAST_ERROR.error_type = *ierr = a;}

void
iRel_dtor(iRel_Instance instance, int *ierr) 
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  delete lasso;

  RETURN(iBase_SUCCESS);
}

void iRel_createAssociation (
  iRel_Instance instance,
  iBase_Instance iface1,
  const int ent_or_set1,
  const int type1,
  iBase_Instance iface2,
  const int ent_or_set2,
  const int type2,
  iRel_RelationHandle *rel,
  int *ierr)
{
    // error check - can't have both interfaces be 'both'-type
    // associations
  if (ent_or_set1 == 2 && ent_or_set2 == 2) {
    iRel_processError(iBase_INVALID_ARGUMENT, 
                      "Can't have both parts of association be 'both'-type associations.");
    RETURN(iBase_INVALID_ARGUMENT);
  }
    
    // assume it's an AssocPairC
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  *rel = (iRel_RelationHandle) new AssocPairC(iface1, ent_or_set1, static_cast<IfaceType>(type1),
                                              iface2, ent_or_set2, static_cast<IfaceType>(type2),
                                              lasso);
  
  RETURN(iBase_SUCCESS);
}


void iRel_destroyAssociation (
  iRel_Instance instance,
  iRel_RelationHandle rel,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *assoc_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == assoc_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find association pair.");
    RETURN(iBase_FAILURE);
  }

  int result = lasso->delete_pair(assoc_pair);

  RETURN(result);
}

void iRel_getAssociatedInterfaces (
  iRel_Instance instance,
  iBase_Instance iface,
  iBase_Instance **interfaces,
  int *interfaces_allocated,
  int *interfaces_size,
  int *ierr
  )
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  std::vector<AssocPair*> tmp_ifaces;
  lasso->find_pairs(iface, tmp_ifaces);
  if (tmp_ifaces.empty()) {
    *interfaces_size = 0;
    iRel_LAST_ERROR.error_type = iBase_SUCCESS;
  }
  else {
    CHECK_SIZE(iBase_Instance, interfaces, interfaces_allocated, 
               (int)tmp_ifaces.size());
    *interfaces_size = tmp_ifaces.size();
    std::copy(tmp_ifaces.begin(), tmp_ifaces.end(), *interfaces);
  }
  RETURN(iRel_LAST_ERROR.error_type);
}

void iRel_setEntEntAssociation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle ent1,
  iBase_EntityHandle ent2,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find association pair.");
    RETURN(iBase_FAILURE);
  }

  int result = this_pair->set_assoc_tags(ent1, ent2);
  
    // xxx - need to check whether either is a set, and act accordingly!
  RETURN(result);
}
void iRel_setEntSetAssociation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle ent1,
  iBase_EntitySetHandle ent2,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find association pair.");
    RETURN(iBase_FAILURE);
  }

  int result = this_pair->set_assoc_tags(ent1, ent2);
  
    // xxx - need to check whether either is a set, and act accordingly!
  RETURN(result);
}
void iRel_setSetEntAssociation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle ent1,
  iBase_EntityHandle ent2,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find association pair.");
    RETURN(iBase_FAILURE);
  }

  int result = this_pair->set_assoc_tags(ent1, ent2);
  
    // xxx - need to check whether either is a set, and act accordingly!
  RETURN(result);
}
void iRel_setSetSetAssociation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle ent1,
  iBase_EntitySetHandle ent2,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find association pair.");
    RETURN(iBase_FAILURE);
  }

  int result = this_pair->set_assoc_tags(ent1, ent2);
  
    // xxx - need to check whether either is a set, and act accordingly!
  RETURN(result);
}

void iRel_setEntEntArrAssociation (
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
      iRel_setEntEntAssociation(instance, rel,
                                ent_array_2[i], ent1, &tmp_result);
    }
    else {
      iRel_setEntEntAssociation(instance, rel,
                                ent1, ent_array_2[i], &tmp_result);
    }
      
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}
void iRel_setEntSetArrAssociation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle ent1,
  int switch_order,
  iBase_EntitySetHandle *ent_array_2,
  int num_entities,
  int *ierr)
{
  int result = iBase_SUCCESS;
  for (int i = 0; i < num_entities; i++) {
    int tmp_result;
    if (switch_order) {
      iRel_setSetEntAssociation(instance, rel,
                                ent_array_2[i], ent1, &tmp_result);
    }
    else {
      iRel_setEntSetAssociation(instance, rel,
                                ent1, ent_array_2[i], &tmp_result);
    }
      
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}
void iRel_setSetEntArrAssociation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle ent1,
  int switch_order,
  iBase_EntityHandle *ent_array_2,
  int num_entities,
  int *ierr)
{
  int result = iBase_SUCCESS;
  for (int i = 0; i < num_entities; i++) {
    int tmp_result;
    if (switch_order) {
      iRel_setEntSetAssociation(instance, rel,
                                ent_array_2[i],  
                                ent1,  
                                &tmp_result);
    }
    else {
      iRel_setSetEntAssociation(instance, rel,
                                ent1,  
                                ent_array_2[i],  
                                &tmp_result);
    }
      
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}
void iRel_setSetSetArrAssociation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle ent1,
  int switch_order,
  iBase_EntitySetHandle *ent_array_2,
  int num_entities,
  int *ierr)
{
  int result = iBase_SUCCESS;
  for (int i = 0; i < num_entities; i++) {
    int tmp_result;
    if (switch_order) {
      iRel_setSetSetAssociation(instance, rel,
                                ent_array_2[i],  
                                ent1,  
                                &tmp_result);
    }
    else {
      iRel_setSetSetAssociation(instance, rel,
                                ent1,  
                                ent_array_2[i],  
                                &tmp_result);
    }
      
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
//
}

void iRel_setEntArrEntArrAssociation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle *ent_array_1,
  int num_ent1,
  iBase_EntityHandle *ent_array_2,
  int num_ent2,
  int *ierr)
{
  if (num_ent1 != num_ent2) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "setArrAssocation doesn't support different #'s of entities.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  int result = iBase_SUCCESS;
  for (int i = 0; i < num_ent1; i++) {
    int tmp_result;
    iRel_setEntEntAssociation(instance, rel,
                              ent_array_1[i],  
                              ent_array_2[i],  &tmp_result);
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}
void iRel_setEntArrSetArrAssociation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle *ent_array_1,
  int num_ent1,
  iBase_EntitySetHandle *ent_array_2,
  int num_ent2,
  int *ierr)
{
  if (num_ent1 != num_ent2) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "setArrAssocation doesn't support different #'s of entities.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  int result = iBase_SUCCESS;
  for (int i = 0; i < num_ent1; i++) {
    int tmp_result;
    iRel_setEntSetAssociation(instance, rel,
                              ent_array_1[i],  
                              ent_array_2[i],  &tmp_result);
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}
void iRel_setSetArrEntArrAssociation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle *ent_array_1,
  int num_ent1,
  iBase_EntityHandle *ent_array_2,
  int num_ent2,
  int *ierr)
{
  if (num_ent1 != num_ent2) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "setArrAssocation doesn't support different #'s of entities.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  int result = iBase_SUCCESS;
  for (int i = 0; i < num_ent1; i++) {
    int tmp_result;
    iRel_setSetEntAssociation(instance, rel,
                              ent_array_1[i],  
                              ent_array_2[i],  &tmp_result);
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}
void iRel_setSetArrSetArrAssociation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle *ent_array_1,
  int num_ent1,
  iBase_EntitySetHandle *ent_array_2,
  int num_ent2,
  int *ierr)
{
  if (num_ent1 != num_ent2) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "setArrAssocation doesn't support different #'s of entities.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  int result = iBase_SUCCESS;
  for (int i = 0; i < num_ent1; i++) {
    int tmp_result;
    iRel_setSetSetAssociation(instance, rel,
                              ent_array_1[i],  
                              ent_array_2[i], &tmp_result);
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}

void iRel_getEntEntAssociation (
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
    iRel_processError(iBase_FAILURE, "Didn't find association pair.");
    RETURN(iBase_FAILURE);
  }
  
  int iface_no = (switch_order ? 1 : 0);
  
  if (this_pair->entOrSet[1-iface_no] > 0) { // iface2 is sets
    iRel_processError(iBase_INVALID_ENTITY_HANDLE, "Expected EntitySet, got Entity");
    RETURN(iBase_INVALID_ENTITY_HANDLE);
  }
  
  int result = this_pair->get_eh_tags(iface_no, 
                           &ent1, 1, this_pair->assocTags[iface_no], ent2);
  
  RETURN(result);
}
void iRel_getEntSetAssociation (
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
    iRel_processError(iBase_FAILURE, "Didn't find association pair.");
    RETURN(iBase_FAILURE);
  }
  
  int iface_no = (switch_order ? 1 : 0);
  
  if (this_pair->entOrSet[1-iface_no] == 0) { // iface2 is not sets
    iRel_processError(iBase_INVALID_ENTITY_HANDLE, "Expected Entity, got EntitySet");
    RETURN(iBase_INVALID_ENTITY_HANDLE);
  }
  
  int result = this_pair->get_eh_tags(iface_no, 
                           &ent1, 1, this_pair->assocTags[iface_no], 
                           reinterpret_cast<iBase_EntityHandle*>(ent2));
  
  RETURN(result);
}
void iRel_getSetEntAssociation (
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
    iRel_processError(iBase_FAILURE, "Didn't find association pair.");
    RETURN(iBase_FAILURE);
  }
  
  int iface_no = (switch_order ? 1 : 0);
  
  if (this_pair->entOrSet[1-iface_no] > 0) { // iface2 is sets
    iRel_processError(iBase_INVALID_ENTITY_HANDLE, "Expected EntitySet, got Entity");
    RETURN(iBase_INVALID_ENTITY_HANDLE);
  }
  
  int result = this_pair->get_eh_tags(iface_no, 
                           &ent1, 1, this_pair->assocTags[iface_no], ent2);
  
  RETURN(result);
}

void iRel_getSetSetAssociation (
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
    iRel_processError(iBase_FAILURE, "Didn't find association pair.");
    RETURN(iBase_FAILURE);
  }
  
  int iface_no = (switch_order ? 1 : 0);
  
  if (this_pair->entOrSet[1-iface_no] == 0) { // iface2 is not sets
    iRel_processError(iBase_INVALID_ENTITY_HANDLE, "Expected Entity, got EntitySet");
    RETURN(iBase_INVALID_ENTITY_HANDLE);
  }
  
  int result = this_pair->get_eh_tags(iface_no, 
                           &ent1, 1, this_pair->assocTags[iface_no], 
                           reinterpret_cast<iBase_EntityHandle*>(ent2));
  
  RETURN(result);
}

void iRel_getEntEntArrAssociation (
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
    iRel_processError(iBase_FAILURE, "Didn't find association pair.");
    RETURN(iBase_FAILURE);
  }
  int iface_no = (switch_order ? 1 : 0);
  
  if (this_pair->entOrSet[1-iface_no] > 0) { // iface2 is sets
    iBase_EntitySetHandle ent_set2 = 0;

    iRel_getEntSetAssociation(instance, rel,
                              ent1, switch_order,
                              &ent_set2, &result);
    if (iBase_SUCCESS != result) RETURN(result);

    result = this_pair->get_entities(1-iface_no, -1, 
                                     ent_set2,
                                     ent_array_2,
                                     ent_array_2_allocated,
                                     ent_array_2_size);
  }
  else {
      // otherwise put the set into the ent list
    CHECK_SIZE(iBase_EntityHandle, ent_array_2,
               ent_array_2_allocated, 1);
    *ent_array_2_size = 1;
    iRel_getEntEntAssociation(instance, rel, ent1, switch_order, 
                              *ent_array_2, &result );
  }
  
  RETURN(result);
}
void iRel_getSetEntArrAssociation (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle ent1,
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
    iRel_processError(iBase_FAILURE, "Didn't find association pair.");
    RETURN(iBase_FAILURE);
  }
  int iface_no = (switch_order ? 1 : 0);
  
  if (this_pair->entOrSet[1-iface_no] > 0) { // iface2 is sets
    iBase_EntitySetHandle ent_set2 = 0;

    iRel_getSetSetAssociation(instance, rel,
                              ent1, switch_order,
                              &ent_set2, &result);
    if (iBase_SUCCESS != result) RETURN(result);

    result = this_pair->get_entities(1-iface_no, -1, 
                                     ent_set2,
                                     ent_array_2,
                                     ent_array_2_allocated,
                                     ent_array_2_size);
  }
  else {
      // otherwise put the set into the ent list
    CHECK_SIZE(iBase_EntityHandle, ent_array_2,
               ent_array_2_allocated, 1);
    *ent_array_2_size = 1;
    iRel_getSetEntAssociation(instance, rel, ent1, switch_order, 
                              *ent_array_2, &result );
  }
  
  RETURN(result);
}

void iRel_getEntArrEntArrAssociation(
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

    // get assocpair so we can check whether ents are associated to
    // entities or sets or both
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  bool switched = false;
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find association pair.");
    RETURN(iBase_FAILURE);
  }
  for (int i = 0; i < ent_array_1_size; i++) {
    (*offset)[i] = tmp_ents.size();
    int tmp_result;
    iBase_EntityHandle *tmp_array = NULL;
    int tmp_array_size, tmp_array_allocated;
    iRel_getEntEntArrAssociation(instance, rel,
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
void iRel_getEntArrSetArrAssociation(
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle *ent_array_1,
  int ent_array_1_size,
  int switch_order,
  iBase_EntitySetHandle **ent_array_2,
  int *ent_array_2_allocated,
  int *ent_array_2_size,
  int *ierr)
{
  std::vector<iBase_EntityHandle> tmp_ents;
  int tmp_result, result = iBase_SUCCESS;
  CHECK_SIZE(iBase_EntitySetHandle, ent_array_2, ent_array_2_allocated, ent_array_1_size );
  *ent_array_2_size = ent_array_1_size;

    // get assocpair so we can check whether ents are associated to
    // entities or sets or both
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  bool switched = false;
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find association pair.");
    RETURN(iBase_FAILURE);
  }
  for (int i = 0; i < ent_array_1_size; i++) {
    iRel_getEntSetAssociation( instance, rel, ent_array_1[i],
                               switch_order, (*ent_array_2)+i,
                               &tmp_result );
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}
void iRel_getSetArrEntArrAssociation(
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle *ent_array_1,
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

    // get assocpair so we can check whether ents are associated to
    // entities or sets or both
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  bool switched = false;
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find association pair.");
    RETURN(iBase_FAILURE);
  }
  for (int i = 0; i < ent_array_1_size; i++) {
    (*offset)[i] = tmp_ents.size();
    int tmp_result;
    iBase_EntityHandle *tmp_array = NULL;
    int tmp_array_size, tmp_array_allocated;
    iRel_getSetEntArrAssociation(instance, rel,
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
void iRel_getSetArrSetArrAssociation(
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle *ent_array_1,
  int ent_array_1_size,
  int switch_order,
  iBase_EntitySetHandle **ent_array_2,
  int *ent_array_2_allocated,
  int *ent_array_2_size,
  int *ierr)
{
  std::vector<iBase_EntityHandle> tmp_ents;
  int tmp_result, result = iBase_SUCCESS;
  CHECK_SIZE(iBase_EntitySetHandle, ent_array_2, ent_array_2_allocated, ent_array_1_size );
  *ent_array_2_size = ent_array_1_size;

    // get assocpair so we can check whether ents are associated to
    // entities or sets or both
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  bool switched = false;
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find association pair.");
    RETURN(iBase_FAILURE);
  }
  for (int i = 0; i < ent_array_1_size; i++) {
    iRel_getSetSetAssociation( instance, rel, ent_array_1[i],
                               switch_order, (*ent_array_2)+i,
                               &tmp_result );
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}

void iRel_createVtxAndAssociate (
  iRel_Instance instance,
  double x,
  double y,
  double z,
  iBase_EntityHandle associatedGeomEnt,
  iBase_EntityHandle *new_entity_handle,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = lasso->find_pair(iRel_IGEOM_IFACE, iRel_IMESH_IFACE);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE,
                      "Didn't find a geometry-mesh association.");
    RETURN(iBase_FAILURE);
  }
  
    // first create the vtx
  int result = 
    this_pair->create_mesh_vertex(x, y, z, *new_entity_handle);
  if (iBase_SUCCESS != result) {
    iRel_processError(iBase_FAILURE, "Mesh vertex creation failed.");
    RETURN(iBase_FAILURE);
  }
  
    // now associate
  result = this_pair->set_assoc_tags(associatedGeomEnt, *new_entity_handle);
  RETURN(result);
}

void iRel_createEntAndAssociate (
  iRel_Instance instance,
  iMesh_EntityTopology new_entity_topology,
  iBase_EntityHandle *lower_order_entity_handles,
  int lower_order_entity_handles_size,
  iBase_EntityHandle associatedGeomEnt,
  iBase_EntityHandle *new_entity_handle,
  iBase_CreationStatus *status,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = lasso->find_pair(iRel_IGEOM_IFACE, iRel_IMESH_IFACE);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE,
                      "Didn't find a geometry-mesh association.");
    RETURN(iBase_FAILURE);
  }
  
    // first create the vtx
  int result = 
    this_pair->create_mesh_entity(new_entity_topology,
                                  lower_order_entity_handles,
                                  lower_order_entity_handles_size,
                                  *new_entity_handle,
                                  (int&)*status);
  if (iBase_SUCCESS != result) {
    iRel_processError(iBase_FAILURE, "Mesh entity creation failed.");
    RETURN(iBase_FAILURE);
  }
  
    // now associate
  result = this_pair->set_assoc_tags(associatedGeomEnt, *new_entity_handle);
  RETURN(result);
//
}

void iRel_createVtxArrAndAssociate (
  iRel_Instance instance,
  int num_verts,
  iBase_StorageOrder storage_order,
  double *new_coords,
  int new_coords_size,
  iBase_EntityHandle *associatedGeomEnts,
  iBase_EntityHandle **new_vertex_handles,
  int *new_vertex_handles_allocated,
  int *new_vertex_handles_size,
  int *ierr)
{
  if (num_verts != new_coords_size/3) {
    iRel_processError(iBase_INVALID_ARGUMENT, 
                      "Number of vertices should be 1/3 number of coords.");
    RETURN(iBase_INVALID_ARGUMENT);
  }
  
  CHECK_SIZE(iBase_EntityHandle, new_vertex_handles,
             new_vertex_handles_allocated,
             num_verts);

  int result = iBase_SUCCESS;
  for (int i = 0; i < num_verts; i++) {
    int tmp_result;
    iRel_createVtxAndAssociate(instance, 
                               new_coords[3*i],
                               new_coords[3*i+1],
                               new_coords[3*i+2],
                               associatedGeomEnts[i],
                               *new_vertex_handles+i, &tmp_result);
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}

void iRel_createEntArrAndAssociate (
  iRel_Instance instance,
  iMesh_EntityTopology new_entity_topology,
  iBase_EntityHandle *lower_order_entity_handles,
  int lower_order_entity_handles_size,
  int *offsets,
  int offsets_size,
  iBase_EntityHandle *associatedGeomEnts,
  iBase_EntityHandle **new_entity_handles,
  int *new_entity_handles_allocated,
  int *new_entity_handles_size,
  int **status,
  int *status_allocated,
  int *status_size,
  int *ierr)
{
  if (offsets[offsets_size-1] != lower_order_entity_handles_size-1) {
    iRel_processError(iBase_INVALID_ARGUMENT, 
                      "Last offset should be index of last element in lower_order_entity_handles.");
    RETURN(iBase_INVALID_ARGUMENT);
  }
  
  CHECK_SIZE(iBase_EntityHandle, new_entity_handles,
             new_entity_handles_allocated,
             offsets_size-1);
  *new_entity_handles_size = offsets_size-1;

  CHECK_SIZE(int, status,
             status_allocated,
             offsets_size-1);
  *status_size = offsets_size-1;

  int result = iBase_SUCCESS;
  for (int i = 0; i < offsets_size-1; i++) {
    int tmp_result;
    iRel_createEntAndAssociate(instance, 
                               new_entity_topology,
                               lower_order_entity_handles+offsets[i],
                               offsets[i+1]-offsets[i],
                               associatedGeomEnts[i],
                               *new_entity_handles+i,
                               (iBase_CreationStatus*)*status+i,
                               &tmp_result);
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  
  RETURN(result);
}

static void
iRel_inferArrArrAssociations(Lasso *lasso,
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
  int ents_size[2] = {ents1_size, ents2_size};
  int is_set[2] = {is_set1, is_set2};
  
  int result;


  std::vector<int> ents_gids, ents_dims;
  std::map<const int, iBase_EntityHandle> ents_gid_map[4];

  for (int i = 0; i < 2; i++) {
    ents_gids.resize(ents_size[i]);
    if (is_set[i]) 
      result = assoc_pair->get_int_tags(ents_index[i], 
                                        reinterpret_cast<iBase_EntitySetHandle*>(ents[i]), 
                                        ents_size[i], 
                                        assoc_pair->gidTags[ents_index[i]],
                                        &ents_gids[0]);
    else
      result = assoc_pair->get_int_tags(ents_index[i], ents[i], ents_size[i], 
                                        assoc_pair->gidTags[ents_index[i]],
                                        &ents_gids[0]);
    if (iBase_SUCCESS != result && iBase_TAG_NOT_FOUND != result) {
      RETURN(result);
    }

    ents_dims.resize(ents_size[i], -1);
    if (is_set[i]) {
      result = assoc_pair->get_int_tags(ents_index[i], 
                                        reinterpret_cast<iBase_EntitySetHandle*>(ents[i]), 
                                        ents_size[i], 
                                        assoc_pair->dimTags[ents_index[i]],
                                        &ents_dims[0]);
      if (iBase_SUCCESS != result && iBase_TAG_NOT_FOUND != result) {
        RETURN(result);
      }
    }

    else {
      int* ents_dims_ptr = &ents_dims[0], ents_dims_size = 0, ents_dims_alloc = ents_dims.size();
      result = 
        assoc_pair->get_ents_dims(ents_index[i], ents[i], ents_size[i], 
                                  &ents_dims_ptr, &ents_dims_alloc, &ents_dims_size);
      if (iBase_SUCCESS != result && iBase_TAG_NOT_FOUND != result) {
        RETURN(result);
      }
    }

    if (0 == i) {
      for (int j = 0; j < ents_size[i]; j++) {
        int dim = ents_dims[j];
        if (0 <= dim && 3 >= dim)
          ents_gid_map[dim][ents_gids[j]] = ents[i][j];
      }
    }
  }

  int dim;
  int c = !!is_set[0] + 2*!!is_set[1];
  for (int j = 0; j < ents_size[1]; j++) {
    dim = ents_dims[j];

      // only check entities for which the dimension entry is in a reasonable
      // range
    if (0 > dim || 3 < dim) continue;

      // there's a match if there's an entity with that dimension with matching id
    std::map<const int, iBase_EntityHandle>::iterator iter =
      ents_gid_map[dim].find(ents_gids[j]);

      // if it matches, set the association tags for those entities
    if (iter != ents_gid_map[dim].end()) {
      
      switch (c) {
        case 0: // !is_set[0] && !is_set[1]
          if (ents_index[0] == 0)
            result = assoc_pair->set_assoc_tags((*iter).second, ents[1][j]);
          else
            result = assoc_pair->set_assoc_tags(ents[1][j], (*iter).second);
          break;
        case 1: // is_set[0] && !is_set[1]
          if (ents_index[0] == 0)
            result = assoc_pair->set_assoc_tags((iBase_EntitySetHandle)(*iter).second, 
                                                 ents[1][j]);
          else
            result = assoc_pair->set_assoc_tags(ents[1][j],
                                                (iBase_EntitySetHandle)(*iter).second);
          break;
        case 2: // !is_set[0] && is_set[1]
          if (ents_index[0] == 0)
            result = assoc_pair->set_assoc_tags((*iter).second,
                                                (iBase_EntitySetHandle)ents[1][j]);
          else
            result = assoc_pair->set_assoc_tags((iBase_EntitySetHandle)ents[1][j],
                                                (*iter).second);
          break;
        case 3: // is_set[0] && is_set[1]
          if (ents_index[0] == 0)
            result = assoc_pair->set_assoc_tags((iBase_EntitySetHandle)(*iter).second,
                                                (iBase_EntitySetHandle)ents[1][j]);
          else
            result = assoc_pair->set_assoc_tags((iBase_EntitySetHandle)ents[1][j],
                                                (iBase_EntitySetHandle)(*iter).second);
          break;
      }
      if (iBase_SUCCESS != result) {
        RETURN(result);
      }
    }
  }

  RETURN(iBase_SUCCESS);
}
void
iRel_inferEntArrEntArrAssociations(Lasso *lasso,
                             AssocPair *assoc_pair,
                             iBase_EntityHandle *ents1,
                             const int ents1_size,
                             iBase_EntityHandle *ents2,
                             const int ents2_size,
                             int *ierr)
{
  iRel_inferArrArrAssociations( lasso, assoc_pair, ents1, ents1_size, 0,
                                ents2, ents2_size, 0, ierr );
}
void
iRel_inferEntArrSetArrAssociations(Lasso *lasso,
                             AssocPair *assoc_pair,
                             iBase_EntityHandle *ents1,
                             const int ents1_size,
                             iBase_EntitySetHandle *ents2,
                             const int ents2_size,
                             int *ierr)
{
  iRel_inferArrArrAssociations( lasso, assoc_pair, ents1, ents1_size, 0,
                                (iBase_EntityHandle*)ents2, ents2_size, 
                                1, ierr );
}
void
iRel_inferSetArrEntArrAssociations(Lasso *lasso,
                             AssocPair *assoc_pair,
                             iBase_EntitySetHandle *ents1,
                             const int ents1_size,
                             iBase_EntityHandle *ents2,
                             const int ents2_size,
                             int *ierr)
{
  iRel_inferArrArrAssociations( lasso, assoc_pair, 
                                (iBase_EntityHandle*)ents1, ents1_size, 
                                1, ents2, ents2_size, 0, ierr );
}
void
iRel_inferSetArrSetArrAssociations(Lasso *lasso,
                             AssocPair *assoc_pair,
                             iBase_EntitySetHandle *ents1,
                             const int ents1_size,
                             int is_set1,
                             iBase_EntitySetHandle *ents2,
                             const int ents2_size,
                             int is_set2,
                             int *ierr)

{
  iRel_inferArrArrAssociations( lasso, assoc_pair, 
                                (iBase_EntityHandle*)ents1, ents1_size, 
                                1, (iBase_EntityHandle*)ents2, 
                                ents2_size, 1, ierr );
}

void iRel_inferAllAssociations (
  iRel_Instance instance,
  iRel_RelationHandle rel,
  int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(instance);
  AssocPair *this_pair = reinterpret_cast<AssocPair*>(rel);
  if (NULL == this_pair) {
    iRel_processError(iBase_FAILURE, "Didn't find association pair.");
    RETURN(iBase_FAILURE);
  }

    // get all entities in those interfaces
  iBase_EntityHandle *ents1 = NULL, *ents2 = NULL;
  int ents1_size, ents2_size, ents1_alloc = 0, ents2_alloc = 0;
  int result;
  if (this_pair->entOrSet[0] > 0)
    result = this_pair->get_all_sets(0, (iBase_EntitySetHandle**)&ents1, &ents1_alloc, &ents1_size);
  else
    result = this_pair->get_all_entities(0, -1, &ents1, &ents1_alloc, &ents1_size);
  if (iBase_SUCCESS != result) {
    RETURN(result);
  }
  
  if (this_pair->entOrSet[1] > 0)
    result = this_pair->get_all_sets(1, (iBase_EntitySetHandle**)&ents2, &ents2_alloc, &ents2_size);
  else
  result = this_pair->get_all_entities(1, -1, &ents2, &ents2_alloc, &ents2_size);
  if (iBase_SUCCESS != result) {
    RETURN(result);
  }
  
  iRel_inferArrArrAssociations(lasso, this_pair, 
                               ents1, ents1_size, this_pair->entOrSet[0],
                               ents2, ents2_size, this_pair->entOrSet[1], &result);

  free(ents1); free(ents2);
  RETURN(result);
}

void iRel_inferEntAssociations (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle entity,
  int is_set,
  int iface_no,
  int *ierr)
{
  iRel_inferEntArrAssociations(instance, rel,
                            &entity, 1, iface_no, 
                            ierr);
}
void iRel_inferSetAssociations (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle entity,
  int iface_no,
  int *ierr)
{
  iRel_inferSetArrAssociations(instance, rel,
                            &entity, 1, iface_no, 
                            ierr);
}

static void iRel_inferArrAssociations (
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
    iRel_processError(iBase_FAILURE, "Didn't find association pair.");
    RETURN(iBase_FAILURE);
  }

  if (0 > iface_no || 1 < iface_no) {
    iRel_processError(iBase_INVALID_ARGUMENT, "Interface number must be 0 or 1");
    RETURN(iBase_INVALID_ARGUMENT);
  }
  else if ((is_set && this_pair->entOrSet[iface_no] == 0) ||
      (!is_set && this_pair->entOrSet[iface_no] > 0)) {
    iRel_processError(iBase_INVALID_ARGUMENT, "is_set must match entOrSet in call to inferArrAssociations");
    RETURN(iBase_INVALID_ARGUMENT);
  }
  
    // get all entities in iface2
  iBase_EntityHandle *ents2 = NULL, *ents1 = entities;
  int ents1_size = entities_size, ents2_size, ents2_alloc = 0;
  int result;
  if (this_pair->entOrSet[1-iface_no] > 0)
    result = this_pair->get_all_sets(1-iface_no, 
                                (iBase_EntitySetHandle**)&ents2, 
                                &ents2_alloc, &ents2_size);
  else
    result = this_pair->get_all_entities(1-iface_no, -1, 
                                &ents2, &ents2_alloc, &ents2_size);
  if (iBase_SUCCESS != result) {
    RETURN(result);
  }

    // switch so that entity lists always go into inferArrArrAssoc in
    // forward order wrt assoc_pair
  if (1 == iface_no) {
    ents1_size = ents2_size;
    ents2_size = entities_size;
    ents1 = ents2;
    ents2 = entities;
  }
  
  iRel_inferArrArrAssociations(lasso, this_pair,
                               ents1, ents1_size, 
                               this_pair->entOrSet[0],
                               ents2, ents2_size, 
                               this_pair->entOrSet[1], &result);

  if (1 == iface_no) free(ents1);
  else free(ents2);
  
  RETURN(result);
}

void iRel_inferEntArrAssociations (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntityHandle *entities,
  int entities_size,
  int iface_no,
  int *ierr)
{
  iRel_inferArrAssociations( instance, rel, entities, entities_size, 0, iface_no, ierr );
}
void iRel_inferSetArrAssociations (
  iRel_Instance instance,
  iRel_RelationHandle rel,    
  iBase_EntitySetHandle *entities,
  int entities_size,
  int iface_no,
  int *ierr)
{
  iRel_inferArrAssociations( instance, rel, (iBase_EntityHandle*)entities, entities_size, 1, iface_no, ierr );
}

void 
iRel_moveTo(iRel_Instance assoc,
            /* in */ iGeom_Instance geom, /* in */ iMesh_Instance mesh,
            /* in */ iBase_EntityHandle gent,
            int *ierr)
{
  Lasso *lasso = reinterpret_cast<Lasso*>(assoc);
  bool switched = false;
  AssocPair *this_pair = lasso->find_pair(geom, mesh, &switched);

    // get mesh entity coresponding to the given geom entity
  int result;

    // get the vertices in the set, inclusive
  iBase_EntityHandle *vertices = NULL;
  int vertices_size = 0, vertices_allocated = 0;
  if (this_pair->entOrSet[1-switched]) {
    iBase_EntitySetHandle ent_set;
    iRel_getEntSetAssociation(assoc, (iRel_RelationHandle)this_pair,
                              gent, switched, &ent_set, 
                              &result);
    if (iBase_SUCCESS != result) {
      iRel_processError(result, "Trouble getting associated set for a gentity.");
      RETURN(result);
    }
  
    result = this_pair->get_entities((switched ? 0 : 1), 0, 
                      reinterpret_cast<iBase_EntitySetHandle>(ent_set),
                                     &vertices, &vertices_allocated, 
                                     &vertices_size);
    if (iBase_SUCCESS != result) {
      iRel_processError(result, "Trouble getting vertices in set for a gentity.");
      RETURN(result);
    }
  }
  else {
      // just a vertex, we hope - put that in array
    CHECK_SIZE(iBase_EntityHandle, &vertices, &vertices_allocated, 1);
    vertices_size = 1;
    iRel_getEntEntAssociation(assoc, (iRel_RelationHandle)this_pair,
                              gent, switched, vertices, 
                              &result);
    if (iBase_SUCCESS != result) {
      iRel_processError(result, "Trouble getting associated set for a gentity.");
      RETURN(result);
    }
  }

    // check to see if there's nothing to do; if so, THIS IS NOT AN ERROR - e.g. curve with
    // 1 edge but no vertices
  if (vertices_size == 0) {
    RETURN(iBase_SUCCESS);
  }
  
    // get the coordinates for these vertices
  double *coords = NULL;
  int coords_size = 0, coords_allocated = 0;
  result = this_pair->get_mesh_coords(vertices, vertices_size,
                                      &coords, &coords_allocated, &coords_size, 
                                      iBase_INTERLEAVED);
  if (iBase_SUCCESS != result) {
    char this_descr[120];
    iMesh_getDescription(mesh, this_descr, ierr, 120);
    std::string this_str(this_descr);
    this_str += "Trouble getting vertex coordinates.";
    iRel_processError(result, this_str.c_str());
    RETURN(result);
  }
  
    // project these points to the gentity; reuse coords array
  std::vector<iBase_EntityHandle> gents(vertices_size);
  std::fill(gents.begin(), gents.end(), gent);

  result = this_pair->get_closest_pt(&gents[0], vertices_size, iBase_INTERLEAVED,
                                     coords, coords_size,
                                     &coords, &coords_allocated, &coords_size);
  if (iBase_SUCCESS != result) {
    char this_descr[120];
    iMesh_getDescription(mesh, this_descr, ierr, 120);
    std::string this_str(this_descr);
    this_str += "Trouble getting closest point for vertices.";
    iRel_processError(result, this_str.c_str());
    RETURN(result);
  }

    // now set the coordinates on the vertices
    //order = iBase::StorageOrder_INTERLEAVED;
  result = this_pair->set_mesh_coords(vertices, vertices_size, 
                                      coords, iBase_INTERLEAVED);
  if (iBase_SUCCESS != result) {
    char this_descr[120];
    iMesh_getDescription(mesh, this_descr, ierr, 120);
    std::string this_str(this_descr);
    this_str += "Trouble setting vertex coordinates.";
    iRel_processError(result, this_str.c_str());
    RETURN(result);
  }

  free(vertices);
  free(coords);
  
  RETURN(iBase_SUCCESS);
}

void iRel_newAssoc(/* in */ const char * /* options */,
                   iRel_Instance *instance,
                   int *ierr,
                   const int options_len) 
{
  if (0 != options_len) {
    iRel_processError(iBase_NOT_SUPPORTED, "No options for iRel factory have been implemented.");
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
AssocPair *Lasso::find_pair(IfaceType type1, IfaceType type2) 
{
  if (type1 > type2) {
    IfaceType tmp_type = type1;
    type1 = type2;
    type2 = tmp_type;
  }
  
  for (std::vector<AssocPair*>::iterator vit = assocPairs.begin();
       vit != assocPairs.end(); vit++) {
    if ((*vit)->ifaceTypes[0] == type1 && 
        (*vit)->ifaceTypes[0] == type2)
      return *vit;
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


  
