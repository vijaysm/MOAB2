#include "AssocPair.hpp"
#include "Lasso.hpp"
#include <assert.h>
#include <stdlib.h>

const char *AssocPair::GLOBAL_ID_TAG_NAME = "GLOBAL_ID";
const char *AssocPair::GEOM_DIMENSION_TAG_NAME = "GEOM_DIMENSION";
const char *AssocPair::ASSOCIATION_TAG_NAME = "ASSOCIATION";

AssocPair::AssocPair(RelationType ent_or_set0, IfaceType type0,
                     RelationType ent_or_set1, IfaceType type1,
                     Lasso *lasso) : myLasso(lasso)
{
  ifaceTypes[0] = type0;
  ifaceTypes[1] = type1;
  entOrSet[0] = ent_or_set0;
  entOrSet[1] = ent_or_set1;
  assocTags[0] = 0;
  assocTags[1] = 0;
  gidTags[0] = 0;
  gidTags[1] = 0;
  dimTags[0] = 0;
  dimTags[1] = 0;
  
  pairId = myLasso->assocPairs.size();
  myLasso->assocPairs.push_back(this);
}

AssocPair::~AssocPair() 
{
    // don't take out of Lasso's list here - that's done by
    // Lasso's destructor
}

int AssocPair::create_tags() 
{
    // first the association tag name
  std::stringstream tag_name_str;
  tag_name_str << ASSOCIATION_TAG_NAME << pairId;
  std::string tmp_name(tag_name_str.str());
  
    // create the tag in each interface
  iBase_EntityHandle def_handle_value = 0;
  int def_int_value = 0;

  for (int i = 0; i <= 1; i++) {
    assocTags[i] = tag_get_handle(i, tmp_name.c_str(), 1,
                                  iBase_ENTITY_HANDLE, true, &def_handle_value);

    gidTags[i] = tag_get_handle(i, GLOBAL_ID_TAG_NAME, 1, iBase_INTEGER, 
                                true, &def_int_value);

    if (1 == i && ifaceTypes[0] == iRel_IGEOM_IFACE && ifaceTypes[1] == iRel_IMESH_IFACE) {
      dimTags[i] = tag_get_handle(i, GEOM_DIMENSION_TAG_NAME, 1, iBase_INTEGER, 
                                  false);
    }
  }
  
  return iBase_SUCCESS;
}

int AssocPair::set_assoc_tags(iBase_EntityHandle ent1, iBase_EntityHandle ent2)
{
    // check that is_setx is consistent with entOrSet
  assert(entOrSet[0] == 0 && entOrSet[1] == 0);

    // check that if we're passing in an ent for a 'both'-type
    // assoc, there's already a set associated to the other ent
  assert((entOrSet[0] < 2 || 
          get_eh_tags(1, &ent2, 1, assocTags[1], &tmp_ent) == iBase_SUCCESS) &&
         (entOrSet[1] < 2 || 
          get_eh_tags(0, &ent1, 1, assocTags[0], &tmp_ent) == iBase_SUCCESS));
  
  int result = iBase_SUCCESS;

    // set ent1 assoc tag to point to ent2
  if (entOrSet[1] == 0) {
    result = set_tags(0, &ent1, 1, assocTags[0], &ent2,
                      sizeof(iBase_EntityHandle));
    if (iBase_SUCCESS != result) return result;
  }
  
    // set ent2 assoc tag to point to ent1
  if (entOrSet[0] == 0) {
    result = set_tags(1, &ent2, 1, assocTags[1], &ent1,
                       sizeof(iBase_EntityHandle));
    if (iBase_SUCCESS != result) return result;
  }

  return result;
}

int AssocPair::set_assoc_tags(iBase_EntityHandle ent1,
                              iBase_EntitySetHandle set2) 
{
    // check that is_setx is consistent with entOrSet
  assert(entOrSet[0] == 0 && entOrSet[1] > 0);

    // check that if we're passing in an ent for a 'both'-type
    // assoc, there's already a set associated to the other ent
  assert( entOrSet[0] < 2 ||
          get_eh_tags(1, &set2, 1, assocTags[1], &tmp_ent) == iBase_SUCCESS);
  
  int result = iBase_SUCCESS;

    // set ent1 assoc tag to point to ent2
  result = set_tags(0, &ent1, 1, assocTags[0], &set2,
                    sizeof(iBase_EntityHandle));
  if (iBase_SUCCESS != result) return result;
  
    // set ent2 assoc tag to point to ent1
  if (entOrSet[0] == 0) {
    result = set_tags(1, &set2, 1, assocTags[1], &ent1,
                      sizeof(iBase_EntityHandle));
    if (iBase_SUCCESS != result) return result;
  }

    // if either are sets and 'both' type association, get the indiv
    // entities & set associations for them too
  if (entOrSet[1] == 2) {
    iBase_EntityHandle *entities = NULL, to_ent;
    int entities_alloc = 0, entities_size, iface_no;

      // get ents from set2 & associate to ent1
    result = get_entities(1, -1, set2, &entities, &entities_alloc,
                          &entities_size);
    if (iBase_SUCCESS != result) return result;
    to_ent = ent1;
    iface_no = 1;

    for (int i = 0; i < entities_size; i++) {
      result = set_tags(iface_no, entities+i, 1, assocTags[iface_no], &to_ent,
                        sizeof(iBase_EntityHandle));
      if (iBase_SUCCESS != result) return result;
    }

    free(entities);
  }

  return result;
}

int AssocPair::set_assoc_tags(iBase_EntitySetHandle set1,
                              iBase_EntityHandle ent2) 
{
    // check that is_setx is consistent with entOrSet
  assert(entOrSet[0] > 0 && entOrSet[1] == 0);

    // check that if we're passing in an ent for a 'both'-type
    // assoc, there's already a set associated to the other ent
  assert( entOrSet[1] < 2 ||
          get_eh_tags(0, &ent1, 1, assocTags[0], &tmp_ent) == iBase_SUCCESS);
  
  int result = iBase_SUCCESS;

    // set ent1 assoc tag to point to ent2
  if (entOrSet[1] == 0) {
    result = set_tags(0, &set1, 1, assocTags[0], &ent2,
                      sizeof(iBase_EntityHandle));
    if (iBase_SUCCESS != result) return result;
  }
  
    // set ent2 assoc tag to point to ent1
  result = set_tags(1, &ent2, 1, assocTags[1], &set1,
                    sizeof(iBase_EntityHandle));
  if (iBase_SUCCESS != result) return result;

    // if either are sets and 'both' type association, get the indiv
    // entities & set associations for them too
  if (entOrSet[0] == 2) {
    iBase_EntityHandle *entities = NULL, to_ent;
    int entities_alloc, entities_size, iface_no;
      // get ents from set1 & associate to ent2
    result = get_entities(0, -1, set1, &entities, &entities_alloc,
                          &entities_size);
    if (iBase_SUCCESS != result) return result;
    to_ent = ent2;
    iface_no = 0;

    for (int i = 0; i < entities_size; i++) {
      result = set_tags(iface_no, entities+i, 1, assocTags[iface_no], &to_ent,
                        sizeof(iBase_EntityHandle));
      if (iBase_SUCCESS != result) return result;
    }

    free(entities);
  }
  
  return result;
}

int AssocPair::set_assoc_tags(iBase_EntitySetHandle set1,
                              iBase_EntitySetHandle set2) 
{
    // check that is_setx is consistent with entOrSet
  assert(entOrSet[0] > 0 && entOrSet[1] > 0);

  int result = iBase_SUCCESS;

    // set ent1 assoc tag to point to ent2
  result = set_tags(0, &set1, 1, assocTags[0], &set2,
                    sizeof(iBase_EntityHandle));
  if (iBase_SUCCESS != result) return result;
  
    // set ent2 assoc tag to point to ent1
  result = set_tags(1, &set2, 1, assocTags[1], &set1,
                    sizeof(iBase_EntityHandle));
  if (iBase_SUCCESS != result) return result;
  
  return result;
}

int AssocPair::get_assoc_tags(const int iface_no, iBase_EntityHandle *entities,
                              int num_entities, iBase_EntityHandle *tag_values)
{
  return get_tags(iface_no, entities, num_entities, assocTags[iface_no],
                  tag_values, sizeof(iBase_EntityHandle));
}

int AssocPair::get_assoc_tags(const int iface_no, iBase_EntitySetHandle *sets,
                              int num_sets, iBase_EntityHandle *tag_values)
{
  return get_tags(iface_no, sets, num_sets, assocTags[iface_no], tag_values,
                  sizeof(iBase_EntityHandle));
}

int AssocPair::get_assoc_tags(const int iface_no, iBase_EntityHandle *entities,
                              int num_entities,
                              iBase_EntitySetHandle *tag_values)
{
  return get_assoc_tags(iface_no, entities, num_entities,
                        reinterpret_cast<iBase_EntityHandle*>(tag_values));
}

int AssocPair::get_assoc_tags(const int iface_no, iBase_EntitySetHandle *sets,
                              int num_sets, iBase_EntitySetHandle *tag_values)
{
  return get_assoc_tags(iface_no, sets, num_sets,
                        reinterpret_cast<iBase_EntityHandle*>(tag_values));
}

int AssocPair::get_gid_tags(const int iface_no, iBase_EntityHandle *entities,
                            int num_entities, int *tag_values)
{
  return get_tags(iface_no, entities, num_entities, gidTags[iface_no],
                  tag_values, sizeof(int));
}

int AssocPair::get_gid_tags(const int iface_no, iBase_EntitySetHandle *sets,
                            int num_sets, int *tag_values)
{
  return get_tags(iface_no, sets, num_sets, gidTags[iface_no], tag_values,
                  sizeof(int));
}

int AssocPair::get_dim_tags(const int iface_no, iBase_EntityHandle *entities,
                            int num_entities, int *tag_values)
{
  return get_tags(iface_no, entities, num_entities, dimTags[iface_no],
                  tag_values, sizeof(int));
}

int AssocPair::get_dim_tags(const int iface_no, iBase_EntitySetHandle *sets,
                            int num_sets, int *tag_values)
{
  return get_tags(iface_no, sets, num_sets, dimTags[iface_no], tag_values,
                  sizeof(int));
}

IfaceType AssocPair::iface_type(const int iface_no) 
{
  return ifaceTypes[iface_no];
}

RelationType AssocPair::ent_or_set(const int iface_no) 
{
  return entOrSet[iface_no];
}
  
bool AssocPair::equivalent(IfaceType type1, 
                           IfaceType type2,
                           bool *order_switched) 
{
  bool equiv = false;
  if (type1 == ifaceTypes[0] && 
      type2 == ifaceTypes[1]) {
    if (NULL != order_switched) *order_switched = false;
    equiv = true;
  }
  else if (type1 == ifaceTypes[1] && 
           type2 == ifaceTypes[0]) {
    equiv = true;
    if (NULL != order_switched) *order_switched = true;
  }

  return equiv;
}
