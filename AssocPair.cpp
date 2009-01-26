#include "AssocPair.hpp"
#include "Lasso.hpp"
#include <assert.h>
#include <stdlib.h>

const char *AssocPair::GLOBAL_ID_TAG_NAME = "GLOBAL_ID";
const char *AssocPair::GEOM_DIMENSION_TAG_NAME = "GEOM_DIMENSION";
const char *AssocPair::ASSOCIATION_TAG_NAME = "ASSOCIATION";

AssocPair::AssocPair(const int ent_or_set0, const int ent_or_set1,
                            Lasso *lasso)
    : 
    myLasso(lasso)
{
  ifaceTypes[0] = iRel_IBASE_IFACE;
  ifaceTypes[1] = iRel_IBASE_IFACE;
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

int AssocPair::set_assoc_tags(iBase_EntityHandle ent1,
                                          int is_set1,
                                          iBase_EntityHandle ent2,
                                          int is_set2) 
{
  iBase_EntityHandle tmp_ent;
    // check that is_setx is consistent with entOrSet
  assert(((is_set1 && entOrSet[0] > 0) ||
          (!is_set1 && entOrSet[0] == 0)) &&
         ((is_set2 && entOrSet[1] > 0) ||
          (!is_set2 && entOrSet[1] == 0)));

    // check that if we're passing in an ent for a 'both'-type
    // assoc, there's already a set associated to the other ent
  assert((entOrSet[0] < 2 || is_set1 ||
          get_eh_tags(1, &ent2, 1, assocTags[1],
                      is_set2, &tmp_ent) == iBase_SUCCESS) &&
         (entOrSet[1] < 2 || is_set2 ||
          get_eh_tags(0, &ent1, 1, assocTags[0],
                      is_set1, &tmp_ent) == iBase_SUCCESS));
  
  int result = iBase_SUCCESS;

    // set ent1 assoc tag to point to ent2
  if (is_set2 || entOrSet[1] == 0) {
    result = set_eh_tags(0, &ent1, 1, assocTags[0],
                         (is_set1 > 0), &ent2);
    if (iBase_SUCCESS != result) return result;
  }
  
    // set ent2 assoc tag to point to ent1
  if (is_set1 || entOrSet[0] == 0) {
    result = set_eh_tags(1, &ent2, 1, assocTags[1],
                         (is_set2 > 0), &ent1);
    if (iBase_SUCCESS != result) return result;
  }

    // if either are sets and 'both' type association, get the indiv
    // entities & set associations for them too
  iBase_EntityHandle *entities = NULL, to_ent;
  int entities_alloc, entities_size, iface_no;
  if (!is_set1 && is_set2 && entOrSet[1] == 2) {
      // get ents from set2 & associate to ent1
    result = get_entities(1, -1, ent2, 
                          &entities, &entities_alloc, &entities_size);
    if (iBase_SUCCESS != result) return result;
    to_ent = ent1;
    iface_no = 1;
  }
  else if (!is_set2 && is_set1 && entOrSet[0] == 2) {
      // get ents from set1 & associate to ent2
    result = get_entities(0, -1, ent1, 
                          &entities, &entities_alloc, &entities_size);
    if (iBase_SUCCESS != result) return result;
    to_ent = ent2;
    iface_no = 0;
  }
  else return result;

  for (int i = 0; i < entities_size; i++) {
    result = set_eh_tags(iface_no, entities+i, 1,
                         assocTags[iface_no],
                         false, &to_ent);
    if (iBase_SUCCESS != result) return result;
  }

  free(entities);
  
  return result;
}

IfaceType AssocPair::iface_type(const int iface_no) 
{
  return ifaceTypes[iface_no];
}

  
