#include "AssocPair.hpp"
#include "Lasso.hpp"
#include <stdlib.h>

const char *AssocPair::GLOBAL_ID_TAG_NAME = "GLOBAL_ID";
const char *AssocPair::GEOM_DIMENSION_TAG_NAME = "GEOM_DIMENSION";
const char *AssocPair::ASSOCIATION_TAG_NAME = "ASSOCIATION";
int AssocPair::currId = 0;

AssocPair::AssocPair(iRel_Instance instance,
                     RelationType ent_or_set0, IfaceType type0,
                     RelationType ent_or_set1, IfaceType type1)
  : instance(instance)
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

  pairId = currId++;
}

AssocPair::~AssocPair()
{}

int AssocPair::create_tags()
{
  // first the association tag name
  std::stringstream tag_name_str;
  tag_name_str << ASSOCIATION_TAG_NAME << pairId;
  std::string tmp_name(tag_name_str.str());

  // create the tag in each interface
  iBase_EntityHandle def_handle_value = 0;
  int def_int_value = 0;

  for (int i = 0; i < 2; i++) {
    CHK_ERRORR( tag_get_handle(i, tmp_name.c_str(), 1, iBase_ENTITY_HANDLE,
                               true, &def_handle_value, assocTags+i) );
    CHK_ERRORR( tag_get_handle(i, GLOBAL_ID_TAG_NAME, 1, iBase_INTEGER,
                               true, &def_int_value, gidTags+i) );

    if (ifaceTypes[ i] == iRel_IMESH_IFACE &&
        ifaceTypes[!i] == iRel_IGEOM_IFACE) {
      CHK_ERRORR( tag_get_handle(i, GEOM_DIMENSION_TAG_NAME, 1, iBase_INTEGER,
                                 false, NULL, dimTags+i) );
    }
  }

  RETURNR(iBase_SUCCESS);
}

int AssocPair::destroy_tags()
{
  // We assume that subclasses have called create_tags
  for (int i = 0; i < 2; i++) {
    tag_destroy(i, assocTags[i]);
  }

  RETURNR(iBase_SUCCESS);
}

int AssocPair::set_assoc_tags(iBase_EntityHandle ent1, iBase_EntityHandle ent2)
{
  // check that is_setx is consistent with entOrSet
  if (entOrSet[0] != iRel_ENTITY || entOrSet[1] != iRel_ENTITY)
    ERRORR(iBase_FAILURE, "Invalid relation type");

  // check that if we're passing in an ent for a 'both'-type
  // assoc, there's already a set associated to the other ent
  iBase_EntityHandle tmp_ent;
  if (entOrSet[0] == iRel_BOTH && get_tags(1, &ent2, 1, assocTags[1], &tmp_ent,
                                           sizeof(tmp_ent)) != iBase_SUCCESS)
    ERRORR(iBase_FAILURE, "Couldn't find associated set on left side");
  if (entOrSet[1] == iRel_BOTH && get_tags(0, &ent1, 1, assocTags[0], &tmp_ent,
                                           sizeof(tmp_ent)) != iBase_SUCCESS)
    ERRORR(iBase_FAILURE, "Couldn't find associated set on right side");

  // set ent1 assoc tag to point to ent2
  if (entOrSet[1] == iRel_ENTITY) {
    CHK_ERRORR( set_tags(0, &ent1, 1, assocTags[0], &ent2, sizeof(ent2)) );
  }

  // set ent2 assoc tag to point to ent1
  if (entOrSet[0] == iRel_ENTITY) {
    CHK_ERRORR( set_tags(1, &ent2, 1, assocTags[1], &ent1, sizeof(ent1)) );
  }

  RETURNR(iBase_SUCCESS);
}

int AssocPair::set_assoc_tags(iBase_EntityHandle ent1,
                              iBase_EntitySetHandle set2)
{
  // check that is_setx is consistent with entOrSet
  if (entOrSet[0] != iRel_ENTITY || entOrSet[1] == iRel_ENTITY)
    ERRORR(iBase_FAILURE, "Invalid relation type");

  // check that if we're passing in an ent for a 'both'-type
  // assoc, there's already a set associated to the other ent
  iBase_EntityHandle tmp_ent;
  if (entOrSet[0] == iRel_BOTH && get_tags(1, &set2, 1, assocTags[1], &tmp_ent,
                                           sizeof(tmp_ent)) != iBase_SUCCESS)
    ERRORR(iBase_FAILURE, "Couldn't find associated set on left side");

  // set ent1 assoc tag to point to ent2
  CHK_ERRORR( set_tags(0, &ent1, 1, assocTags[0], &set2, sizeof(set2)) );

  // set ent2 assoc tag to point to ent1
  if (entOrSet[0] == iRel_ENTITY) {
    CHK_ERRORR( set_tags(1, &set2, 1, assocTags[1], &ent1, sizeof(ent1)) );
  }

  // if either are sets and 'both' type association, get the indiv
  // entities & set associations for them too
  if (entOrSet[1] == iRel_BOTH) {
    iBase_EntityHandle *entities = NULL, to_ent;
    int entities_alloc = 0, entities_size, iface_no;

    // get ents from set2 & associate to ent1
    CHK_ERRORR( get_entities(1, -1, set2, &entities, &entities_alloc,
                             &entities_size) );
    to_ent = ent1;
    iface_no = 1;

    for (int i = 0; i < entities_size; i++) {
      CHK_ERRORR( set_tags(iface_no, entities+i, 1, assocTags[iface_no],
                           &to_ent, sizeof(to_ent)) );
    }

    free(entities);
  }

  RETURNR(iBase_SUCCESS);
}

int AssocPair::set_assoc_tags(iBase_EntitySetHandle set1,
                              iBase_EntityHandle ent2)
{
  // check that is_setx is consistent with entOrSet
  if (entOrSet[0] == iRel_ENTITY || entOrSet[1] != iRel_ENTITY)
    ERRORR(iBase_FAILURE, "Invalid relation type");

  // check that if we're passing in an ent for a 'both'-type
  // assoc, there's already a set associated to the other ent
  iBase_EntityHandle tmp_ent;
  if (entOrSet[1] == iRel_BOTH && get_tags(0, &set1, 1, assocTags[0], &tmp_ent,
                                           sizeof(tmp_ent)) != iBase_SUCCESS)
    ERRORR(iBase_FAILURE, "Couldn't find associated set on right side");

  // set ent1 assoc tag to point to ent2
  if (entOrSet[1] == iRel_ENTITY) {
    CHK_ERRORR( set_tags(0, &set1, 1, assocTags[0], &ent2, sizeof(ent2)) );
  }

  // set ent2 assoc tag to point to ent1
  CHK_ERRORR( set_tags(1, &ent2, 1, assocTags[1], &set1, sizeof(set1)) );

  // if either are sets and 'both' type association, get the indiv
  // entities & set associations for them too
  if (entOrSet[0] == iRel_BOTH) {
    iBase_EntityHandle *entities = NULL, to_ent;
    int entities_alloc, entities_size, iface_no;
      // get ents from set1 & associate to ent2
    CHK_ERRORR( get_entities(0, -1, set1, &entities, &entities_alloc,
                             &entities_size) );
    to_ent = ent2;
    iface_no = 0;

    for (int i = 0; i < entities_size; i++) {
      CHK_ERRORR( set_tags(iface_no, entities+i, 1, assocTags[iface_no],
                           &to_ent, sizeof(to_ent)) );
    }

    free(entities);
  }

  RETURNR(iBase_SUCCESS);
}

int AssocPair::set_assoc_tags(iBase_EntitySetHandle set1,
                              iBase_EntitySetHandle set2)
{
  // check that is_setx is consistent with entOrSet
  if (entOrSet[0] == iRel_ENTITY || entOrSet[1] == iRel_ENTITY)
    ERRORR(iBase_FAILURE, "Invalid relation type");

  // set ent1 assoc tag to point to ent2
  CHK_ERRORR( set_tags(0, &set1, 1, assocTags[0], &set2, sizeof(set2)) );

  // set ent2 assoc tag to point to ent1
  CHK_ERRORR( set_tags(1, &set2, 1, assocTags[1], &set1, sizeof(set1)) );

  RETURNR(iBase_SUCCESS);
}

int AssocPair::get_assoc_tags(const int iface_no, iBase_EntityHandle *entities,
                              int num_entities, iBase_EntityHandle *tag_values)
{
  if (entOrSet[!iface_no] != iRel_ENTITY) { // other iface is sets
    ERRORR(iBase_INVALID_ENTITY_HANDLE, "Expected EntitySet, got Entity");
  }

  return get_tags(iface_no, entities, num_entities, assocTags[iface_no],
                  tag_values, sizeof(iBase_EntityHandle));
}

int AssocPair::get_assoc_tags(const int iface_no, iBase_EntitySetHandle *sets,
                              int num_sets, iBase_EntityHandle *tag_values)
{
  if (entOrSet[!iface_no] != iRel_ENTITY) { // other iface is sets
    ERRORR(iBase_INVALID_ENTITY_HANDLE, "Expected EntitySet, got Entity");
  }

  return get_tags(iface_no, sets, num_sets, assocTags[iface_no], tag_values,
                  sizeof(iBase_EntityHandle));
}

int AssocPair::get_assoc_tags(const int iface_no, iBase_EntityHandle *entities,
                              int num_entities,
                              iBase_EntitySetHandle *tag_values)
{
  if (entOrSet[!iface_no] == iRel_ENTITY) { // other iface is not sets
    ERRORR(iBase_INVALID_ENTITY_HANDLE, "Expected Entity, got EntitySet");
  }

  return get_tags(iface_no, entities, num_entities, assocTags[iface_no],
                  reinterpret_cast<iBase_EntityHandle*>(tag_values),
                  sizeof(iBase_EntityHandle));
}

int AssocPair::get_assoc_tags(const int iface_no, iBase_EntitySetHandle *sets,
                              int num_sets, iBase_EntitySetHandle *tag_values)
{
  if (entOrSet[!iface_no] == iRel_ENTITY) { // other iface is not sets
    ERRORR(iBase_INVALID_ENTITY_HANDLE, "Expected Entity, got EntitySet");
  }

  return get_tags(iface_no, sets, num_sets, assocTags[iface_no],
                  reinterpret_cast<iBase_EntityHandle*>(tag_values),
                  sizeof(iBase_EntityHandle));
}

int AssocPair::get_assoc_tags(const int iface_no, iBase_EntityHandle *entities,
                              int num_entities,
                              iBase_EntityIterator *tag_values)
{
  std::vector<iBase_EntitySetHandle> sets(num_entities);
  CHK_ERRORR( get_assoc_tags(iface_no, entities, num_entities, &sets[0]) );

  for(int i=0; i<num_entities; i++) {
    CHK_ERRORR( get_iterator(iface_no, sets[i], &tag_values[i]) );
  }

  RETURNR(iBase_SUCCESS);
}

int AssocPair::get_assoc_tags(const int iface_no, iBase_EntitySetHandle *sets,
                              int num_sets, iBase_EntityIterator *tag_values)
{
  std::vector<iBase_EntitySetHandle> sets2(num_sets);
  CHK_ERRORR( get_assoc_tags(iface_no, sets, num_sets, &sets2[0]) );

  for(int i=0; i<num_sets; i++) {
    CHK_ERRORR( get_iterator(iface_no, sets2[i], &tag_values[i]) );
  }

  RETURNR(iBase_SUCCESS);
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
