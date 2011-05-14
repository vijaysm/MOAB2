#ifndef ASSOCPAIR_HPP
#define ASSOCPAIR_HPP

#include "iRel_Lasso.hpp"
#include <sstream>

class Lasso;

class AssocPair
{
public:
  friend class Lasso;

  iRel_IfaceType iface_type(const int iface_no);
  iRel_RelationType ent_or_set(const int iface_no);

  virtual iBase_Instance iface_instance(const int iface_no) = 0;

  bool equivalent(iRel_IfaceType type1, iRel_IfaceType type2,
                  bool *order_switched = NULL);

  virtual bool equivalent(iBase_Instance iface1, iBase_Instance iface2,
                          bool *order_switched = NULL) = 0;

  virtual bool contains(iBase_Instance iface) = 0;

  virtual bool same_interface(iBase_Instance iface1, iBase_Instance iface2) = 0;

  virtual int get_all_entities(const int iface_no,
                               const int dimension,
                               iBase_EntityHandle **entities,
                               int *entities_alloc,
                               int *entities_size) = 0;

  virtual int get_all_sets(const int iface_no,
                           iBase_EntitySetHandle **entities,
                           int *entities_alloc,
                           int *entities_size) = 0;

  virtual int get_entities(const int iface_no,
                           const int dimension,
                           iBase_EntitySetHandle set_handle,
                           iBase_EntityHandle **entities,
                           int *entities_allocated,
                           int *entities_size) = 0;

  virtual int get_ents_dims(const int iface_no,
                            iBase_EntityHandle *entities,
                            int entities_size,
                            int **ent_types,
                            int *ent_types_alloc,
                            int *ent_types_size) = 0;

  int set_assoc_tags(iBase_EntityHandle    ent1, iBase_EntityHandle    ent2);
  int set_assoc_tags(iBase_EntitySetHandle set1, iBase_EntityHandle    ent2);
  int set_assoc_tags(iBase_EntityHandle    ent1, iBase_EntitySetHandle set2);
  int set_assoc_tags(iBase_EntitySetHandle set1, iBase_EntitySetHandle set2);

  int get_assoc_tags(const int iface_no, iBase_EntityHandle *entities,
                     int num_entities, iBase_EntityHandle *tag_values);
  int get_assoc_tags(const int iface_no, iBase_EntitySetHandle *sets,
                     int num_sets, iBase_EntityHandle *tag_values);
  int get_assoc_tags(const int iface_no, iBase_EntityHandle *entities,
                     int num_entities, iBase_EntitySetHandle *tag_values);
  int get_assoc_tags(const int iface_no, iBase_EntitySetHandle *sets,
                     int num_sets, iBase_EntitySetHandle *tag_values);
  int get_assoc_tags(const int iface_no, iBase_EntityHandle *entities,
                     int num_entities, iBase_EntityIterator *tag_values);
  int get_assoc_tags(const int iface_no, iBase_EntitySetHandle *sets,
                     int num_sets, iBase_EntityIterator *tag_values);

  int get_gid_tags(const int iface_no, iBase_EntityHandle *entities,
                   int num_entities, int *tag_values);
  int get_gid_tags(const int iface_no, iBase_EntitySetHandle *sets,
                   int num_sets, int *tag_values);

  int get_dim_tags(const int iface_no, iBase_EntityHandle *entities,
                   int num_entities, int *tag_values);
  int get_dim_tags(const int iface_no, iBase_EntitySetHandle *sets,
                   int num_sets, int *tag_values);
protected:
  AssocPair(iRel_Instance instance, 
            iRel_RelationType ent_or_set0, iRel_IfaceType type0,
            iRel_RelationType ent_or_set1, iRel_IfaceType type1);

  virtual ~AssocPair();

  int create_tags();
  int destroy_tags();

  virtual int get_tags(const int iface_no,
                       iBase_EntityHandle *entities,
                       const int num_entities,
                       iBase_TagHandle tag_handle,
                       void *tag_values,
                       const int tag_size) = 0;

  virtual int get_tags(const int iface_no,
                       iBase_EntitySetHandle *sets,
                       const int num_sets,
                       iBase_TagHandle tag_handle,
                       void *tag_values,
                       const int tag_size) = 0;

  virtual int set_tags(const int iface_no,
                       iBase_EntityHandle *entities,
                       const int num_entities,
                       iBase_TagHandle tag_handle,
                       const void *tag_values,
                       const int tag_size) = 0;

  virtual int set_tags(const int iface_no,
                       iBase_EntitySetHandle *sets,
                       const int num_sets,
                       iBase_TagHandle tag_handle,
                       const void *tag_values,
                       const int tag_size) = 0;

  virtual int rmv_tags(const int iface_no,
                       iBase_EntityHandle *entities,
                       const int num_entities,
                       iBase_TagHandle tag_handle) = 0;

  virtual int rmv_tags(const int iface_no,
                       iBase_EntitySetHandle *sets,
                       const int num_sets,
                       iBase_TagHandle tag_handle) = 0;

  virtual int get_iterator(const int iface_no,
                           iBase_EntitySetHandle set,
                           iBase_EntityIterator *iter) = 0;

  virtual int tag_get_handle(const int iface_no,
                             const char *tag_name,
                             const int tag_size_bytes,
                             const int tag_data_type,
                             const bool create_if_missing,
                             void *default_val,
                             iBase_TagHandle *tag_handle) = 0;

  virtual int tag_destroy(const int iface_no, iBase_TagHandle tag_handle) = 0;

  iRel_Instance instance;
private:
  AssocPair();

  int pairId;

  iRel_RelationType entOrSet[2];
  iRel_IfaceType ifaceTypes[2];
  iBase_TagHandle assocTags[2], gidTags[2], dimTags[2];

  static const char *GLOBAL_ID_TAG_NAME;
  static const char *GEOM_DIMENSION_TAG_NAME;
  static const char *ASSOCIATION_TAG_NAME;
  static int currId;
};

static inline AssocPair *assocpair_handle(iRel_PairHandle pair)
{
  return reinterpret_cast<AssocPair*>(pair);
}
#define ASSOCPAIRI assocpair_handle(pair)


#endif
