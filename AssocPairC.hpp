#ifndef ASSOCPAIRC_HPP
#define ASSOCPAIRC_HPP

#include "AssocPair.hpp"

class AssocPairC : public AssocPair
{
public:
  AssocPairC(iBase_Instance iface0,
             RelationType ent_or_set0,
             IfaceType type0,
             iBase_Instance iface1,
             RelationType ent_or_set1,
             IfaceType type1);

  virtual ~AssocPairC();

  virtual iBase_Instance iface_instance(const int iface_no);

  virtual bool equivalent(iBase_Instance iface0, iBase_Instance iface1,
                          bool *order_switched = NULL);

  virtual bool contains(iBase_Instance iface);

  virtual bool same_interface(iBase_Instance iface0, iBase_Instance iface1);

  virtual int get_all_entities(const int iface_no,
                               const int dimension,
                               iBase_EntityHandle **entities,
                               int *entities_alloc,
                               int *entities_size);

  virtual int get_all_sets(const int iface_no,
                               iBase_EntitySetHandle **entities,
                               int *entities_alloc,
                               int *entities_size);

  virtual int get_entities(const int iface_no,
                           const int dimension,
                           iBase_EntitySetHandle set_handle,
                           iBase_EntityHandle **entities,
                           int *entities_allocated,
                           int *entities_size);

  virtual int get_ents_dims(const int iface_no,
                            iBase_EntityHandle *entities,
                            int entities_size,
                            int **ent_types,
                            int *ent_types_alloc,
                            int *ent_types_size);
protected:
  virtual int get_tags(const int iface_no,
                       iBase_EntityHandle *entities,
                       const int num_entities,
                       iBase_TagHandle tag_handle,
                       void *tag_values,
                       const int tag_size);

  virtual int get_tags(const int iface_no,
                       iBase_EntitySetHandle *sets,
                       const int num_sets,
                       iBase_TagHandle tag_handle,
                       void *tag_values,
                       const int tag_size);

  virtual int set_tags(const int iface_no,
                       iBase_EntityHandle *entities,
                       const int num_entities,
                       iBase_TagHandle tag_handle,
                       const void *tag_values,
                       const int tag_size);

  virtual int set_tags(const int iface_no,
                       iBase_EntitySetHandle *sets,
                       const int num_sets,
                       iBase_TagHandle tag_handle,
                       const void *tag_values,
                       const int tag_size);

  virtual int get_iterator(const int iface_no,
                           iBase_EntitySetHandle set,
                           iBase_EntityIterator *iter);

  virtual int tag_get_handle(const int iface_no,
                             const char *tag_name,
                             const int tag_size_bytes,
                             const int tag_data_type,
                             const bool create_if_missing,
                             void *default_val,
                             iBase_TagHandle *tag_handle);

  virtual int tag_destroy(const int iface_no, iBase_TagHandle tag_handle);
private:
  iBase_Instance ifaceInstances[2];
};

#endif
