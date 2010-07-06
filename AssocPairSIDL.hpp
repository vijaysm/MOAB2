#ifndef ASSOCPAIRSIDL_HPP
#define ASSOCPAIRSIDL_HPP

#include "AssocPair.hpp"
#include "sidl.h"
#include "sidl_BaseInterface_IOR.h"

class AssocPairSIDL : public AssocPair 
{
public:
  AssocPairSIDL(::sidl::BaseInterface iface1,
                const int ent_or_set1,
                ::sidl::BaseInterface iface2,
                const int ent_or_set2,
                Lasso *lasso);

  virtual ~AssocPairSIDL();
  
  virtual bool equivalent(iBase_Instance iface1, iBase_Instance iface2,
                          bool *order_switched = NULL);
  
  virtual bool contains(iBase_Instance iface);
  
  virtual bool same_interface(iBase_Instance iface1, iBase_Instance iface2);
  
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

  virtual iBase_TagHandle tag_get_handle(const int iface_no,
                                         const char *tag_name,
                                         const int tag_size_bytes,
                                         const int tag_data_type,
                                         const bool create_if_missing,
                                         void *default_val = NULL);
private:
  static IFaceType get_type(::sidl::BaseInterface iface);

  ::sidl::BaseInterface ifaceInstances[2];
};

#endif
