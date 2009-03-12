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

  virtual bool equivalent(iBase_Instance iface1, iBase_Instance iface2,
                          bool *order_switched = NULL);
  
  virtual bool contains(iBase_Instance iface);
  
  virtual bool same_interface(iBase_Instance iface1, iBase_Instance iface2);
  
  virtual int get_int_tags(const int iface_no,
                                       iBase_EntityHandle *entities,
                                       int num_entities,
                                       iBase_TagHandle tag_handle,
                                       int *tag_values);
  
  virtual int get_eh_tags(const int iface_no,
                                      iBase_EntityHandle *entities,
                                      int num_entities,
                                      iBase_TagHandle tag_handle,
                                      iBase_EntityHandle *tag_values);
  
  virtual int set_int_tags(const int iface_no,
                                       iBase_EntityHandle *entities,
                                       int num_entities,
                                       iBase_TagHandle tag_handle,
                                       int *tag_values);
  
  virtual int set_eh_tags(const int iface_no,
                                      iBase_EntityHandle *entities,
                                      int num_entities,
                                      iBase_TagHandle tag_handle,
                                      iBase_EntityHandle *tag_values);
  
  virtual int get_int_tags(const int iface_no,
                                       iBase_EntitySetHandle *entities,
                                       int num_entities,
                                       iBase_TagHandle tag_handle,
                                       int *tag_values);
  
  virtual int get_eh_tags(const int iface_no,
                                      iBase_EntitySetHandle *entities,
                                      int num_entities,
                                      iBase_TagHandle tag_handle,
                                      iBase_EntityHandle *tag_values);
  
  virtual int set_int_tags(const int iface_no,
                                       iBase_EntitySetHandle *entities,
                                       int num_entities,
                                       iBase_TagHandle tag_handle,
                                       int *tag_values);
  
  virtual int set_eh_tags(const int iface_no,
                                      iBase_EntitySetHandle *entities,
                                      int num_entities,
                                      iBase_TagHandle tag_handle,
                                      iBase_EntityHandle *tag_values);

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

  virtual iBase_TagHandle tag_get_handle(const int iface_no,
                                         const char *tag_name,
                                         const int tag_size_bytes,
                                         const int tag_data_type,
                                         const bool create_if_missing,
                                         void *default_val = NULL);
  
  virtual int create_mesh_vertex(const double x,
                                             const double y,
                                             const double z,
                                             iBase_EntityHandle &vertex);
  
  virtual int create_mesh_entity(iMesh_EntityTopology ent_topology,
                                             iBase_EntityHandle *lower_handles,
                                             int lower_handles_size,
                                             iBase_EntityHandle &new_ent,
                                 int &status);

  virtual int set_mesh_coords(iBase_EntityHandle *verts,
                                          int verts_size,
                                          double *coords,
                                          iBase_StorageOrder order);
  

  virtual int get_mesh_coords(iBase_EntityHandle *verts,
                                          int verts_size,
                                          double **coords,
                                          int *coords_alloc,
                                          int *coords_size,
                                          iBase_StorageOrder order);
  
  virtual int get_closest_pt(iBase_EntityHandle *gents, 
                                         int gents_size,
                                         iBase_StorageOrder storage_order,
                                         double *near_coords, 
                                         int near_coords_size,
                                         double **on_coords,
                                         int *on_coords_alloc, 
                                         int *on_coords_size);
  
  virtual ~AssocPairSIDL();
  
private:

  bool is_type(::sidl::BaseInterface iface,
               IfaceType this_type);

  ::sidl::BaseInterface ifaceInstances[2];
  
};

#endif
