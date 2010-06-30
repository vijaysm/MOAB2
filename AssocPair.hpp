#ifndef ASSOCPAIR_HPP
#define ASSOCPAIR_HPP

#include "iRel.h"
#include <sstream>

class Lasso;

class AssocPair 
{
public:
  friend class Lasso;
  
  IfaceType iface_type(const int iface_no);
  RelationType ent_or_set(const int iface_no);

  virtual iBase_Instance iface_instance(const int iface_no) = 0;
  
  bool equivalent(IfaceType type1, IfaceType type2,
                  bool *order_switched = NULL);

  virtual bool equivalent(iBase_Instance iface1, iBase_Instance iface2,
                          bool *order_switched = NULL) = 0;
  
  virtual bool contains(iBase_Instance iface) = 0;
  
  virtual bool same_interface(iBase_Instance iface1, iBase_Instance iface2) = 0;
  
  virtual int get_int_tags(const int iface_no,
                                       iBase_EntityHandle *entities,
                                       int num_entities,
                                       iBase_TagHandle tag_handle,
                                       int *tag_values) = 0;
  
  virtual int get_eh_tags(const int iface_no,
                                      iBase_EntityHandle *entities,
                                      int num_entities,
                                      iBase_TagHandle tag_handle,
                                      iBase_EntityHandle *tag_values) = 0;
  
  virtual int set_int_tags(const int iface_no,
                                       iBase_EntityHandle *entities,
                                       int num_entities,
                                       iBase_TagHandle tag_handle,
                                       int *tag_values) = 0;
  
  virtual int set_eh_tags(const int iface_no,
                                      iBase_EntityHandle *entities,
                                      int num_entities,
                                      iBase_TagHandle tag_handle,
                                      iBase_EntityHandle *tag_values) = 0;
  
  virtual int get_int_tags(const int iface_no,
                                       iBase_EntitySetHandle *entities,
                                       int num_entities,
                                       iBase_TagHandle tag_handle,
                                       int *tag_values) = 0;
  
  virtual int get_eh_tags(const int iface_no,
                                      iBase_EntitySetHandle *entities,
                                      int num_entities,
                                      iBase_TagHandle tag_handle,
                                      iBase_EntityHandle *tag_values) = 0;
  
  virtual int set_int_tags(const int iface_no,
                                       iBase_EntitySetHandle *entities,
                                       int num_entities,
                                       iBase_TagHandle tag_handle,
                                       int *tag_values) = 0;
  
  virtual int set_eh_tags(const int iface_no,
                                      iBase_EntitySetHandle *entities,
                                      int num_entities,
                                      iBase_TagHandle tag_handle,
                                      iBase_EntityHandle *tag_values) = 0;

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

  virtual iBase_TagHandle tag_get_handle(const int iface_no,
                                         const char *tag_name,
                                         const int tag_size_values,
                                         const int tag_data_type,
                                         const bool create_if_missing,
                                         void *default_val = NULL) = 0;

  virtual int create_mesh_vertex(const double x,
                                             const double y,
                                             const double z,
                                             iBase_EntityHandle &vertex) = 0;
  
  virtual int create_mesh_entity(iMesh_EntityTopology ent_topology,
                                 iBase_EntityHandle *lower_handles,
                                 int lower_handles_size,
                                 iBase_EntityHandle &new_ent,
                                 int &status) = 0;

  virtual int set_mesh_coords(iBase_EntityHandle *verts,
                                          int verts_size,
                                          double *coords,
                                          iBase_StorageOrder order) = 0;
  

  virtual int get_mesh_coords(iBase_EntityHandle *verts,
                                          int verts_size,
                                          double **coords,
                                          int *coords_alloc,
                                          int *coords_size,
                                          iBase_StorageOrder order) = 0;
  
  virtual int get_closest_pt(iBase_EntityHandle *gents, 
                                         int gents_size,
                                         iBase_StorageOrder storage_order,
                                         double *near_coords, 
                                         int near_coords_size,
                                         double **on_coords,
                                         int *on_coords_alloc, 
                                         int *on_coords_size) = 0;
  
  int set_assoc_tags(iBase_EntityHandle    ent1, iBase_EntityHandle    ent2);
  int set_assoc_tags(iBase_EntitySetHandle ent1, iBase_EntityHandle    ent2);
  int set_assoc_tags(iBase_EntityHandle    ent1, iBase_EntitySetHandle ent2);
  int set_assoc_tags(iBase_EntitySetHandle ent1, iBase_EntitySetHandle ent2);
  
  RelationType entOrSet[2];

  iBase_TagHandle assocTags[2], gidTags[2], dimTags[2];

  static const char *GLOBAL_ID_TAG_NAME;
  static const char *GEOM_DIMENSION_TAG_NAME;
  static const char *ASSOCIATION_TAG_NAME;

protected:
  AssocPair(RelationType ent_or_set0, IfaceType type0,
            RelationType ent_or_set1, IfaceType type1,
            Lasso *lasso);
  
  int create_tags();
  
  IfaceType ifaceTypes[2];

  virtual ~AssocPair();

private:
  AssocPair();
  
  Lasso *myLasso;
  
  int pairId;
  
};

#endif
