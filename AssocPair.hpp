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
  
  int ent_or_set(const int iface_no);

  virtual bool equivalent(iBase_Instance iface1, iBase_Instance iface2,
                          bool *order_switched = NULL) = 0;
  
  virtual bool contains(iBase_Instance iface) = 0;
  
  virtual bool same_interface(iBase_Instance iface1, iBase_Instance iface2) = 0;
  
  virtual int get_int_tags(const int iface_no,
                                       iBase_EntityHandle *entities,
                                       int num_entities,
                                       iBase_TagHandle tag_handle,
                                       bool are_sets,
                                       int *tag_values) = 0;
  
  virtual int get_eh_tags(const int iface_no,
                                      iBase_EntityHandle *entities,
                                      int num_entities,
                                      iBase_TagHandle tag_handle,
                                      bool are_sets,
                                      iBase_EntityHandle *tag_values) = 0;
  
  virtual int set_int_tags(const int iface_no,
                                       iBase_EntityHandle *entities,
                                       int num_entities,
                                       iBase_TagHandle tag_handle,
                                       bool are_sets,
                                       int *tag_values) = 0;
  
  virtual int set_eh_tags(const int iface_no,
                                      iBase_EntityHandle *entities,
                                      int num_entities,
                                      iBase_TagHandle tag_handle,
                                      bool are_sets,
                                      iBase_EntityHandle *tag_values) = 0;

  virtual int get_all_entities(const int iface_no,
                                           const int dimension,
                                           const bool are_sets,
                                           iBase_EntityHandle **entities,
                                           int *entities_alloc,
                                           int *entities_size) = 0;

  virtual int get_entities(const int iface_no,
                                       const int dimension,
                                       iBase_EntityHandle set_handle,
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
                                          int *order) = 0;
  
  virtual int get_closest_pt(iBase_EntityHandle *gents, 
                                         int gents_size,
                                         iBase_StorageOrder storage_order,
                                         double *near_coords, 
                                         int near_coords_size,
                                         double **on_coords,
                                         int *on_coords_alloc, 
                                         int *on_coords_size) = 0;
  
  int set_assoc_tags(iBase_EntityHandle ent1,
                                 int is_set1,
                                 iBase_EntityHandle ent2,
                                 int is_set2);
  
  int entOrSet[2];

  iBase_TagHandle assocTags[2], gidTags[2], dimTags[2];

  static char *GLOBAL_ID_TAG_NAME;
  static char *GEOM_DIMENSION_TAG_NAME;
  static char *ASSOCIATION_TAG_NAME;

protected:
  AssocPair(const int ent_or_set0, const int ent_or_set1,
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
