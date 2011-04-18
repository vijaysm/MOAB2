#ifndef FBENGINE_HPP_
#define FBENGINE_HPP_
#include <stdlib.h>

#include <vector>
#include <map>

#include "moab/Types.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/CartVect.hpp"

namespace moab {
class GeomTopoTool;

// some forward declarations
class SmoothFace;
class SmoothCurve;

/*
 *  Facet Based engine class for mesh-based geometry
 */
class FBEngine {
public:
  FBEngine(Interface *impl, GeomTopoTool* geomTopoTool = NULL,
      const bool smooth = false);

  ~FBEngine();

  ErrorCode Init();

  ErrorCode getRootSet(EntityHandle * root_set);

  ErrorCode getNumEntSets(EntityHandle set, int num_hops, int * all_sets);

  ErrorCode createEntSet(int isList, EntityHandle * pSet);

  ErrorCode addEntSet(EntityHandle entity_set_to_add,
      EntityHandle entity_set_handle);

  ErrorCode getEntities(EntityHandle root_set, int ent_type, Range & gentities);

  ErrorCode addEntArrToSet(Range entities, EntityHandle set);

  ErrorCode getNumOfType(EntityHandle set, int ent_type, int * pNum);

  ErrorCode getEntType(EntityHandle gent, int * type);

  ErrorCode getEntBoundBox(EntityHandle this_gent, double * x0, double * y0,
      double * z0, double * x1, double *y1, double * z1);
  ErrorCode getEntClosestPt(EntityHandle this_gent, double x, double y,
      double z, double * x1, double * y1, double *y3);

  ErrorCode getVtxCoord(EntityHandle this_gent, double * x0, double * y0, double * z0);

  ErrorCode gsubtract(EntityHandle entity_set_1, EntityHandle entity_set_2,
      EntityHandle result_entity_set);

  ErrorCode getEntNrmlXYZ( EntityHandle entity_handle, double x, double y, double z,
      double* nrml_i, double* nrml_j, double* nrml_k);

  ErrorCode getPntRayIntsct( double x, double y, double z,
      double dir_x, double dir_y, double dir_z,
      std::vector<EntityHandle> &intersect_entity_handles,
      /* int storage_order,*/
      std::vector<double> & intersect_coords,
      std::vector<double> & param_coords);

  // some new methods, that are needed

  ErrorCode createTag( const char* tag_name,
                    int tag_num_type_values,
                    int tag_type,
                    Tag & tag_handle_out );

  Interface * moab_instance () { return _mbImpl; }

  ErrorCode getArrData(const moab::EntityHandle* entity_handles,
      int entity_handles_size,
      Tag tag_handle,
      void* tag_values_out);

  ErrorCode setArrData(const EntityHandle* entity_handles,
        int entity_handles_size,
        Tag tag_handle,
        const void* tag_values);

  ErrorCode getEntAdj(EntityHandle handle,
        int type_requested, Range & adjEnts  );

  ErrorCode getEgFcSense(EntityHandle mbedge, EntityHandle mbface, int & sense );

  ErrorCode measure(const EntityHandle * moab_entities, int entities_size,
      double * measures);

  // to do
  ErrorCode getEntNrmlSense( EntityHandle face, EntityHandle region,
      int& sense );

  ErrorCode getEgEvalXYZ( EntityHandle edge,
                                 double x, double y, double z,
                                 double& on_x, double& on_y, double& on_z,
                                 double& tngt_i, double& tngt_j, double& tngt_k,
                                 double& cvtr_i, double& cvtr_j, double& cvtr_k );
  ErrorCode getFcEvalXYZ( EntityHandle face,
                                 double x, double y, double z,
                                 double& on_x, double& on_y, double& on_z,
                                 double& nrml_i, double& nrml_j, double& nrml_k,
                                 double& cvtr1_i, double& cvtr1_j, double& cvtr1_k,
                                 double& cvtr2_i, double& cvtr2_j, double& cvtr2_k );

  ErrorCode getEgVtxSense( EntityHandle edge, EntityHandle vtx1, EntityHandle vtx2, int& sense );

  ErrorCode getEntURange( EntityHandle edge,
                                 double& u_min, double& u_max );

  ErrorCode getEntUtoXYZ( EntityHandle edge, double u,
                                 double& x, double& y, double& z );
  ErrorCode isEntAdj( EntityHandle entity1, EntityHandle entity2,
      bool& adjacent_out );

  ErrorCode split_surface_with_direction(EntityHandle face, std::vector<double> & xyz, double * direction,
      EntityHandle & newFace);
  // these new points will be on edges or triangles, if in interior of triangles
  ErrorCode split_surface(EntityHandle face, std::vector<double> & points,
      std::vector<EntityHandle> & entities, std::vector<EntityHandle> & triangles,
      EntityHandle & newFace);

  // helper for cleaning the stuff
  // will be called if the topology is modified
  void clean();

  void delete_smooth_tags();
private:

  ErrorCode initializeSmoothing();

  ErrorCode getAdjacentEntities(const EntityHandle from, const int to_dim,
      Range &adj_ents);

  ErrorCode compute_intersection_points(EntityHandle & face, CartVect & p1, CartVect & p2,
      EntityHandle from, EntityHandle to, CartVect & Dir, std::vector<double> & points,
      std::vector<EntityHandle> & entities, std::vector<EntityHandle> & triangles);

  ErrorCode BreakTriangle(EntityHandle tri, EntityHandle e1, EntityHandle e3, EntityHandle n1,
      EntityHandle n2, EntityHandle n3, Range & new_triangles);// nodesAlongPolyline are on entities!

  ErrorCode BreakTriangle2(EntityHandle tri, EntityHandle e1, EntityHandle e2, EntityHandle n1,
        EntityHandle n2, Range & new_triangles);// nodesAlongPolyline are on entities!

  void print_debug_triangle(EntityHandle triangle);

  ErrorCode  create_new_gedge(std::vector<EntityHandle> &nodesAlongPolyline,
      EntityHandle & new_geo_edge, Range & geo_vertices);

  // used for splitting surfaces
  ErrorCode separate (Range & iniTriangles, Range & triToDelete,
      Range & new_triangles, EntityHandle new_geo_edge, Range & first,
      Range & second);

  ErrorCode smooth_new_intx_points(EntityHandle face,
      std::vector<EntityHandle> & nodesAlongPolyline);

  Interface * _mbImpl;

  GeomTopoTool* _my_geomTopoTool;
  bool _t_created;
  bool _smooth;
  bool _initialized;
  Range _my_gsets[4];
  // these are used only for smooth evaluations
  // these smooth faces and edges will be initialized after reading the file
  // the maps keep the link between EH in moab (geom sets) and
  //   their corresponding smooth counterparts
  std::map<EntityHandle, SmoothFace*> _faces;
  std::map<EntityHandle, SmoothCurve*> _edges;
  SmoothFace ** _smthFace;
  SmoothCurve ** _smthCurve;

};

} // namespace moab
#endif /* FBENGINE_HPP_ */
