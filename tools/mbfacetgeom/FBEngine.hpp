#ifndef FBENGINE_HPP_
#define FBENGINE_HPP_
#include <stdlib.h>

#include <vector>
#include <map>

#include "moab/Types.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"

namespace moab {
class GeomTopoTool;
}

using namespace moab;

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

private:

  ErrorCode initializeSmoothing();

  ErrorCode getAdjacentEntities(const EntityHandle from, const int to_dim,
      Range &adj_ents);

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

#endif /* FBENGINE_HPP_ */
