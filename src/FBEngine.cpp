#include <iostream>
#include <map>

#include "moab/FBEngine.hpp"
#include "moab/Interface.hpp"
#include "moab/GeomTopoTool.hpp"
#include "moab/OrientedBoxTreeTool.hpp"

#include <stdlib.h>
#include <cstring>
#include <map>
#include <set>
#include <queue>
#include "assert.h"

#include "SmoothCurve.hpp"
#include "SmoothFace.hpp"

// this is just to replace MBI with moab interface, which is _mbImpl in this class
#define MBI _mbImpl
#define MBERRORR(rval, STR) { if (MB_SUCCESS != rval) { std::cout<<STR<<std::endl; return rval; } }

namespace moab {

// some tolerances for ray tracing and geometry intersections
// these are involved in ray tracing, at least

unsigned min_tolerace_intersections = 1000;
double tolerance = 0.01; // TODO: how is this used ????
double tolerance_segment = 0.000001; // for segments intersection, points collapse
const bool Debug_surf_eval = false;
bool debug_splits = false;

// will compute intersection between a segment and slice of a plane
// output is the intersection point
bool intersect_segment_and_plane_slice(CartVect & from, CartVect & to,
    CartVect & p1, CartVect & p2, CartVect & Dir, CartVect & normPlane,
    CartVect & intx_point, double & parPos)
{
  //
  // plane eq is normPlane % r + d = 0, or normPlane % r - normPlane%p1 = 0
  double dd = -normPlane % p1;
  double valFrom = normPlane % from + dd;
  double valTo = normPlane % to + dd;

  if (fabs(valFrom) < tolerance_segment) {
    intx_point = from;
    parPos = 0.;
    double proj1 = (intx_point - p1) % (p2 - p1);
    double proj2 = (intx_point - p2) % (p1 - p2);
    if (proj1 <= -tolerance_segment || proj2 <= -tolerance_segment)
      return false;
    if (debug_splits)
      std::cout << "intx : " << intx_point << "\n";
    return true;
  }
  if (fabs(valTo) < tolerance_segment) {
    intx_point = to;
    parPos = 1;
    double proj1 = (intx_point - p1) % (p2 - p1);
    double proj2 = (intx_point - p2) % (p1 - p2);
    if (proj1 <= -tolerance_segment || proj2 <= -tolerance_segment)
      return false;
    if (debug_splits)
      std::cout << "intx : " << intx_point << "\n";
    return true;
  }
  if (valFrom * valTo > 0)
    return false; // no intersection, although it could be very close
  // else, it could intersect the plane; check for the slice too.
  parPos = valFrom / (valFrom - valTo);// this is 0 for valFrom 0, 1 for valTo 0
  intx_point = from + (to - from) * parPos;
  // now check if the intx_point is indeed between p1 and p2 in the slice.
  double proj1 = (intx_point - p1) % (p2 - p1);
  double proj2 = (intx_point - p2) % (p1 - p2);
  if (proj1 <= -tolerance_segment || proj2 <= -tolerance_segment)
    return false;

  if (debug_splits)
    std::cout << "intx : " << intx_point << "\n";
  return true;
}

ErrorCode area_coordinates(Interface * mbi, EntityHandle tri, CartVect & pnt,
    double * area_coord, EntityHandle & boundary_handle, bool & onBoundary)
{

  int nnodes;
  const EntityHandle * conn3;
  ErrorCode rval = mbi->get_connectivity(tri, conn3, nnodes);
  MBERRORR(rval, "Failed to get connectivity");
  assert(3 == nnodes);
  CartVect P[3];
  rval = mbi->get_coords(conn3, nnodes, (double*) &P[0]);
  MBERRORR(rval, "Failed to get coordinates");

  CartVect r0(P[0] - pnt);
  CartVect r1(P[1] - pnt);
  CartVect r2(P[2] - pnt);
  if (r0.length() < tolerance_segment) {
    area_coord[0] = 1.;
    area_coord[1] = 0.;
    area_coord[2] = 0.;
    boundary_handle = conn3[0];
    onBoundary = true;
    return MB_SUCCESS;
  }
  if (r1.length() < tolerance_segment) {
    area_coord[0] = 0.;
    area_coord[1] = 1.;
    area_coord[2] = 0.;
    boundary_handle = conn3[1];
    onBoundary = true;
    return MB_SUCCESS;
  }
  if (r1.length() < tolerance_segment) {
    area_coord[0] = 0.;
    area_coord[1] = 0.;
    area_coord[2] = 1.;
    boundary_handle = conn3[2];
    onBoundary = true;
    return MB_SUCCESS;
  }

  CartVect v1(P[1] - P[0]);
  CartVect v2(P[2] - P[0]);

  double areaDouble = (v1 * v2).length();// the same for CartVect
  if (areaDouble < tolerance_segment * tolerance_segment) {
    MBERRORR(MB_FAILURE, "area of triangle too small");
  }
  area_coord[0] = (r1 * r2).length() / areaDouble;
  area_coord[1] = (r2 * r0).length() / areaDouble;
  area_coord[2] = (r0 * r1).length() / areaDouble;

  if (fabs(area_coord[0] + area_coord[1] + area_coord[2] - 1)
      > tolerance_segment) {
    MBERRORR(MB_FAILURE, "point outside triangle");
  }
  onBoundary = false;
  bool side0 = (area_coord[0] < tolerance_segment);
  bool side1 = (area_coord[1] < tolerance_segment);
  bool side2 = (area_coord[2] < tolerance_segment);
  if (side0 || side1 || side2) {
    onBoundary = true;
  } else
    return MB_SUCCESS; // interior point
  // now, find out what boundary is in question
  // first, get all edges, in order
  std::vector<EntityHandle> edges;
  EntityHandle nn2[2];
  for (int i = 0; i < 3; i++) {
    nn2[0] = conn3[(i + 1) % 3];
    nn2[1] = conn3[(i + 2) % 3];
    std::vector<EntityHandle> adjacent;
    rval = mbi->get_adjacencies(nn2, 2, 1, false, adjacent,
        Interface::INTERSECT);
    MBERRORR(rval, "Failed to get edges");
    if (adjacent.size() != 1)
      MBERRORR(MB_FAILURE, "Failed to get adjacent edges");
    // should be only one edge here
    edges.push_back(adjacent[0]);
  }

  if (side0)
    boundary_handle = edges[0];
  if (side1)
    boundary_handle = edges[1];
  if (side2)
    boundary_handle = edges[2];

  return MB_SUCCESS;
}

FBEngine::FBEngine(Interface *impl, GeomTopoTool * topoTool, const bool smooth) :
  _mbImpl(impl), _my_geomTopoTool(topoTool), _t_created(false),
      _smooth(smooth), _initialized(false), _smthFace(NULL), _smthCurve(NULL)
{
  if (!_my_geomTopoTool) {
    _my_geomTopoTool = new GeomTopoTool(_mbImpl);
    _t_created = true;
  }
  // should this be part of the constructor or not?
  //Init();
}
FBEngine::~FBEngine()
{
  clean();
  _smooth = false;
}

void FBEngine::clean()
{
  if (_smooth) {
    _faces.clear();
    _edges.clear();
    int size1 = _my_gsets[1].size();
    int i = 0;
    for (i = 0; i < size1; i++)
      delete _smthCurve[i];
    delete[] _smthCurve;
    _smthCurve = NULL;
    size1 = _my_gsets[2].size();
    for (i = 0; i < size1; i++)
      delete _smthFace[i];
    delete[] _smthFace;
    _smthFace = NULL;
    //_smooth = false;
  }

  for (int j = 0; j < 4; j++)
    _my_gsets[j].clear();
  if (_t_created)
    delete _my_geomTopoTool;
  _my_geomTopoTool = NULL;
  _t_created = false;
}

ErrorCode FBEngine::Init()
{
  if (!_initialized) {
    if (!_my_geomTopoTool)
      return MB_FAILURE;

    ErrorCode rval = _my_geomTopoTool->find_geomsets(_my_gsets);
    assert(rval == MB_SUCCESS);

    rval = split_quads();
    assert (rval == MB_SUCCESS);

    rval = _my_geomTopoTool->construct_obb_trees();
    assert(rval == MB_SUCCESS);

    if (_smooth)
      rval = initializeSmoothing();
    assert(rval == MB_SUCCESS);

    _initialized = true;
  }
  return MB_SUCCESS;
}
ErrorCode FBEngine::initializeSmoothing()
{
  //
  /*ErrorCode rval = Init();
   MBERRORR(rval, "failed initialize");*/
  // first of all, we need to retrieve all the surfaces from the (root) set
  // in icesheet_test we use iGeom, but maybe that is a stretch
  // get directly the sets with geom dim 2, and from there create the SmoothFace
  Tag geom_tag, gid_tag;
  ErrorCode rval = MBI->tag_get_handle(GEOM_DIMENSION_TAG_NAME, geom_tag);
  MBERRORR(rval, "can't get geom tag");
  rval = MBI->tag_get_handle(GLOBAL_ID_TAG_NAME, gid_tag);
  MBERRORR(rval, "can't get id tag");
  int numSurfaces = _my_gsets[2].size();
  //SmoothFace ** smthFace = new SmoothFace *[numSurfaces];
  _smthFace = new SmoothFace *[numSurfaces];

  // there should also be a map from surfaces to evaluators
  //std::map<MBEntityHandle, SmoothFace*> mapSurfaces;

  int i = 0;
  Range::iterator it;
  for (it = _my_gsets[2].begin(); it != _my_gsets[2].end(); it++, i++) {
    EntityHandle face = *it;
    _smthFace[i] = new SmoothFace(MBI, face, _my_geomTopoTool);// geom topo tool will be used for searching,
    // among other things; also for senses in edge sets...
    _faces[face] = _smthFace[i];
  }

  int numCurves = _my_gsets[1].size();//csets.size();
  //SmoothCurve ** smthCurve = new SmoothCurve *[numCurves];
  _smthCurve = new SmoothCurve *[numCurves];
  // there should also be a map from surfaces to evaluators
  //std::map<MBEntityHandle, SmoothCurve*> mapCurves;

  i = 0;
  for (it = _my_gsets[1].begin(); it != _my_gsets[1].end(); it++, i++) {
    EntityHandle curve = *it;
    _smthCurve[i] = new SmoothCurve(MBI, curve);
    _edges[curve] = _smthCurve[i];
  }

  for (i = 0; i < numSurfaces; i++) {
    _smthFace[i]->init_gradient();// this will also retrieve the triangles in each surface
    _smthFace[i]->compute_tangents_for_each_edge();// this one will consider all edges internal, so the
    // tangents are all in the direction of the edge; a little bit of waste, as we store
    // one tangent for each edge node , even though they are equal here...
    // no loops are considered
  }

  // this will be used to mark boundary edges, so for them the control points are computed earlier
  unsigned char value = 0; // default value is "not used"=0 for the tag
  // unsigned char def_data_bit = 1;// valid by default
  // rval = mb->tag_create("valid", 1, MB_TAG_BIT, validTag, &def_data_bit);
  Tag markTag;
  rval = MBI->tag_create("MARKER", 1, MB_TAG_BIT, markTag, &value); // default value : 0 = not computed yet
  // each feature edge will need to have a way to retrieve at every moment the surfaces it belongs to
  // from edge sets, using the sense tag, we can get faces, and from each face, using the map, we can get
  // the SmoothFace (surface evaluator), that has everything, including the normals!!!
  assert(rval==MB_SUCCESS);

  // create the tag also for control points on the edges
  double defCtrlPoints[9] = { 0., 0., 0., 0., 0., 0., 0., 0., 0. };
  Tag edgeCtrlTag;
  rval = MBI->tag_create("CONTROLEDGE", 9 * sizeof(double), MB_TAG_DENSE,
      edgeCtrlTag, &defCtrlPoints);
  assert(rval == MB_SUCCESS);

  Tag facetCtrlTag;
  double defControls[18] = { 0. };
  rval = MBI->tag_create("CONTROLFACE", 18 * sizeof(double), MB_TAG_DENSE,
      facetCtrlTag, &defControls);
  assert(rval == MB_SUCCESS);

  Tag facetEdgeCtrlTag;
  double defControls2[27] = { 0. }; // corresponding to 9 control points on edges, in order from edge 0, 1, 2 ( 1-2, 2-0, 0-1 )
  rval = MBI->tag_create("CONTROLEDGEFACE", 27 * sizeof(double), MB_TAG_DENSE,
      facetEdgeCtrlTag, &defControls2);
  assert(rval == MB_SUCCESS);
  // if the
  double min_dot = -1.0; // depends on _angle, but now we ignore it, for the time being
  for (i = 0; i < numCurves; i++) {
    _smthCurve[i]->compute_tangents_for_each_edge();// do we need surfaces now? or just the chains?
    // the computed edges will be marked accordingly; later one, only internal edges to surfaces are left
    _smthCurve[i]->compute_control_points_on_boundary_edges(min_dot, _faces,
        edgeCtrlTag, markTag);
  }

  // when done with boundary edges, compute the control points on all edges in the surfaces

  for (i = 0; i < numSurfaces; i++) {
    // also pass the tags for
    _smthFace[i]->compute_control_points_on_edges(min_dot, edgeCtrlTag, markTag);
  }

  // now we should be able to compute the control points for the facets

  for (i = 0; i < numSurfaces; i++) {
    // also pass the tags for edge and facet control points
    _smthFace[i]->compute_internal_control_points_on_facets(min_dot,
        facetCtrlTag, facetEdgeCtrlTag);
  }
  // we will need to compute the tangents for all edges in the model
  // they will be needed for control points for each edge
  // the boundary edges and the feature edges are more complicated
  // the boundary edges need to consider their loops, but feature edges need to consider loops and the normals
  // on each connected surface

  // some control points
  if (Debug_surf_eval)
    for (i = 0; i < numSurfaces; i++)
      _smthFace[i]->DumpModelControlPoints();

  return MB_SUCCESS;
}

// clean up the smooth tags data if created, so the files will be smaller
// if saved
// also, recompute the tags if topology is modified
void FBEngine::delete_smooth_tags()
{
  // get all tags from database that are created for smooth data, and
  // delete them; it will delete all data associated with them
  // first tags from faces, edges:
  std::vector<Tag> smoothTags;
  int size1 = (int)_my_gsets[2].size();

  for (int i=0; i<size1; i++)
  {
    // these 2 will append gradient tag and plane tag
    _smthFace[i]->append_smooth_tags(smoothTags);
  }
  // then , get other tags:
  // "TANGENTS", "MARKER", "CONTROLEDGE", "CONTROLFACE", "CONTROLEDGEFACE"
  std::string t1("TANGENTS");
  Tag tag_handle;
  ErrorCode rval = _mbImpl->tag_get_handle( t1.c_str(), tag_handle );
  if (rval != MB_TAG_NOT_FOUND)
    smoothTags.push_back(tag_handle);

  std::string t2("MARKER");
  rval = _mbImpl->tag_get_handle( t2.c_str(), tag_handle );
  if (rval != MB_TAG_NOT_FOUND)
    smoothTags.push_back(tag_handle);

  std::string t3("CONTROLEDGE");
  rval = _mbImpl->tag_get_handle( t3.c_str(), tag_handle );
  if (rval != MB_TAG_NOT_FOUND)
    smoothTags.push_back(tag_handle);

  std::string t4("CONTROLFACE");
  rval = _mbImpl->tag_get_handle( t4.c_str(), tag_handle );
  if (rval != MB_TAG_NOT_FOUND)
    smoothTags.push_back(tag_handle);

  std::string t5("CONTROLEDGEFACE");
  rval = _mbImpl->tag_get_handle( t5.c_str(), tag_handle );
  if (rval != MB_TAG_NOT_FOUND)
    smoothTags.push_back(tag_handle);

  // a lot of tags, delete them
  for (unsigned int k = 0; k<smoothTags.size(); k++ )
  {
    // could be a lot of data
    _mbImpl->tag_delete(smoothTags[k]);
  }
}
#define COPY_RANGE(r, vec) {                      \
    EntityHandle *tmp_ptr = reinterpret_cast<EntityHandle*>(vec);	\
    std::copy(r.begin(), r.end(), tmp_ptr);}

/*static inline void
 ProcessError(const char* desc);*/

ErrorCode FBEngine::getRootSet(EntityHandle * root_set)
{
  *root_set =  _my_geomTopoTool-> get_root_model_set();
  return MB_SUCCESS;
}

ErrorCode FBEngine::getNumEntSets(EntityHandle set, int num_hops,
    int * all_sets)
{
  ErrorCode rval = MBI->num_contained_meshsets(set, all_sets, num_hops + 1);
  return rval;
}

ErrorCode FBEngine::createEntSet(int isList, EntityHandle * pSet)
{
  ErrorCode rval;

  if (isList)
    rval = MBI->create_meshset(MESHSET_ORDERED, *pSet);
  else
    rval = MBI->create_meshset(MESHSET_SET, *pSet);

  return rval;
}

ErrorCode FBEngine::getEntities(EntityHandle set_handle, int entity_type,
    Range & gentities)
{
  int i;
  if (0 > entity_type || 4 < entity_type) {
    return MB_FAILURE;
  } else if (entity_type < 4) {// 4 means all entities
    gentities = _my_gsets[entity_type];// all from root set!
  } else {
    gentities.clear();
    for (i = 0; i < 4; i++) {
      gentities.merge(_my_gsets[i]);
    }
  }
  Range sets;
  // see now if they are in the set passed as input or not
  ErrorCode rval = MBI->get_entities_by_type(set_handle, MBENTITYSET, sets);
  MBERRORR(rval, "can't get sets in the initial set");
  gentities = intersect(gentities, sets);

  return MB_SUCCESS;
}

ErrorCode FBEngine::addEntArrToSet(Range entities, EntityHandle set)
{
  return MBI->add_entities(set, entities);
}

ErrorCode FBEngine::addEntSet(EntityHandle entity_set_to_add,
    EntityHandle entity_set_handle)
{
  return MBI->add_entities(entity_set_handle, &entity_set_to_add, 1);
}

ErrorCode FBEngine::getNumOfType(EntityHandle set, int ent_type, int * pNum)
{
  if (0 > ent_type || 3 < ent_type) {
    std::cout << "Invalid type\n";
    return MB_FAILURE;
  }
  int num_sets;
  ErrorCode rval = MBI->get_number_entities_by_type(set, MBENTITYSET, num_sets);

  MBERRORR(rval, "Failed to get number of sets in the original set.");

  // get also all sets in the set
  Range sets;
  rval = _mbImpl->get_entities_by_type(set, MBENTITYSET, sets, false); // nonrecursive

  // see how many are in the range
  sets = intersect(sets, _my_gsets[ent_type]);
  *pNum = sets.size();
  // we do not really check if it is in the set or not;
  // _my_gsets[i].find(gent) != _my_gsets[i].end()
  return MB_SUCCESS;
}

ErrorCode FBEngine::getEntType(EntityHandle gent, int * type)
{
  for (int i = 0; i < 4; i++) {
    if (_my_gsets[i].find(gent) != _my_gsets[i].end()) {
      *type = i;
      return MB_SUCCESS;
    }
  }
  *type = -1; // failure
  return MB_FAILURE;
}
ErrorCode FBEngine::getEntBoundBox(EntityHandle gent, double* min_x,
    double* min_y, double* min_z, double* max_x, double* max_y, double* max_z)
{
  ErrorCode rval;
  int type;
  rval = getEntType(gent, &type);
  MBERRORR(rval, "Failed to get entity type.");

  if (type == 0) {
    rval = getVtxCoord(gent, min_x, min_y, min_z);
    MBERRORR(rval, "Failed to get vertex coordinates.");
    max_x = min_x;
    max_y = min_y;
    max_z = min_z;
  } else if (type == 1) {
    rval = MB_FAILURE;
    MBERRORR(rval, "iGeom_getEntBoundBox is not supported for Edge entity type.");
  } else if (type == 2 || type == 3) {

    EntityHandle root;
    CartVect center, axis[3];
    rval = _my_geomTopoTool->get_root(gent, root);
    MBERRORR(rval, "Failed to get tree root in iGeom_getEntBoundBox.");
    rval = _my_geomTopoTool->obb_tree()->box(root, center.array(),
        axis[0].array(), axis[1].array(), axis[2].array());
    MBERRORR(rval, "Failed to get closest point in iGeom_getEntBoundBox.");

    CartVect absv[3];
    for (int i = 0; i < 3; i++) {
      absv[i] = CartVect(fabs(axis[i][0]), fabs(axis[i][1]), fabs(axis[i][2]));
    }
    CartVect min, max;
    min = center - absv[0] - absv[1] - absv[2];
    max = center + absv[0] + absv[1] + absv[2];
    *min_x = min[0];
    *min_y = min[1];
    *min_z = min[2];
    *max_x = max[0];
    *max_y = max[1];
    *max_z = max[2];
  } else
    return MB_FAILURE;

  return MB_SUCCESS;
}
ErrorCode FBEngine::getEntClosestPt(EntityHandle this_gent, double near_x,
    double near_y, double near_z, double* on_x, double* on_y, double* on_z)
{
  ErrorCode rval;
  int type;
  rval = getEntType(this_gent, &type);
  MBERRORR(rval, "Failed to get entity type.");

  if (type == 0) {
    rval = getVtxCoord(this_gent, on_x, on_y, on_z);
    MBERRORR(rval, "Failed to get vertex coordinates.");
  } else if (_smooth && type == 1) {
    *on_x = near_x;
    *on_y = near_y;
    *on_z = near_z;
    SmoothCurve * smthcurve = _edges[this_gent];
    // call the new method from smooth edge
    smthcurve->move_to_curve( *on_x, *on_y, *on_z);

  } else if (type == 2 || type == 3) {
    double point[3] = { near_x, near_y, near_z };
    double point_out[3];
    EntityHandle root, facet_out;
    if (_smooth && 2 == type) {
      SmoothFace* smthFace = _faces[this_gent];
      *on_x = near_x;
      *on_y = near_y;
      *on_z = near_z;
      smthFace->move_to_surface(*on_x, *on_y, *on_z);
    } else {
      rval = _my_geomTopoTool->get_root(this_gent, root);
      MBERRORR(rval, "Failed to get tree root in iGeom_getEntClosestPt.");
      rval = _my_geomTopoTool->obb_tree()->closest_to_location(point, root,
          point_out, facet_out);
      MBERRORR(rval, "Failed to get closest point in iGeom_getEntClosestPt.");

      *on_x = point_out[0];
      *on_y = point_out[1];
      *on_z = point_out[2];
    }
  } else
    return MB_TYPE_OUT_OF_RANGE;

  return MB_SUCCESS;
}

ErrorCode FBEngine::getVtxCoord(EntityHandle vertex_handle, double * x0,
    double * y0, double * z0)
{
  int type;
  ErrorCode rval = getEntType(vertex_handle, &type);
  MBERRORR(rval, "Failed to get entity type in getVtxCoord.");

  if (type != 0) {
    rval = MB_FAILURE;
    MBERRORR(rval, "Entity is not a vertex type.");
  }

  Range entities;
  rval = MBI->get_entities_by_type(vertex_handle, MBVERTEX, entities);
  MBERRORR(rval, "can't get nodes in vertex set.");

  if (entities.size() != 1) {
    MBERRORR(MB_FAILURE, "Vertex has multiple points.");
  }
  double coords[3];
  EntityHandle node = entities[0];
  rval = MBI->get_coords(&node, 1, coords);
  MBERRORR(rval, "can't get coordinates.");
  *x0 = coords[0];
  *y0 = coords[1];
  *z0 = coords[2];

  return MB_SUCCESS;
}

ErrorCode FBEngine::gsubtract(EntityHandle entity_set_1,
    EntityHandle entity_set_2, EntityHandle result_entity_set)
{
  /*result_entity_set = subtract(entity_set_1, entity_set_2);*/
  Range ents1, ents2;
  ErrorCode rval = MBI->get_entities_by_type(entity_set_1, MBENTITYSET, ents1);
  MBERRORR(rval, "can't get entities from set 1.");

  rval = MBI->get_entities_by_type(entity_set_2, MBENTITYSET, ents2);
  MBERRORR(rval, "can't get entities from set 2.");

  ents1 = subtract(ents1, ents2);
  rval = MBI->clear_meshset(&result_entity_set, 1);
  MBERRORR(rval, "can't empty set.");

  rval = MBI->add_entities(result_entity_set, ents1);
  MBERRORR(rval, "can't add result to set.");

  return rval;
}

ErrorCode FBEngine::getEntNrmlXYZ(EntityHandle entity_handle, double x,
    double y, double z, double* nrml_i, double* nrml_j, double* nrml_k)
{
  // just do for surface and volume
  int type;
  ErrorCode rval = getEntType(entity_handle, &type);
  MBERRORR(rval, "Failed to get entity type in iGeom_getEntNrmlXYZ.");

  if (type != 2 && type != 3) {
    MBERRORR(MB_FAILURE, "Entities passed into gentityNormal must be face or volume.");
  }

  if (_smooth && 2 == type) {
    SmoothFace* smthFace = _faces[entity_handle];
    //*on_x = near_x; *on_y = near_y; *on_z = near_z;
    smthFace-> normal_at(x, y, z, *nrml_i, *nrml_j, *nrml_k);

  } else {
    // get closest location and facet
    double point[3] = { x, y, z };
    double point_out[3];
    EntityHandle root, facet_out;
    _my_geomTopoTool->get_root(entity_handle, root);
    rval = _my_geomTopoTool->obb_tree()->closest_to_location(point, root,
        point_out, facet_out);
    MBERRORR(rval , "Failed to get closest location in iGeom_getEntNrmlXYZ.");

    // get facet normal
    const EntityHandle* conn;
    int len;
    CartVect coords[3], normal;
    rval = MBI->get_connectivity(facet_out, conn, len);
    MBERRORR(rval, "Failed to get triangle connectivity in iGeom_getEntNrmlXYZ.");
    if (len != 3)
      MBERRORR(MB_FAILURE, " not a triangle, error ");

    rval = MBI->get_coords(conn, len, coords[0].array());
    MBERRORR(rval, "Failed to get triangle coordinates in iGeom_getEntNrmlXYZ.");

    coords[1] -= coords[0];
    coords[2] -= coords[0];
    normal = coords[1] * coords[2];
    normal.normalize();
    *nrml_i = normal[0];
    *nrml_j = normal[1];
    *nrml_k = normal[2];
  }
  return MB_SUCCESS;
}

ErrorCode FBEngine::getPntRayIntsct(double x, double y, double z, double dir_x,
    double dir_y, double dir_z,
    std::vector<EntityHandle> &intersect_entity_handles,
    /* int storage_order,*/
    std::vector<double> & intersect_coords, std::vector<double> & param_coords)
{
  // this is pretty cool
  // we will return only surfaces (gentities )
  //
  ErrorCode rval;

  unsigned int numfaces = _my_gsets[2].size();
  // do ray fire
  const double point[] = { x, y, z };
  const double dir[] = { dir_x, dir_y, dir_z };
  CartVect P(point);
  CartVect V(dir);

  //std::vector<double> distances;
  std::vector<EntityHandle> facets;
  //std::vector<EntityHandle> sets;
  unsigned int i;
  for (i = 0; i < numfaces; i++) {
    EntityHandle face = _my_gsets[2][i];
    EntityHandle rootForFace;
    rval = _my_geomTopoTool->get_root(face, rootForFace);
    MBERRORR(rval, "Failed to get root of face.");
    std::vector<double> distances_out;
    std::vector<EntityHandle> sets_out;
    std::vector<EntityHandle> facets_out;
    rval = _my_geomTopoTool->obb_tree()-> ray_intersect_sets(distances_out,
        sets_out, facets_out, rootForFace, tolerance,
        min_tolerace_intersections, point, dir);
    unsigned int j;
    for (j = 0; j < distances_out.size(); j++)
      param_coords.push_back(distances_out[j]);
    for (j = 0; j < sets_out.size(); j++)
      intersect_entity_handles.push_back(sets_out[j]);
    for (j = 0; j < facets_out.size(); j++)
      facets.push_back(facets_out[j]);

    MBERRORR(rval, "Failed to get ray intersections.");
  }
  // facets.size == distances.size()!!
  for (i = 0; i < param_coords.size(); i++) {
    CartVect intx = P + param_coords[i] * V;
    for (int j = 0; j < 3; j++)
      intersect_coords.push_back(intx[j]);

  }
  if (_smooth) {
    // correct the intersection point and the distance for smooth surfaces
    for (i = 0; i < intersect_entity_handles.size(); i++) {
      //EntityHandle geoSet = MBH_cast(sets[i]);
      SmoothFace* sFace = _faces[intersect_entity_handles[i]];
      // correct coordinates and distance from point
      /*moab::ErrorCode ray_intersection_correct(moab::EntityHandle facet, // (IN) the facet where the patch is defined
       moab::CartVect &pt, // (IN) shoot from
       moab::CartVect &ray, // (IN) ray direction
       moab::CartVect &eval_pt, // (INOUT) The intersection point
       double & distance, // (IN OUT) the new distance
       bool &outside);*/
      CartVect pos(&(intersect_coords[3 * i]));
      double dist = param_coords[i];
      bool outside = false;
      rval = sFace->ray_intersection_correct(facets[i], P, V, pos, dist,
          outside);
      MBERRORR(rval, "Failed to get better point on ray.");
      param_coords[i] = dist;

      for (int j = 0; j < 3; j++)
        intersect_coords[3 * i + j] = pos[j];
    }
  }
  return MB_SUCCESS;
}

ErrorCode FBEngine::getAdjacentEntities(const EntityHandle from,
    const int to_dim, Range &adjs)
{
  int this_dim = -1;
  for (int i = 0; i < 4; i++) {
    if (_my_gsets[i].find(from) != _my_gsets[i].end()) {
      this_dim = i;
      break;
    }
  }

  // check target dimension
  if (-1 == this_dim) {
    //ProcessError(iBase_FAILURE, "Entity not a geometry entity.");
    return MB_FAILURE;
  } else if (0 > to_dim || 3 < to_dim) {
    //ProcessError(iBase_FAILURE, "To dimension must be between 0 and 3.");
    return MB_FAILURE;
  } else if (to_dim == this_dim) {
    //ProcessError(iBase_FAILURE,
    //      "To dimension must be different from entity dimension.");
    return MB_FAILURE;
  }

  ErrorCode rval;
  adjs.clear();
  if (to_dim > this_dim) {
    int number;
    rval = MBI->num_parent_meshsets(from, &number, 0);
    rval = MBI->get_parent_meshsets(from, adjs);
    adjs.clear();
    rval = MBI->get_parent_meshsets(from, adjs, to_dim - this_dim);
  } else {
    int number;
    rval = MBI->num_child_meshsets(from, &number, 0);
    rval = MBI->get_child_meshsets(from, adjs);
    adjs.clear();
    rval = MBI->get_child_meshsets(from, adjs, this_dim - to_dim);
  }

  return MB_SUCCESS;
}

/*// new methods needed
 ErrorCode FBEngine::getTagHandle( const char* name, moab::Tag & handle_out )
 {
 return _mbImpl->tag_get_handle(name, handle_out);
 }*/

ErrorCode FBEngine::createTag(const char* /*tag_name*/, int /*tag_num_type_values*/,
                              int /*tag_type*/, Tag & /*tag_handle_out*/)
{
  // not implemented yet, some mapping needed for tag type
  return MB_FAILURE;
}

ErrorCode FBEngine::getArrData(const EntityHandle* entity_handles,
    int entity_handles_size, Tag tag_handle, void* tag_values_out)
{
  // responsibility of the user to have tag_values_out properly allocated
  // only some types of Tags are possible (double, int, etc)
  int tag_size;
  ErrorCode rval = MBI->tag_get_size(tag_handle, tag_size);
  if (MB_SUCCESS != rval)
    return rval;
  rval = MBI->tag_get_data(tag_handle, entity_handles, entity_handles_size,
      tag_values_out);
  return rval;
}

ErrorCode FBEngine::setArrData(const EntityHandle* entity_handles,
    int entity_handles_size, Tag tag_handle, const void* tag_values)
{
  // responsibility of the user to have tag_values_out properly allocated
  // only some types of Tags are possible (double, int, etc)
  int tag_size;
  ErrorCode rval = MBI->tag_get_size(tag_handle, tag_size);
  if (MB_SUCCESS != rval)
    return rval;
  rval = MBI->tag_set_data(tag_handle, entity_handles, entity_handles_size,
      tag_values);
  return rval;
}

ErrorCode FBEngine::getEntAdj(EntityHandle handle, int type_requested,
    Range & adjEnts)
{
  return getAdjacentEntities(handle, type_requested, adjEnts);
}

ErrorCode FBEngine::getEgFcSense(EntityHandle mbedge, EntityHandle mbface,
    int & sense_out)
{

  // this one is important, for establishing the orientation of the edges in faces
  // use senses
  std::vector<EntityHandle> faces;
  std::vector<int> senses; // 0 is forward and 1 is backward
  ErrorCode rval = _my_geomTopoTool->get_senses(mbedge, faces, senses);
  if (MB_SUCCESS != rval)
    return rval;

  for (unsigned int i = 0; i < faces.size(); i++) {
    if (faces[i] == mbface) {
      sense_out = senses[i];
      return MB_SUCCESS;
    }
  }
  return MB_FAILURE;

}
// we assume the measures array was allocated correctly
ErrorCode FBEngine::measure(const EntityHandle * moab_entities,
    int entities_size, double * measures)
{
  ErrorCode rval;
  for (int i = 0; i < entities_size; i++) {
    measures[i] = 0.;

    int type;
    EntityHandle gset = moab_entities[i];
    rval = getEntType(gset, &type);
    if (MB_SUCCESS != rval)
      return rval;
    if (type == 1) { // edge: get all edges part of the edge set
      Range entities;
      rval = MBI->get_entities_by_type(gset, MBEDGE, entities);
      if (MB_SUCCESS != rval)
        return rval;

      for (Range::iterator it = entities.begin(); it != entities.end(); it++) {
        EntityHandle edge = *it;
        CartVect vv[2];
        const EntityHandle *conn2 = NULL;
        int num_nodes;
        rval = MBI->get_connectivity(edge, conn2, num_nodes);
        if (MB_SUCCESS != rval || num_nodes != 2)
          return MB_FAILURE;
        rval = MBI->get_coords(conn2, 2, (double *) &(vv[0][0]));
        if (MB_SUCCESS != rval)
          return rval;

        vv[0] = vv[1] - vv[0];
        measures[i] += vv[0].length();
      }
    }
    if (type == 2) { // surface
      // get triangles in surface; TODO: quads!
      Range entities;
      rval = MBI->get_entities_by_type(gset, MBTRI, entities);
      if (MB_SUCCESS != rval)
        return rval;

      for (Range::iterator it = entities.begin(); it != entities.end(); it++) {
        EntityHandle tri = *it;
        CartVect vv[3];
        const EntityHandle *conn3 = NULL;
        int num_nodes;
        rval = MBI->get_connectivity(tri, conn3, num_nodes);
        if (MB_SUCCESS != rval || num_nodes != 3)
          return MB_FAILURE;
        rval = MBI->get_coords(conn3, 3, (double *) &(vv[0][0]));
        if (MB_SUCCESS != rval)
          return rval;

        vv[1] = vv[1] - vv[0];
        vv[2] = vv[2] - vv[0];
        vv[0] = vv[1] * vv[2];
        measures[i] += vv[0].length() / 2;// area of triangle
      }

    }
  }
  return MB_SUCCESS;
}

ErrorCode FBEngine::getEntNrmlSense(EntityHandle /*face*/, EntityHandle /*region*/,
                                    int& /*sense*/)
{
  return MB_NOT_IMPLEMENTED; // not implemented
}

ErrorCode FBEngine::getEgEvalXYZ(EntityHandle /*edge*/, double /*x*/, double /*y*/,
                                 double /*z*/, double& /*on_x*/, double& /*on_y*/, double& /*on_z*/, double& /*tngt_i*/,
                                 double& /*tngt_j*/, double& /*tngt_k*/, double& /*cvtr_i*/, double& /*cvtr_j*/,
                                 double& /*cvtr_k*/)
{
  return MB_NOT_IMPLEMENTED; // not implemented
}
ErrorCode FBEngine::getFcEvalXYZ(EntityHandle /*face*/, double /*x*/, double /*y*/,
                                 double /*z*/, double& /*on_x*/, double& /*on_y*/, double& /*on_z*/, double& /*nrml_i*/,
                                 double& /*nrml_j*/, double& /*nrml_k*/, double& /*cvtr1_i*/, double& /*cvtr1_j*/,
                                 double& /*cvtr1_k*/, double& /*cvtr2_i*/, double& /*cvtr2_j*/, double& /*cvtr2_k*/)
{
  return MB_NOT_IMPLEMENTED; // not implemented
}

ErrorCode FBEngine::getEgVtxSense(EntityHandle edge, EntityHandle vtx1,
    EntityHandle vtx2, int& sense)
{
  // need to decide first or second vertex
  // important for moab
  int type;

  EntityHandle v1, v2;
  ErrorCode rval = getEntType(vtx1, &type);
  if (MB_SUCCESS != rval || type != 0)
    return MB_FAILURE;
  // edge: get one vertex as part of the vertex set
  Range entities;
  rval = MBI->get_entities_by_type(vtx1, MBVERTEX, entities);
  if (MB_SUCCESS != rval)
    return rval;
  if (entities.size() < 1)
    return MB_FAILURE;
  v1 = entities[0]; // the first vertex
  entities.clear();
  rval = getEntType(vtx2, &type);
  if (MB_SUCCESS != rval || type != 0)
    return MB_FAILURE;
  rval = MBI->get_entities_by_type(vtx2, MBVERTEX, entities);
  if (MB_SUCCESS != rval)
    return rval;
  if (entities.size() < 1)
    return MB_FAILURE;
  v2 = entities[0]; // the first vertex
  entities.clear();
  // now get the edges, and get the first node and the last node in sequence of edges
  // the order is important...
  // these are ordered sets !!
  std::vector<EntityHandle> ents;
  rval = MBI->get_entities_by_type(edge, MBEDGE, ents);
  if (MB_SUCCESS != rval)
    return rval;
  if (ents.size() < 1)
    return MB_FAILURE;

  const EntityHandle* conn = NULL;
  int len;
  EntityHandle startNode, endNode;
  rval = MBI->get_connectivity(ents[0], conn, len);
  if (MB_SUCCESS != rval)
    return rval;
  startNode = conn[0];
  rval = MBI->get_connectivity(ents[ents.size() - 1], conn, len);
  if (MB_SUCCESS != rval)
    return rval;

  endNode = conn[1];
  sense = 1; //
  if ((startNode == endNode) && (v1 == startNode)) {
    sense = 0; // periodic
  }
  if ((startNode == v1) && (endNode == v2)) {
    sense = 1; // forward
  }
  if ((startNode == v2) && (endNode == v1)) {
    sense = -1; // reverse
  }
  return MB_SUCCESS;
}

ErrorCode FBEngine::getEntURange(EntityHandle edge, double& u_min,
    double& u_max)
{
  SmoothCurve * smoothCurve = _edges[edge];// this is a map
  // now, call smoothCurve methods
  smoothCurve -> get_param_range(u_min, u_max);
  return MB_SUCCESS;
}

ErrorCode FBEngine::getEntUtoXYZ(EntityHandle edge, double u, double& x,
    double& y, double& z)
{
  SmoothCurve * smoothCurve = _edges[edge];// this is a map
  // now, call smoothCurve methods
  smoothCurve -> position_from_u(u, x, y, z);
  return MB_SUCCESS;
}
ErrorCode FBEngine::isEntAdj(EntityHandle entity1, EntityHandle entity2,
    bool& adjacent_out)
{
  int type1, type2;
  ErrorCode rval = getEntType(entity1, &type1);
  if (MB_SUCCESS != rval)
    return rval;
  rval = getEntType(entity2, &type2);
  if (MB_SUCCESS != rval)
    return rval;

  Range adjs;
  if (type1 < type2) {
    rval = MBI->get_parent_meshsets(entity1, adjs, type2 - type1);
    if (MB_SUCCESS != rval)
      return rval;// MBERRORR("Failed to get parent meshsets in iGeom_isEntAdj.");
  } else {
    rval = MBI->get_child_meshsets(entity2, adjs, type2 - type1);
    if (MB_SUCCESS != rval)
      return rval;//MBERRORR("Failed to get child meshsets in iGeom_isEntAdj.");
  }

  adjacent_out = adjs.find(entity2) != _my_gsets[type2].end();

  return MB_SUCCESS;
}

ErrorCode FBEngine::split_surface_with_direction(EntityHandle face, std::vector<double> & xyz,
    double * direction, EntityHandle & newFace, int closed)
{

  // first of all, find all intersection points (piercing in the face along the direction)
  // assume it is robust; what if it is not sufficiently robust?

  ErrorCode rval;
  std::vector<CartVect> points; // if the point is interior, it will be matched to a triangle,
  // otherwise to edges or even nodes
  std::vector<EntityHandle> entities;
  //start with initial points, intersect along the direction, find the facets
  // then find the position
  int numIniPoints = (int) xyz.size() / 3;
  EntityHandle rootForFace;

  rval = _my_geomTopoTool->get_root(face, rootForFace);
  MBERRORR(rval, "Failed to get root of face.");

  const double dir[] = { direction[0], direction[1], direction[2] };
  std::vector<EntityHandle> nodes; // get the nodes closest to the ray traces of interest
  std::vector<EntityHandle> trianglesAlong;

  int i = 0;
  for (; i < numIniPoints; i++) {
    const double point[] = { xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2] };// or even point( &(xyz[3*i]) ); //
    std::vector<double> distances_out;
    std::vector<EntityHandle> sets_out;
    std::vector<EntityHandle> facets_out;
    rval = _my_geomTopoTool->obb_tree()-> ray_intersect_sets(distances_out,
        sets_out, facets_out, rootForFace, tolerance,
        min_tolerace_intersections, point, dir);
    MBERRORR(rval, "Failed to get ray intersections.");
    if (distances_out.size() < 1)
      MBERRORR(MB_FAILURE, "Failed to get one intersection point, bad direction.");

    if (distances_out.size() > 1) {
      std::cout
          << " too many intersection points. Only the first one considered\n";
    }
    std::vector<EntityHandle>::iterator pFace = std::find(sets_out.begin(), sets_out.end(), face);

    if (pFace == sets_out.end())
      MBERRORR(MB_FAILURE, "Failed to intersect given face, bad direction.");
    unsigned int index = pFace-sets_out.begin();
    // get the closest node of the triangle, and modify
    CartVect P(point);
    CartVect Dir(dir);
    CartVect newPoint = P + distances_out[index] * Dir;
    // get the triangle coordinates
    //
    int nnodes;
    const EntityHandle * conn3;
    rval = MBI->get_connectivity(facets_out[index], conn3, nnodes);
    MBERRORR(rval, "Failed to get connectivity");

    CartVect PP[3];
    rval = _mbImpl->get_coords(conn3, nnodes, (double*) &PP[0]);
    MBERRORR(rval, "Failed to get coordinates");

    EntityHandle vertex=conn3[0];
    double minD2=(newPoint-PP[0]).length_squared();
    for (int j=1; j<nnodes; j++) // nnodes should be 3, actually
    {
      double d2=(newPoint-PP[j]).length_squared();
      if ( d2 < minD2)
      {
        minD2 = d2;
        vertex = conn3[j];
      }
    }
    nodes.push_back(vertex);
  }
  // now, we have to find more intersection points, either interior to triangles, or on edges, or on vertices
  // use the same tolerance as before
  // starting from 2 points on 2 triangles, and having the direction, get more intersection points
  // between the plane formed by direction and those 2 points, and edges from triangulation (the triangles
  // involved will be part of the same gentity , original face ( moab set)
  //int closed = 1;// closed = 0 if the polyline is not closed

  CartVect Dir(direction);
  for (i = 0; i < numIniPoints - 1 + closed; i++) {
    int nextIndex = (i + 1) % numIniPoints;
    rval = compute_intersection_points(face, nodes[i], nodes[nextIndex], Dir, points,
        entities, trianglesAlong);
    MBERRORR(rval, "can't get intersection points");
  }
  // the segment between point_i and point_i+1 is in trianglesAlong_i
  // points_i is on entities_i
  rval = split_surface(face, points, entities, trianglesAlong, newFace);

  //
  return rval;
}
/**
 *  this method splits along the polyline defined by points and entities
 *  the polyline will be defined with
 *  // the entities are now only nodes and edges, no triangles!!!
 *  the first and last ones are also nodes for sure
 */
ErrorCode FBEngine::split_surface(EntityHandle face,
    std::vector<CartVect> & points, std::vector<EntityHandle> & entities,
    std::vector<EntityHandle> & triangles, EntityHandle & newFace)
{
  // if the last point is the same as the first point, assume a loop defined
  // in this case, crop
  // 2 successive entities are in the same triangle, split the triangle in 2, 3, or even
  // 4 triangles; do not care about quality first, for the quads that get generated.
  // in the first pass, split triangles , create new ones, and put the affected ones in
  // a group / set that will be deleted / removed; or is it better to modify an existing triangle?
  // plan: list of triangles, parallel to points, entities
  // modify each triangle, split it
  // use a fill to determine the new sets, up to the polyline
  // points on the polyline will be moved to the closest point location, with some constraints
  // then the sets will be reset, geometry recomputed. new vertices, new edges, etc.

  Range iniTris;
  ErrorCode rval;
  rval = _mbImpl -> get_entities_by_type(face, MBTRI, iniTris);
  MBERRORR(rval, "can't get initial triangles");
  
  int num_points = (int) points.size();
  // if first point is the same as last point, we have a loop for cropping
  // otherwise, we have a trimming line for splitting

  Range triToDelete;
  Range edgesToDelete;
  bool loop = false;
  if ( (points[0]-points[num_points-1]).length()<tolerance_segment )
    loop = true;

  // go along entities, find the triangles they belong to, and create new triangles if needed
  // (also new vertices; edges will be created later)

  // we need to start at an edge or at a vertex; we could start in an interior of a triangle, but
  // let's skip that
  std::vector<EntityHandle> nodesAlongPolyline;
  nodesAlongPolyline.push_back(entities[0]); // it is for sure a node
  for (int i = 0; i < num_points-1; i++) {
    EntityHandle tri = triangles[i]; // this is happening in triangle i
    EntityHandle e1 = entities[i];
    EntityHandle e2 = entities[i + 1];
    EntityType et1 = _mbImpl->type_from_handle(e1);
    //EntityHandle vertex1 = nodesAlongPolyline[i];// irrespective of the entity type i,
    // we already have the vertex there
    EntityType et2 = _mbImpl->type_from_handle(e2);
    if (et2 == MBVERTEX) {
      nodesAlongPolyline.push_back(e2);
    }
    else // if (et2==MBEDGE || et2==MBTRI)
    {
      CartVect coord_vert=points[i+1];
      EntityHandle newVertex;
      rval = _mbImpl->create_vertex((double*)&coord_vert, newVertex);
      MBERRORR(rval, "can't create vertex");
      nodesAlongPolyline.push_back(newVertex);
    }
    // if vertices, do not split anything, just get the edge for polyline
    if (et2 == MBVERTEX && et1 == MBVERTEX) {
      // nothing to do, just continue;
      continue; // continue the for loop
    }
    // create some new triangles, with common edge nodesAlongPolyline[i], nodesAlongPolyline[i+1]
    // first get the nodes of current triangle

    // at least one is edge; et1 cannot be tri, but et2 can
    // when et2 triangle, find next entity, (node or edge), and break all at the same time
    // there will be at most 5 triangles formed, and at least 4 in that case (also advance i)
    // many cases, if, then, else etc.
    // initially form only triangles, worry about edges later (or not?)
    /*if (MBTRI == et2)
    {

      // triangles [i+1] has to be the same, and point to a different entity (edge, or vertex)
      if (tri != triangles[i] || tri!=triangles[i+1])
      {
        // it should be the same triangle
        MBERRORR(MB_FAILURE, "redefine polyline, it is almost self-intersecting");
      }
      EntityHandle e3=entities[i+2];
      EntityHandle et3 = _mbImpl->type_from_handle(e3);
      // et3 should be opposed to et1
      //
      if (et3 == MBVERTEX) {
        nodesAlongPolyline.push_back(e3);
      }
      else // if (et2==MBEDGE || et2==MBTRI)
      {
        double coord_vert[3] = { points[3 *i+3], points[3 * i + 4], points[3 * i
            + 5] };
        EntityHandle newVertex;
        rval = _mbImpl->create_vertex(coord_vert, newVertex);
        MBERRORR(rval, "can't create vertex");
        nodesAlongPolyline.push_back(newVertex);
      }
      // so now we have the same triangle, make sure it is another edge or vertex
      rval = BreakTriangle( tri, e1, e3, nodesAlongPolyline[i], nodesAlongPolyline[i+1],
          nodesAlongPolyline[i+2], new_triangles);// nodesAlongPolyline are on entities!
      MBERRORR(rval, "can't break triangle");
      if (et3==MBEDGE)
          edgesToDelete.insert(e3);
    }
    else
    {*/
    if (debug_splits)
    {
      std::cout <<"tri: type: " << _mbImpl->type_from_handle(tri) << " id:" <<
          _mbImpl->id_from_handle(tri) << " e1:" << e1 << " e2:" << e2 << "\n";
    }
    rval = BreakTriangle2( tri, e1, e2, nodesAlongPolyline[i], nodesAlongPolyline[i+1]);
    MBERRORR(rval, "can't break triangle 2");
    if (et2==MBEDGE)
      edgesToDelete.insert(e2);
    triToDelete.insert(tri);

  }
  // nodesAlongPolyline will define a new geometric edge
  if (debug_splits)
  {
    std::cout<<"nodesAlongPolyline: " << nodesAlongPolyline.size() << "\n";
    std::cout << "points: " << num_points << "\n";
  }
  // if needed, create edges along polyline, or revert the existing ones, to
  // put them in a new edge set
  EntityHandle new_geo_edge;
  Range geo_vertices;
  rval = create_new_gedge(nodesAlongPolyline, new_geo_edge, geo_vertices);
  MBERRORR(rval, "can't create a new edge");

  // start from a triangle that is not in the triangles to delete
  // flood fill 

  if (!loop)
  {
    // we will have to split the boundary edges
    // first, find the actual boundary, and try to split with the 2 end points (nodes)
    // get the adjacent edges, and see which one has the end nodes

    rval = split_boundary(face, nodesAlongPolyline[0]);
    MBERRORR(rval, "can't split with first node");
    rval = split_boundary(face, nodesAlongPolyline[nodesAlongPolyline.size()-1]);
    MBERRORR(rval, "can't split with second node)");
  }
  // we will separate triangles to delete, unaffected, new_triangles,
  //  nodesAlongPolyline, 
  Range first, second;
  rval = separate (face, triToDelete, new_geo_edge,
      first, second);

  // now, we are done with the computations;
  // we need to put the new nodes on the smooth surface
  if (this->_smooth)
  {
    rval = smooth_new_intx_points(face, nodesAlongPolyline);
    MBERRORR(rval, "can't smooth new points");
  }

  // create the new set
  rval = _mbImpl->create_meshset(MESHSET_SET, newFace);
  MBERRORR(rval, "can't create a new face");

  _my_geomTopoTool->add_geo_set(newFace, 2);
  _my_geomTopoTool->add_geo_set(new_geo_edge, 1);
  for (unsigned int k=0; k< geo_vertices.size(); k++)
    _my_geomTopoTool->add_geo_set(geo_vertices[k], 0);

  // the new face will have the first set (positive sense triangles, to the left)
  rval = _mbImpl->add_entities(newFace, first);
  MBERRORR(rval, "can't add first range triangles to new face");
  // both faces will have the edge now
  rval = _mbImpl ->add_parent_child( face, new_geo_edge);
  MBERRORR(rval, "can't add parent child relations for new edge");

  rval = _mbImpl ->add_parent_child( newFace, new_geo_edge);
  MBERRORR(rval, "can't add parent child relations for new edge");
  // add senses
  // sense newFace is 1, old face is -1
  rval = _my_geomTopoTool-> set_sense( new_geo_edge, newFace, 1);
  MBERRORR(rval, "can't set sense for new edge");

  rval = _my_geomTopoTool-> set_sense( new_geo_edge, face, -1);
  MBERRORR(rval, "can't set sense for new edge in original face");

  rval = set_neumann_tags(face, newFace);
  MBERRORR(rval, "can't set NEUMANN set tags");

  // now, we should remove from the original set all tris, and put the "second" range
  rval = _mbImpl->remove_entities(face, iniTris);
  MBERRORR(rval, "can't remove original tris from initial face set");

  rval = _mbImpl->add_entities(face, second);
  MBERRORR(rval, "can't add second range to the original set");

  if (!loop)
  {
    rval = redistribute_boundary_edges_to_faces(face, newFace, new_geo_edge);
    MBERRORR(rval, "fail to reset the proper boundary faces");
  }

  if (_smooth)
    delete_smooth_tags();
  // this will remove the extra smooth faces and edges
  clean();
  // also, these nodes need to be moved to the smooth surface, sometimes before deleting the old
  // triangles
  // remove the triangles from the set, then delete triangles (also some edges need to be deleted!)
  rval=_mbImpl->delete_entities( triToDelete );
  MBERRORR(rval, "can't delete triangles");

  // delete edges that are broke up in 2
  rval=_mbImpl->delete_entities(edgesToDelete);
  MBERRORR(rval, "can't delete edges");

  if (debug_splits)
  {
    _mbImpl->write_file("newFace.vtk", "vtk", 0, &newFace, 1);
    _mbImpl->write_file("leftoverFace.vtk", "vtk", 0, &face, 1);
  }
  return MB_SUCCESS;
}

ErrorCode FBEngine::smooth_new_intx_points(EntityHandle face,
      std::vector<EntityHandle> & nodesAlongPolyline)
{
  // first of all, find out all nodes from face; do not move them if they are on the
  // original face
  std::vector<double> ini_coords;
  int num_points = (int)nodesAlongPolyline.size();
  ini_coords.resize(3*nodesAlongPolyline.size());
  ErrorCode rval = _mbImpl->get_coords(&(nodesAlongPolyline[0]),
      num_points, &(ini_coords[0]));
  MBERRORR(rval, "can't get coordinates");

  // do not move nodes from the original face
  // first get all triangles, and then all nodes from those triangles

  Range tris;
  rval = _mbImpl->get_entities_by_type(face, MBTRI, tris);
  MBERRORR(rval, "can't get triangles");

  Range ini_nodes;
  rval = _mbImpl->get_connectivity( tris, ini_nodes);
  MBERRORR(rval, "can't get connectivities");

  SmoothFace* smthFace = _faces[face];

  for (int i=0; i<num_points; i++)
  {
    if ( ini_nodes.find(nodesAlongPolyline[i])==ini_nodes.end() )
    {
      int i3=3*i;
      smthFace->move_to_surface(ini_coords[i3], ini_coords[i3+1], ini_coords[i3+2]);
      // reset the coordinates of this node
      rval = _mbImpl->set_coords(&(nodesAlongPolyline[i]), 1,
          &(ini_coords[i3]));
      MBERRORR(rval, "can't set new coordinates for a node");
    }
  }

  return MB_SUCCESS;
}
// we will use the fact that the splitting edge is oriented right now
// to the left will be new face, to the right, old face
// (to the left, positively oriented triangles)
ErrorCode FBEngine::separate (EntityHandle face, Range & triToDelete,
    EntityHandle new_geo_edge, Range & first,  Range & second)
{
  //Range unaffectedTriangles = subtract(iniTriangles, triToDelete);
  // insert in each
  // start with a new triangle, and flood to get the first range; what is left is the
  // second range
  // flood fill is considering edges adjacent to seed triangles; if there is
  //  an edge in the new_geo_edge, it is skipped; triangles in the
  // triangles to delete are not added
  // first, create all edges of the new triangles

  //
  // new face will have the new edge oriented positively
  // get mesh edges from geo edge (splitting gedge);

  Range mesh_edges;
  ErrorCode rval = _mbImpl->get_entities_by_type(new_geo_edge, MBEDGE, mesh_edges);
  MBERRORR(rval, "can't get new polyline edges");

  // get a positive triangle adjacent to mesh_edge[0]
  // add to first triangles to the left, second triangles to the right of the mesh_edges ;
  std::set<EntityHandle> firstSet;
  std::set<EntityHandle> secondSet;
  for (Range::iterator it = mesh_edges.begin(); it!=mesh_edges.end(); it++)
  {
    EntityHandle meshEdge = *it;
    Range adj_tri;
    rval =  _mbImpl->get_adjacencies(&meshEdge, 1,
            2, false, adj_tri);
    MBERRORR(rval, "can't get adj_tris to mesh edge");

    for ( Range::iterator it2=adj_tri.begin(); it2!=adj_tri.end(); it2++)
    {
      EntityHandle tr=*it2;
      int num1, sense, offset;
      rval = _mbImpl->side_number(tr, meshEdge, num1, sense, offset);
      MBERRORR(rval, "edge not adjacent");
      if (sense==1)
        firstSet.insert(tr);
      else
        secondSet.insert(tr);
    }
  }

  // get the first new triangle: will be part of first set;
  // flood fill first set, the rest will be in second set
  // the edges from new_geo_edge will not be crossed

  // get edges of face (adjacencies)
  // also get the old boundary edges, from face; they will be edges to not cross
  Range bound_edges;
  rval = getAdjacentEntities(face, 1, bound_edges);
  MBERRORR(rval, "can't get boundary edges");

  // add to the do not cross edges range, all edges from initial boundary
  Range initialBoundaryEdges;
  for (Range::iterator it= bound_edges.begin(); it!=bound_edges.end(); it++)
  {
    EntityHandle bound_edge=*it;
    rval = _mbImpl->get_entities_by_dimension(bound_edge, 1, initialBoundaryEdges);
  }

  Range doNotCrossEdges = unite(initialBoundaryEdges, mesh_edges);// add the splitting edges !

  std::queue<EntityHandle> firstQueue;
  for (std::set<EntityHandle>::iterator it3 = firstSet.begin(); it3!=firstSet.end(); it3++)
  {
    EntityHandle firstTri = *it3;
    firstQueue.push(firstTri);
  }
  std::set<EntityHandle> visited=firstSet;// already decided, do not care about them again
  while(!firstQueue.empty())
  {
    EntityHandle currentTriangle=firstQueue.front();
    firstQueue.pop();
    firstSet.insert(currentTriangle);
    // add new triangles that share an edge
    Range currentEdges;
    rval =  _mbImpl->get_adjacencies(&currentTriangle, 1,
        1, false, currentEdges, Interface::UNION);
    MBERRORR(rval, "can't get adjacencies");
    for (Range::iterator it=currentEdges.begin(); it!=currentEdges.end(); it++)
    {
      EntityHandle frontEdge= *it;
      if ( doNotCrossEdges.find(frontEdge)==doNotCrossEdges.end())
      {
        // this is an edge that can be crossed
        Range adj_tri;
        rval =  _mbImpl->get_adjacencies(&frontEdge, 1,
                2, false, adj_tri, Interface::UNION);
        MBERRORR(rval, "can't get adj_tris");
        // if the triangle is not in first range, add it to the queue
        for (Range::iterator it2=adj_tri.begin(); it2!=adj_tri.end(); it2++)
        {
          EntityHandle tri2=*it2;
          if ( (firstSet.find(tri2)==firstSet.end()) &&
                (triToDelete.find(tri2) == triToDelete.end())
              && (visited.find(tri2) == visited.end()) )
          {
            firstQueue.push(tri2);
          }
          visited.insert(tri2);
        }

      }
    }
  }
  // try a second queue
  std::queue<EntityHandle> secondQueue;
  for (std::set<EntityHandle>::iterator it4 = secondSet.begin(); it4!=secondSet.end(); it4++)
  {
    EntityHandle secondTri = *it4;
    secondQueue.push(secondTri);
  }
  visited=secondSet;// already decided, do not care about them again
  /*// now "first" should have one set of triangles
  // second = iniTriangles + new_triangles - triangles to delete - first
  second =  unite (iniTriangles, new_triangles);
  // now subtract the ones to delete and the first set
  Range second2 = subtract (second, triToDelete);
  second = subtract(second2, first);*/

  while(!secondQueue.empty())
  {
    EntityHandle currentTriangle=secondQueue.front();
    secondQueue.pop();
    secondSet.insert(currentTriangle);
    // add new triangles that share an edge
    Range currentEdges;
    rval =  _mbImpl->get_adjacencies(&currentTriangle, 1,
        1, false, currentEdges, Interface::UNION);
    MBERRORR(rval, "can't get adjacencies");
    for (Range::iterator it=currentEdges.begin(); it!=currentEdges.end(); it++)
    {
      EntityHandle frontEdge= *it;
      if ( doNotCrossEdges.find(frontEdge)==doNotCrossEdges.end())
      {
        // this is an edge that can be crossed
        Range adj_tri;
        rval =  _mbImpl->get_adjacencies(&frontEdge, 1,
                2, false, adj_tri, Interface::UNION);
        MBERRORR(rval, "can't get adj_tris");
        // if the triangle is not in first range, add it to the queue
        for (Range::iterator it2=adj_tri.begin(); it2!=adj_tri.end(); it2++)
        {
          EntityHandle tri2=*it2;
          if ( (secondSet.find(tri2)==secondSet.end()) &&
                (triToDelete.find(tri2) == triToDelete.end())
              && (visited.find(tri2) == visited.end()) )
          {
            secondQueue.push(tri2);
          }
          visited.insert(tri2);
        }

      }
    }
  }


  if (debug_splits)
  {
    std::cout << "first size: " << first.size() << "  second size:" << second.size() << "\n";

  }

  // now create first and second ranges, from firstSet and secondSet
  std::copy(firstSet.rbegin(), firstSet.rend(), range_inserter(first));
  std::copy(secondSet.rbegin(), secondSet.rend(), range_inserter(second));
  return MB_SUCCESS;
}
// if there is an edge between 2 nodes, then check it's orientation, and revert it if needed
ErrorCode  FBEngine::create_new_gedge(std::vector<EntityHandle> &nodesAlongPolyline, EntityHandle & new_geo_edge,
    Range & geo_vertices)
{

  ErrorCode rval = _mbImpl->create_meshset(MESHSET_ORDERED, new_geo_edge);
  MBERRORR(rval, "can't create geo edge");

  // now, get the edges, or create if not existing
  std::vector<EntityHandle> mesh_edges;
  for (unsigned int i=0; i<nodesAlongPolyline.size()-1; i++)
  {
    EntityHandle n1 = nodesAlongPolyline[i], n2 = nodesAlongPolyline[i+1];

    EntityHandle nn2[2];
    nn2[0] = n1;
    nn2[1] = n2;

    std::vector<EntityHandle> adjacent;
    rval = _mbImpl->get_adjacencies(nn2, 2, 1, false, adjacent,
                Interface::INTERSECT);
    // see if there is an edge between those 2 already
    if (adjacent.size()>=1)
    {
      // check the orientation, and reverse if needed
      const EntityHandle * conn2;
      int nnodes;
      rval = _mbImpl->get_connectivity(adjacent[0], conn2, nnodes);
      MBERRORR(rval, "can't get connectivity");
      if (conn2[0]==nn2[0] && conn2[1]==nn2[1])
      {
        // everything is fine
        mesh_edges.push_back(adjacent[0]);
      }
      else
      {
        // reset connectivity for the edge!
        rval = _mbImpl->set_connectivity(adjacent[0], nn2, 2);
        MBERRORR(rval, "can't reset connectivity");
        mesh_edges.push_back(adjacent[0]);
      }
    }
    else
    {
      // there is no edge between n1 and n2, create one
      EntityHandle mesh_edge;
      rval = _mbImpl->create_element(MBEDGE, nn2, 2, mesh_edge);
      MBERRORR(rval, "Failed to create a new edge");
      mesh_edges.push_back(mesh_edge);
    }
  }

  // add loops edges to the edge set
  rval = _mbImpl->add_entities(new_geo_edge, &mesh_edges[0], mesh_edges.size());// only one edge
  MBERRORR(rval, "can't add edges to new_geo_set");
  // check vertex sets for vertex 1 and vertex 2?
  // get all sets of dimension 0 from database, and see if our ends are here or not

  Range ends_geo_edge;
  ends_geo_edge.insert(nodesAlongPolyline[0]);
  ends_geo_edge.insert(nodesAlongPolyline[nodesAlongPolyline.size()-1]);

  for (unsigned int k = 0; k<ends_geo_edge.size(); k++ )
  {
    EntityHandle node = ends_geo_edge[k];
    EntityHandle nodeSet;
    bool found=find_vertex_set_for_node(node, nodeSet);

    if (!found)
    {
      // create a node set and add the node

      rval = _mbImpl->create_meshset(MESHSET_SET, nodeSet);
      MBERRORR(rval, "Failed to create a new vertex set");

      rval = _mbImpl->add_entities(nodeSet, &node, 1);
      MBERRORR(rval, "Failed to add the node to the set");

      rval = _my_geomTopoTool->add_geo_set(nodeSet, 0);//
      MBERRORR(rval, "Failed to commit the node set");

      if (debug_splits)
      {
        std::cout<<" create a vertex set " << _mbImpl->id_from_handle(nodeSet) << " global id:"<<
            this->_my_geomTopoTool->global_id(nodeSet) << " for node " << node <<  "\n";
      }

    }
    geo_vertices.insert(nodeSet); //(new or already existing)
    // arrange
      //
    rval = _mbImpl ->add_parent_child( new_geo_edge, nodeSet);
    MBERRORR(rval, "Failed to add parent child relation");
  }

  return rval;
}

void FBEngine::print_debug_triangle(EntityHandle t)
{
  std::cout<< " triangle id:" << _mbImpl->id_from_handle(t) << "\n";
  const EntityHandle * conn3;
  int nnodes;
  _mbImpl->get_connectivity(t, conn3, nnodes);
  // get coords
  CartVect P[3];
  _mbImpl->get_coords(conn3, 3, (double*) &P[0]);
  std::cout <<"  nodes:" << conn3[0] << " " << conn3[1] << " " << conn3[2] << "\n";
  CartVect PP[3];
  PP[0] = P[1]-P[0];
  PP[1] = P[2]-P[1];
  PP[2] = P[0] - P[2];

  std::cout <<"  pos:" <<  P[0] << " " << P[1] << " " << P[2] << "\n";
  std::cout <<"   x,y diffs " <<  PP[0][0] <<" " << PP[0][1]  << ",  " << PP[1][0] <<" " << PP[1][1]
                   << ",  " << PP[2][0] <<" " << PP[2][1]  << "\n";
  return;
}
// actual breaking of triangles
// case 1: n2 interior to triangle
ErrorCode FBEngine::BreakTriangle(EntityHandle tri, EntityHandle e1, EntityHandle e3,
    EntityHandle n1, EntityHandle n2, EntityHandle n3)
{
  std::cout<< "FBEngine::BreakTriangle not implemented yet\n";
  return MB_FAILURE;
}
// case 2, n1 and n2 on boundary
ErrorCode FBEngine::BreakTriangle2(EntityHandle tri, EntityHandle e1, EntityHandle e2, EntityHandle n1,
  EntityHandle n2)// nodesAlongPolyline are on entities!
{
  // we have the nodes, we just need to reconnect to form new triangles
  ErrorCode rval;
  const EntityHandle * conn3;
  int nnodes;
  rval = _mbImpl->get_connectivity(tri, conn3, nnodes);
  MBERRORR(rval, "Failed to get connectivity");
  assert(3 == nnodes);

  EntityType et1= _mbImpl->type_from_handle(e1);
  EntityType et2= _mbImpl->type_from_handle(e2);

  if (MBVERTEX == et1)
  {
    // find the vertex in conn3, and form 2 other triangles
    int index =-1;
    for (index=0; index<3; index++)
    {
      if (conn3[index]==e1)// also n1
        break;
    }
    if (index==3)
      return MB_FAILURE;
    // 2 triangles: n1, index+1, n2, and n1, n2, index+2
    EntityHandle conn[6]={ n1, conn3[(index+1)%3], n2, n1, n2, conn3[(index+2)%3]};
    EntityHandle newTriangle;
    rval = _mbImpl->create_element(MBTRI, conn, 3, newTriangle);
    MBERRORR(rval, "Failed to create a new triangle");
    if (debug_splits)
      print_debug_triangle(newTriangle);
    rval = _mbImpl->create_element(MBTRI, conn+3, 3, newTriangle);// the second triangle
    MBERRORR(rval, "Failed to create a new triangle");
    if (debug_splits)
      print_debug_triangle(newTriangle);
    return MB_SUCCESS;
  }
  else if (MBVERTEX == et2)
  {
    int index =-1;
    for (index=0; index<3; index++)
    {
      if (conn3[index]==e2) // also n2
        break;
    }
    if (index==3)
      return MB_FAILURE;
    // 2 triangles: n1, index+1, n2, and n1, n2, index+2
    EntityHandle conn[6]={ n2, conn3[(index+1)%3], n1, n2, n1, conn3[(index+2)%3]};
    EntityHandle newTriangle;
    rval = _mbImpl->create_element(MBTRI, conn, 3, newTriangle);
    MBERRORR(rval, "Failed to create a new triangle");
    if (debug_splits)
          print_debug_triangle(newTriangle);
    rval = _mbImpl->create_element(MBTRI, conn+3, 3, newTriangle);// the second triangle
    MBERRORR(rval, "Failed to create a new triangle");
    if (debug_splits)
          print_debug_triangle(newTriangle);
    return MB_SUCCESS;
  }
  else
  {
    // both are edges adjacent to triangle tri
    // there are several configurations possible for n1, n2, between conn3 nodes.
    int num1, num2, sense, offset;
    rval = _mbImpl->side_number(tri, e1, num1, sense, offset);
    MBERRORR(rval, "edge not adjacent");

    rval = _mbImpl->side_number(tri, e2, num2, sense, offset);
    MBERRORR(rval, "edge not adjacent");

    const EntityHandle * conn12; // connectivity for edge 1
    const EntityHandle * conn22; // connectivity for edge 2
    //int nnodes;
    rval = _mbImpl->get_connectivity(e1, conn12, nnodes);
    MBERRORR(rval, "Failed to get connectivity of edge 1");
    assert(2 == nnodes);
    rval = _mbImpl->get_connectivity(e2, conn22, nnodes);
    MBERRORR(rval, "Failed to get connectivity of edge 2");
    assert(2 == nnodes);
    // now, having the side number, work through
    if (debug_splits)
    {
      std::cout << "tri conn3:" << conn3[0] << " "<< conn3[1] <<" " << conn3[2] << "\n";
      std::cout << " edge1: conn12:" << conn12[0] << " "<< conn12[1] <<"  side: " << num1 << "\n";
      std::cout << " edge2: conn22:" << conn22[0] << " "<< conn22[1] <<"  side: " << num2 << "\n";
    }
    int unaffectedSide = 3-num1-num2;
    int i3 = (unaffectedSide+2)%3;// to 0 is 2, to 1 is 0, to 2 is 1
    // triangles will be formed with triVertexIndex , n1, n2 (in what order?)
    EntityHandle v1, v2; // to hold the 2 nodes on edges
    if (num1==i3)
    {
      v1 = n1;
      v2 = n2;
    }
    else // if (num2==i3)
    {
      v1 = n2;
      v2 = n1;
    }
    // three triangles are formed
    int i1 = (i3+1)%3;
    int i2 = (i3+2)%3;
    // we could break the surface differently
    EntityHandle conn[9]={ conn3[i3], v1, v2, v1, conn3[i1], conn3[i2],
        v2, v1,  conn3[i2]};
    EntityHandle newTriangle;
    if (debug_splits)
       std::cout << "Split 2 edges :\n";
    rval = _mbImpl->create_element(MBTRI, conn, 3, newTriangle);
    MBERRORR(rval, "Failed to create a new triangle");
    if (debug_splits)
          print_debug_triangle(newTriangle);
    rval = _mbImpl->create_element(MBTRI, conn+3, 3, newTriangle);// the second triangle
    MBERRORR(rval, "Failed to create a new triangle");
    if (debug_splits)
          print_debug_triangle(newTriangle);
    rval = _mbImpl->create_element(MBTRI, conn+6, 3, newTriangle);// the second triangle
    MBERRORR(rval, "Failed to create a new triangle");
    if (debug_splits)
          print_debug_triangle(newTriangle);
    return MB_SUCCESS;
  }

  return MB_SUCCESS;
}

// build the list of intx points and entities from database involved
// vertices, edges, triangles
//it could be just a list of vertices (easiest case to handle after)

ErrorCode FBEngine::compute_intersection_points(EntityHandle & face,
    EntityHandle from, EntityHandle to,
    CartVect & Dir, std::vector<CartVect> & points,
    std::vector<EntityHandle> & entities, std::vector<EntityHandle> & triangles)
{
  // keep a stack of triangles to process, and do not add those already processed
  // either mark them, or maybe keep them in a local set?
  // from and to are now nodes, start from them
  CartVect p1, p2;// the position of from and to
  ErrorCode rval = _mbImpl->get_coords(&from, 1, (double *)&p1);
  MBERRORR(rval, "failed to get 'from' coordinates");
  rval = _mbImpl->get_coords(&to, 1, (double *)&p2);
  MBERRORR(rval, "failed to get 'from' coordinates");

  CartVect vect(p2 - p1);
  double dist2 = vect.length();
  if (dist2 < tolerance_segment) {
    // we are done, return
    return MB_SUCCESS;
  }
  CartVect normPlane = Dir * vect;
  normPlane.normalize();
  std::set<EntityHandle> visitedTriangles;
  CartVect currentPoint = p1;
  // push the current point if it is empty
  if (points.size() == 0) {
    points.push_back(p1);
    entities.push_back(from);// this is a node now
  }

  // these will be used in many loops
  CartVect intx = p1;// somewhere to start
  double param = -1.;

  // first intersection
  EntityHandle currentBoundary = from;// it is a node, in the beginning

  vect = p2 - currentPoint;
  while (vect.length() > 0.) {
    // advance towards "to" node, from boundary handle
    EntityType etype = _mbImpl->type_from_handle(currentBoundary);
    //if vertex, look for other triangles connected which intersect our plane (defined by p1, p2, dir)
    std::vector<EntityHandle> adj_tri;
    rval = _mbImpl->get_adjacencies(&currentBoundary, 1, 2, false, adj_tri);
    unsigned int j = 0;
    EntityHandle tri;
    for (; j < adj_tri.size(); j++) {
      tri = adj_tri[j];
      if (visitedTriangles.find(tri) != visitedTriangles.end())
        continue;// get another triangle, this one was already visited
      // if vertex, look for opposite edge
      // if edge, look for 2 opposite edges
      // get vertices
      int nnodes;
      const EntityHandle * conn3;
      rval = _mbImpl->get_connectivity(tri, conn3, nnodes);
      MBERRORR(rval, "Failed to get connectivity");
      // if one of the nodes is to, stop right there
      {
        if (conn3[0]==to || conn3[1]==to || conn3[2]==to)
        {
          visitedTriangles.insert(tri);
          triangles.push_back(tri);
          currentPoint = p2;
          points.push_back(p2);
          entities.push_back(to);// we are done
          break;// this is break from for loop, we still need to get out of while
          // we will get out, because vect will become 0, (p2-p2)
        }
      }
      EntityHandle nn2[2];
      if (MBVERTEX == etype) {
        nn2[0] = conn3[0];
        nn2[1] = conn3[1];
        if (nn2[0] == currentBoundary)
          nn2[0] = conn3[2];
        if (nn2[1] == currentBoundary)
          nn2[1] = conn3[2];
        // get coordinates
        CartVect Pt[2];

        rval = _mbImpl->get_coords(nn2, 2, (double*) &Pt[0]);
        MBERRORR(rval, "Failed to get coordinates");
        // see the intersection
        if (intersect_segment_and_plane_slice(Pt[0], Pt[1], currentPoint, p2,
            Dir, normPlane, intx, param)) {
          // we should stop for loop, and decide if it is edge or vertex
          if (param == 0.)
            currentBoundary = nn2[0];
          else {
            if (param == 1.)
              currentBoundary = nn2[1];
            else // param between 0 and 1, so edge
            {
              //find the edge between vertices
              std::vector<EntityHandle> edges1;
              rval = _mbImpl->get_adjacencies(nn2, 2, 1, false, edges1,
                  Interface::INTERSECT);
              MBERRORR(rval, "Failed to get edges");
              if (edges1.size() != 1)
                MBERRORR(MB_FAILURE, "Failed to get adjacent edges to 2 nodes");
              currentBoundary = edges1[0];
            }
          }
          visitedTriangles.insert(tri);
          currentPoint = intx;
          points.push_back(intx);
          entities.push_back(currentBoundary);
          triangles.push_back(tri);
          if (debug_splits)
            std::cout << "vtx new tri : " << _mbImpl->id_from_handle(tri)
                << " type bdy:" << _mbImpl->type_from_handle(currentBoundary)
                << "\n";
          break; // out of for loop over triangles

        }
      } else // this is MBEDGE, we have the other triangle to try out
      {
        //first find the nodes from existing boundary
        int nnodes2;
        const EntityHandle * conn2;
        rval = _mbImpl->get_connectivity(currentBoundary, conn2, nnodes2);
        MBERRORR(rval, "Failed to get connectivity");
        EntityHandle thirdNode = conn3[0];
        int thirdIndex = -1;
        for (int j = 0; j < 3; j++) {
          if ((conn3[j] != conn2[0]) && (conn3[j] != conn2[1])) {
            thirdNode = conn3[j];
            thirdIndex = j;
            break;
          }
        }
        if (-1 == thirdIndex)
          MBERRORR(MB_FAILURE, " can't get third node");
        CartVect Pt[3];
        rval = _mbImpl->get_coords(conn3, 3, (double*) &Pt[0]);
        int indexFirst = (thirdIndex + 1) % 3;
        int indexSecond = (thirdIndex + 2) % 3;
        int index[2] = { indexFirst, indexSecond };
        for (int k = 0; k < 2; k++) {
          EntityHandle nn2[2] = { conn3[index[k]], conn3[thirdIndex] };
          if (intersect_segment_and_plane_slice(Pt[index[k]], Pt[thirdIndex],
              currentPoint, p2, Dir, normPlane, intx, param)) {
            // we should stop for loop, and decide if it is edge or vertex
            if (param == 0.)
              currentBoundary = conn3[index[k]];//it is not really possible
            else {
              if (param == 1.)
                currentBoundary = conn3[thirdIndex];
              else // param between 0 and 1, so edge is fine
              {
                //find the edge between vertices
                std::vector<EntityHandle> edges1;
                rval = _mbImpl->get_adjacencies(nn2, 2, 1, false, edges1,
                    Interface::INTERSECT);
                MBERRORR(rval, "Failed to get edges");
                if (edges1.size() != 1)
                  MBERRORR(MB_FAILURE, "Failed to get adjacent edges to 2 nodes");
                currentBoundary = edges1[0];
              }
            }
            visitedTriangles.insert(tri);
            currentPoint = intx;
            points.push_back(intx);
            entities.push_back(currentBoundary);
            triangles.push_back(tri);
            if (debug_splits)
              std::cout << "edge new tri : " << _mbImpl->id_from_handle(tri)
                  << "  type bdy: " << _mbImpl->type_from_handle(
                  currentBoundary) << "\n";
            break; // out of for loop over triangles

          }
          // we should not reach here
        }

      }

    }
    /*if (j==adj_tri.size())
     {
     MBERRORR(MB_FAILURE, "did not advance");
     }*/
    vect = p2 - currentPoint;

  }


  if (debug_splits)
    std::cout << "nb entities: " << entities.size() <<  " triangles:" << triangles.size() <<
     " points.size()/3: " << points.size()/3 <<  "\n";

  return MB_SUCCESS;
}

ErrorCode  FBEngine::split_edge_at_point(EntityHandle edge, CartVect & point,
    EntityHandle & new_edge)
{
  //return MB_NOT_IMPLEMENTED;
  // first, we need to find the closest point on the smooth edge, then
  // split the smooth edge, then call the split_edge_at_mesh_node
  // or maybe just find the closest node??
  // first of all, we need to find a point on the smooth edge, close by
  // then split the mesh edge (if needed)
  // this would be quite a drama, as the new edge has to be inserted in
  // order for proper geo edge definition

  // first of all, find a closest point
  // usually called for
  if (debug_splits)
  {
    std::cout<<"Split edge " << _mbImpl->id_from_handle(edge) << " at point:" <<
        point << "\n";
  }
  int dim = _my_geomTopoTool->dimension(edge);
  if (dim !=1)
    return MB_FAILURE;
  if (!_smooth)
    return MB_FAILURE; // call it only for smooth option...
  // maybe we should do something for linear too

  SmoothCurve * curve = this->_edges[edge];
  EntityHandle closeNode;
  int edgeIndex;
  double u = curve-> u_from_position(point[0], point[1], point[2],
      closeNode, edgeIndex) ;
  if (0==closeNode)
  {
    // we really need to split an existing edge
    // do not do that yet
    // evaluate from u:
    /*double pos[3];
    curve->position_from_u(u,  pos[0], pos[1], pos[2] );*/
    // create a new node here, and split one edge
    // change also connectivity, create new triangles on both sides ...
    std::cout << "not found a close node,  u is: " << u << " edge index: " <<
        edgeIndex << "\n";
    return MB_FAILURE;// not implemented yet
  }
  return split_edge_at_mesh_node(edge, closeNode, new_edge);

}

ErrorCode FBEngine::split_edge_at_mesh_node(EntityHandle edge, EntityHandle node,
    EntityHandle & new_edge)
{
  // the node should be in the list of nodes

  int dim = _my_geomTopoTool->dimension(edge);
  if (dim!=1)
    return MB_FAILURE; // not an edge

  if (debug_splits)
  {
    std::cout<<"Split edge " << _mbImpl->id_from_handle(edge) << " with global id: "<<
        _my_geomTopoTool->global_id(edge) << " at node:" <<
        _mbImpl->id_from_handle(node) << "\n";
  }

  // now get the edges
  // the order is important...
  // these are ordered sets !!
  std::vector<EntityHandle> ents;
  ErrorCode rval = _mbImpl->get_entities_by_type(edge, MBEDGE, ents);
  if (MB_SUCCESS != rval)
    return rval;
  if (ents.size() < 1)
    return MB_FAILURE; // no edges

  const EntityHandle* conn = NULL;
  int len;
  // find the edge connected to the splitting node
  int num_mesh_edges = (int)ents.size();
  int index_edge;
  EntityHandle firstNode;
  for (index_edge = 0; index_edge<num_mesh_edges; index_edge++)
  {
    rval = MBI->get_connectivity(ents[index_edge], conn, len);
    if (MB_SUCCESS != rval)
      return rval;
    if (index_edge == 0)
      firstNode = conn[0];// will be used to decide vertex sets adjacencies
    if (conn[0] == node)
    {
      if (index_edge==0)
      {
        new_edge = 0; // no need to split, it is the first node
        return MB_SUCCESS; // no split
      }
      else
        return MB_FAILURE; // we should have found the node already , wrong
                           // connectivity
    }
    else if (conn[1] == node)
    {
      // we found the index of interest
      break;
    }
  }
  if (index_edge==num_mesh_edges-1)
  {
    new_edge = 0; // no need to split, it is the last node
    return MB_SUCCESS; // no split
  }

  // here, we have 0 ... index_edge edges in the first set,
  // create a vertex set and add the node to it

  if (debug_splits)
  {
    std::cout<<"Split edge with " << num_mesh_edges << " mesh edges, at index (0 based) " <<
        index_edge << "\n";
  }

  // at this moment, the node set should have been already created in new_geo_edge
  EntityHandle nodeSet; // the node set that has the node (find it!)
  bool found=find_vertex_set_for_node(node, nodeSet);

  if (!found) {
    // create a node set and add the node

    // must be an error, but create one nevertheless
    rval = _mbImpl->create_meshset(MESHSET_SET, nodeSet);
    MBERRORR(rval, "Failed to create a new vertex set");

    rval = _mbImpl->add_entities(nodeSet, &node, 1);
    MBERRORR(rval, "Failed to add the node to the set");

    rval = _my_geomTopoTool->add_geo_set(nodeSet, 0);//
    MBERRORR(rval, "Failed to commit the node set");

    if (debug_splits)
    {
      std::cout<<" create a vertex set (this should have been created before!)" << _mbImpl->id_from_handle(nodeSet) << " global id:"<<
          this->_my_geomTopoTool->global_id(nodeSet) <<  "\n";
    }
  }

  // we need to remove the remaining mesh edges from first set, and add it
  // to the second set, in order

  rval = _mbImpl->create_meshset(MESHSET_ORDERED, new_edge);
  MBERRORR(rval, "can't create geo edge");

  int remaining= num_mesh_edges - 1 - index_edge;

  // add edges to the edge set
  rval = _mbImpl->add_entities(new_edge, &ents[index_edge+1], remaining);
  MBERRORR(rval, "can't add edges to the new edge");

  // also, remove the second node set from old edge
  // remove the edges index_edge+1 and up

  rval = _mbImpl->remove_entities(edge, &ents[index_edge+1], remaining);
  MBERRORR(rval, "can't remove edges from the old edge");

  // need to find the adjacent vertex sets
  Range vertexRange;
  rval = getAdjacentEntities(edge, 0, vertexRange);

  EntityHandle firstSet, secondSet;
  if (vertexRange.size() == 1)
  {
    // initially a periodic edge, OK to add the new set to both edges, and the
    // second set
    firstSet = secondSet = vertexRange[0];
  }
  else
  {
    if (vertexRange.size() > 2)
      return MB_FAILURE; // something must be wrong with too many vertices
    // find first node
    int k;
    for (k=0; k<2; k++)
    {
      Range verts;
      rval = _mbImpl->get_entities_by_type(vertexRange[k], MBVERTEX, verts);

      MBERRORR(rval, "can't get vertices from vertex set");
      if (verts.size()!=1)
         MBERRORR(MB_FAILURE, " node set not defined well");
      if (firstNode == verts[0])
      {
        firstSet = vertexRange[k];
        secondSet = vertexRange[1-k]; // the other set; it is 1 or 0
        break;
      }
    }
    if (k>=2)
    {
      // it means we didn't find the right node set
      MBERRORR(MB_FAILURE, " can't find the right vertex");
    }
    // remove the second set from the connectivity with the
    //  edge (parent-child relation)
    //remove_parent_child( EntityHandle parent,
     //                                          EntityHandle child )
    rval = _mbImpl->remove_parent_child(edge, secondSet);
    MBERRORR(rval, " can't remove the second vertex from edge");
  }
  // at this point, we just need to add to both edges the new node set (vertex)
  rval = _mbImpl->add_parent_child(edge, nodeSet);
  MBERRORR(rval, " can't add new vertex to old edge");

  rval = _mbImpl->add_parent_child(new_edge, nodeSet);
  MBERRORR(rval, " can't add new vertex to new edge");

  // one thing that I forgot: add the secondSet as a child to new edge!!!
  // (so far, the new edge has only one end vertex!)
  rval = _mbImpl->add_parent_child(new_edge, secondSet);
  MBERRORR(rval, " can't add second vertex to new edge");

// now, add the edge and vertex to geo tool

  rval = _my_geomTopoTool->add_geo_set(new_edge, 1);
  MBERRORR(rval, " can't add new edge");

  // next, get the adjacent faces to initial edge, and add them as parents
  // to the new edge

  // need to find the adjacent face sets
  Range faceRange;
  rval = getAdjacentEntities(edge, 2, faceRange);

  // these faces will be adjacent to the new edge too!
  // need to set the senses of new edge within faces

  for (Range::iterator it= faceRange.begin(); it!=faceRange.end(); it++)
  {
    EntityHandle face = *it;
    rval = _mbImpl->add_parent_child(face, new_edge);
    MBERRORR(rval, " can't add new edge - face parent relation");
    int sense;
    rval = _my_geomTopoTool->get_sense(edge, face, sense);
    MBERRORR(rval, " can't get initial sense of edge in the adjacent face");
    // keep the same sense for the new edge within the faces
    rval = _my_geomTopoTool->set_sense(new_edge, face, sense);
    MBERRORR(rval, " can't set sense of new edge in the adjacent face");
  }

  return MB_SUCCESS;
}

ErrorCode FBEngine::split_boundary(EntityHandle face, EntityHandle atNode)
{
  // find the boundary edges, and split the one that we find it is a part of
  if (debug_splits)
  {
    std::cout<<"Split face " << _mbImpl->id_from_handle(face) << " at node:" <<
        _mbImpl->id_from_handle(atNode) << "\n";
  }
  Range bound_edges;
  ErrorCode rval = getAdjacentEntities(face, 1, bound_edges);
  MBERRORR(rval, " can't get boundary edges");
  for (Range::iterator it =bound_edges.begin(); it!=bound_edges.end(); it++ )
  {
    EntityHandle b_edge = *it;
    // get all edges in range
    Range mesh_edges;
    rval = _mbImpl->get_entities_by_dimension(b_edge, 1,
        mesh_edges);
    MBERRORR(rval, " can't get mesh edges");
    Range nodes;
    rval = _mbImpl->get_connectivity(mesh_edges, nodes);
    MBERRORR(rval, " can't get nodes from mesh edges");

    if (nodes.find(atNode)!=nodes.end())
    {
      // we found our boundary edge candidate
      EntityHandle new_edge;
      rval = split_edge_at_mesh_node(b_edge, atNode, new_edge);
      return rval;
    }
  }
  MBERRORR(MB_FAILURE, " we did not find an appropriate boundary edge"); ; //
}

bool FBEngine::find_vertex_set_for_node(EntityHandle iNode, EntityHandle & oVertexSet)
{
  bool found = false;
  Range vertex_sets;

  const int zero = 0;
  const void* const zero_val[] = { &zero };
  Tag geom_tag;
  ErrorCode rval = MBI->tag_get_handle(GEOM_DIMENSION_TAG_NAME, geom_tag);
  if (MB_SUCCESS!=rval)
    return false;
  rval = _mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &geom_tag,
        zero_val, 1, vertex_sets);
  if (MB_SUCCESS!=rval)
    return false;
  // local _gsets, as we might have not updated the local lists
  // see if ends of geo edge generated is in a geo set 0 or not

  for( Range::iterator vsit=vertex_sets.begin(); vsit!=vertex_sets.end(); vsit++)
  {
    EntityHandle vset=*vsit;
    // is the node part of this set?
    if (_mbImpl->contains_entities(vset, &iNode, 1))
    {

      found = true;
      oVertexSet = vset;
      break;
    }
  }
  return found;
}
ErrorCode FBEngine::redistribute_boundary_edges_to_faces(EntityHandle face, EntityHandle newFace,
      EntityHandle new_geo_edge)
{

  // so far, original boundary edges are all parent/child relations for face
  // we should get them all, and see if they are truly adjacent to face or newFace
  // decide also on the orientation/sense with respect to the triangles
  Range r1; // range in old face
  Range r2; // range of tris in new face
  ErrorCode rval = _mbImpl->get_entities_by_dimension(face, 2, r1);
  MBERRORR(rval, " can't get triangles from old face");
  rval = _mbImpl->get_entities_by_dimension(newFace, 2, r2);
  MBERRORR(rval, " can't get triangles from new face");
  // right now, all boundary edges are children of face
  // we need to get them all, and verify if they indeed are adjacent to face
  Range children;
  rval = _mbImpl->get_child_meshsets(face, children);// only direct children are of interest
  MBERRORR(rval, " can't get children sets from face");

  for (Range::iterator it = children.begin(); it!=children.end(); it++)
  {
    EntityHandle edge=*it;
    if (new_geo_edge==edge)
      continue; // we already set this one fine
    // get one mesh edge from the edge; we have to get all of them!!
    if (_my_geomTopoTool->dimension(edge)!=1)
      continue; // not an edge
    Range mesh_edges;
    rval = _mbImpl->get_entities_by_handle(edge, mesh_edges);
    MBERRORR(rval, " can't get mesh edges from edge");
    if (mesh_edges.empty())
      MBERRORR(MB_FAILURE, " no mesh edges");
    EntityHandle mesh_edge = mesh_edges[0]; // first one is enough
    //get its triangles; see which one is in r1 or r2; it should not be in both
    Range triangles;
    rval = _mbImpl->get_adjacencies(&mesh_edge, 1, 2, false, triangles);
    MBERRORR(rval, " can't get adjacent triangles");
    Range intx1 = intersect(triangles, r1);
    Range intx2 = intersect(triangles, r2);
    if (!intx1.empty() && !intx2.empty())
      MBERRORR(MB_FAILURE, " at least one should be empty");

    if (intx2.empty())
    {
      // it means it is in the range r1; the sense should be fine, no need to reset
      // the sense should have been fine, also
      continue;
    }
    // so it must be a triangle in r2;
    EntityHandle triangle = intx2[0];// one triangle only
    // remove the edge from boundary of face, and add it to the boundary of newFace
    // remove_parent_child( EntityHandle parent,  EntityHandle child )
    rval = _mbImpl->remove_parent_child(face, edge);
    MBERRORR(rval, " can't remove parent child relation for edge");
    // add to the new face
    rval = _mbImpl->add_parent_child(newFace, edge);
    MBERRORR(rval, " can't add parent child relation for edge");

    // set some sense, based on the sense of the mesh_edge in triangle
    int num1, sense, offset;
    rval = _mbImpl->side_number(triangle, mesh_edge, num1, sense, offset);
    MBERRORR(rval, "mesh edge not adjacent to triangle");

    rval = this->_my_geomTopoTool->set_sense(edge, newFace, sense);
    MBERRORR(rval, "can't set proper sense of edge in face");

  }

  return MB_SUCCESS;
}

ErrorCode FBEngine::set_neumann_tags(EntityHandle face, EntityHandle newFace)
{
  // these are for debugging purposes only
  // check the initial tag, then
  Tag ntag;
  ErrorCode rval = _mbImpl->tag_get_handle(NEUMANN_SET_TAG_NAME, ntag);
  MBERRORR(rval, "can't get tag handle");
  // check the value for face
  int nval;
  rval = _mbImpl->tag_get_data(ntag, &face, 1, &nval);
  if (MB_SUCCESS == rval)
  {
    nval++;
  }
  else
  {
    nval = 1;
    rval = _mbImpl->tag_set_data(ntag, &face, 1, &nval);
    MBERRORR(rval, "can't set tag");
    nval = 2;
  }
  rval = _mbImpl->tag_set_data(ntag, &newFace, 1, &nval);
  MBERRORR(rval, "can't set tag");

  return MB_SUCCESS;
}

// split the quads if needed; it will create a new gtt, which will
// contain triangles instead of quads
ErrorCode FBEngine::split_quads()
{
  // first see if we have any quads in the 2d faces
  //  _my_gsets[2] is a range of surfaces (moab sets of dimension 2)
  int num_quads=0;
  for (Range::iterator it=_my_gsets[2].begin(); it!=_my_gsets[2].end(); it++)
  {
    EntityHandle surface = *it;
    int num_q=0;
    _mbImpl->get_number_entities_by_type(surface, MBQUAD, num_q);
    num_quads+=num_q;
  }

  if (num_quads==0)
    return MB_SUCCESS; // nothing to do

  GeomTopoTool * new_gtt = _my_geomTopoTool->duplicate_model();
  if (this->_t_created)
    delete _my_geomTopoTool;

  _t_created = true; // this one is for sure created here, even if the original gtt was not

  // if we were using smart pointers, we would decrease the reference to the _my_geomTopoTool, at least
  _my_geomTopoTool = new_gtt;

  // replace the _my_gsets with the new ones, from the new set
  _my_geomTopoTool->find_geomsets(_my_gsets);

  // we have duplicated now the model, and the quads are part of the new _my_gsets[2]
  // we will split them now, by creating triangles along the smallest diagonal
  // maybe we will come up with a better criteria, but for the time being, it should be fine.
  // what we will do: we will remove all quads from the surface sets, and split them

  for (Range::iterator it2=_my_gsets[2].begin(); it2!=_my_gsets[2].end(); it2++)
  {
    EntityHandle surface = *it2;
    Range quads;
    ErrorCode rval = _mbImpl->get_entities_by_type(surface, MBQUAD, quads);
    MBERRORR(rval, "can't get quads from the surface set");
    rval = _mbImpl->remove_entities(surface, quads);
    MBERRORR(rval, "can't remove quads from the surface set"); // they are not deleted, just removed from the set
    for (Range::iterator it=quads.begin(); it!=quads.end(); it++)
    {
      EntityHandle quad = *it;
      int nnodes;
      const EntityHandle * conn;
      rval = _mbImpl->get_connectivity(quad, conn, nnodes);
      MBERRORR(rval, "can't get quad connectivity");
      // get all vertices position, to see which one is the shortest diagonal
      CartVect pos[4];
      rval = _mbImpl->get_coords(conn, 4, (double*) &pos[0]);
      MBERRORR(rval, "can't get node coordinates");
      bool diag1 = ( (pos[2]-pos[0]).length_squared() < (pos[3]-pos[1]).length_squared() );
      EntityHandle newTris[2];
      EntityHandle tri1[3]= { conn[0], conn[1], conn[2] };
      EntityHandle tri2[3]= { conn[0], conn[2], conn[3] };
      if (!diag1)
      {
        tri1[2] = conn[3];
        tri2[0] = conn[1];
      }
      rval = _mbImpl->create_element(MBTRI, tri1, 3, newTris[0]);
      MBERRORR(rval, "can't create triangle 1");
      rval = _mbImpl->create_element(MBTRI, tri2, 3, newTris[1]);
      MBERRORR(rval, "can't create triangle 2");
      rval = _mbImpl->add_entities(surface, newTris, 2);
      MBERRORR(rval, "can't add triangles to the set");
    }
    //
  }
  return MB_SUCCESS;
}

} // namespace moab


