#include <iostream>
#include <map>

#include "moab/FBEngine.hpp"
#include "moab/Interface.hpp"
#include "moab/GeomTopoTool.hpp"
#include "moab/OrientedBoxTreeTool.hpp"
#include "moab/CartVect.hpp"
#include <stdlib.h>
#include <cstring>
#include <map>
#include "assert.h"

#include "SmoothCurve.hpp"
#include "SmoothFace.hpp"

// this is just to replace MBI with moab interface, which is _mbImpl in this class
#define MBI _mbImpl
#define MBERRORR(rval, STR) { if (MB_SUCCESS != rval) { std::cout<<STR<<std::endl; return rval; } }

namespace moab {

const bool Debug_surf_eval = false;

double get_edge_length(double* p1, double* p2)
{
  return std::sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1]
      - p2[1]) + (p1[2] - p2[2]) * (p1[2] - p2[2]));
}

FBEngine::FBEngine(Interface *impl, GeomTopoTool * topoTool, const bool smooth) :
  _mbImpl(impl), _my_geomTopoTool(topoTool), _t_created(false),
      _smooth(smooth), _initialized(false), _smthFace(NULL), _smthCurve(NULL)
{
  if (!_my_geomTopoTool)
  {
    _my_geomTopoTool = new GeomTopoTool(_mbImpl);
    _t_created = true;
  }
  // should this be part of the constructor or not?
  //Init();
}
FBEngine::~FBEngine()
{
  if (_smooth)
  {
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
    _smooth = false;
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
  if (!_initialized)
  {
    if (!_my_geomTopoTool)
      return MB_FAILURE;

    ErrorCode rval = _my_geomTopoTool->find_geomsets(_my_gsets);
    assert(rval == MB_SUCCESS);

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
  for (it = _my_gsets[2].begin(); it != _my_gsets[2].end(); it++, i++)
  {
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
  for (it = _my_gsets[1].begin(); it != _my_gsets[1].end(); it++, i++)
  {
    EntityHandle curve = *it;
    _smthCurve[i] = new SmoothCurve(MBI, curve);
    _edges[curve] = _smthCurve[i];
  }

  for (i = 0; i < numSurfaces; i++)
  {
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
  for (i = 0; i < numCurves; i++)
  {
    _smthCurve[i]->compute_tangents_for_each_edge();// do we need surfaces now? or just the chains?
    // the computed edges will be marked accordingly; later one, only internal edges to surfaces are left
    _smthCurve[i]->compute_control_points_on_boundary_edges(min_dot, _faces,
        edgeCtrlTag, markTag);
  }

  // when done with boundary edges, compute the control points on all edges in the surfaces

  for (i = 0; i < numSurfaces; i++)
  {
    // also pass the tags for
    _smthFace[i]->compute_control_points_on_edges(min_dot, edgeCtrlTag, markTag);
  }

  // now we should be able to compute the control points for the facets

  for (i = 0; i < numSurfaces; i++)
  {
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

#define COPY_RANGE(r, vec) {                      \
    EntityHandle *tmp_ptr = reinterpret_cast<EntityHandle*>(vec);	\
    std::copy(r.begin(), r.end(), tmp_ptr);}

/*static inline void
 ProcessError(const char* desc);*/

double get_edge_length(double* p1, double* p2);

ErrorCode FBEngine::getRootSet(EntityHandle * root_set)
{
  *root_set = _mbImpl->get_root_set();
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
  if (0 > entity_type || 4 < entity_type)
  {
    return MB_FAILURE;
  }
  else if (entity_type < 4)
  {// 4 means all entities
    gentities = _my_gsets[entity_type];// all from root set!
  }
  else
  {
    gentities.clear();
    for (i = 0; i < 4; i++)
    {
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
  if (0 > ent_type || 3 < ent_type)
  {
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
  for (int i = 0; i < 4; i++)
  {
    if (_my_gsets[i].find(gent) != _my_gsets[i].end())
    {
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

  if (type == 0)
  {
    rval = getVtxCoord(gent, min_x, min_y, min_z);
    MBERRORR(rval, "Failed to get vertex coordinates.");
    max_x = min_x;
    max_y = min_y;
    max_z = min_z;
  }
  else if (type == 1)
  {
    rval = MB_FAILURE;
    MBERRORR(rval, "iGeom_getEntBoundBox is not supported for Edge entity type.");
  }
  else if (type == 2 || type == 3)
  {

    EntityHandle root;
    CartVect center, axis[3];
    rval = _my_geomTopoTool->get_root(gent, root);
    MBERRORR(rval, "Failed to get tree root in iGeom_getEntBoundBox.");
    rval = _my_geomTopoTool->obb_tree()->box(root, center.array(),
        axis[0].array(), axis[1].array(), axis[2].array());
    MBERRORR(rval, "Failed to get closest point in iGeom_getEntBoundBox.");

    CartVect absv[3];
    for (int i = 0; i < 3; i++)
    {
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
  }
  else
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

  if (type == 0)
  {
    rval = getVtxCoord(this_gent, on_x, on_y, on_z);
    MBERRORR(rval, "Failed to get vertex coordinates.");
  }
  else if (type == 1)
  {
    // just copy over the coordinates
    // should be modified
    *on_x = near_x;
    *on_y = near_y;
    *on_z = near_z;
  }
  else if (type == 2 || type == 3)
  {
    double point[3] = { near_x, near_y, near_z };
    double point_out[3];
    EntityHandle root, facet_out;
    if (_smooth && 2 == type)
    {
      SmoothFace* smthFace = _faces[this_gent];
      *on_x = near_x;
      *on_y = near_y;
      *on_z = near_z;
      smthFace->move_to_surface(*on_x, *on_y, *on_z);
    }
    else
    {
      rval = _my_geomTopoTool->get_root(this_gent, root);
      MBERRORR(rval, "Failed to get tree root in iGeom_getEntClosestPt.");
      rval = _my_geomTopoTool->obb_tree()->closest_to_location(point, root,
          point_out, facet_out);
      MBERRORR(rval, "Failed to get closest point in iGeom_getEntClosestPt.");

      *on_x = point_out[0];
      *on_y = point_out[1];
      *on_z = point_out[2];
    }
  }
  else
    return MB_TYPE_OUT_OF_RANGE;

  return MB_SUCCESS;
}

ErrorCode FBEngine::getVtxCoord(EntityHandle vertex_handle, double * x0,
    double * y0, double * z0)
{
  int type;
  ErrorCode rval = getEntType(vertex_handle, &type);
  MBERRORR(rval, "Failed to get entity type in getVtxCoord.");

  if (type != 0)
  {
    rval = MB_FAILURE;
    MBERRORR(rval, "Entity is not a vertex type.");
  }

  Range entities;
  rval = MBI->get_entities_by_type(vertex_handle, MBVERTEX, entities);
  MBERRORR(rval, "can't get nodes in vertex set.");

  if (entities.size() != 1)
  {
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

  if (type != 2 && type != 3)
  {
    MBERRORR(MB_FAILURE, "Entities passed into gentityNormal must be face or volume.");
  }

  if (_smooth && 2 == type)
  {
    SmoothFace* smthFace = _faces[entity_handle];
    //*on_x = near_x; *on_y = near_y; *on_z = near_z;
    smthFace-> normal_at(x, y, z, *nrml_i, *nrml_j, *nrml_k);

  }
  else
  {
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
  unsigned min_tolerace_intersections = 1000;
  double tolerance = 0.01; // TODO: how is this used ????
  //std::vector<double> distances;
  std::vector<EntityHandle> facets;
  //std::vector<EntityHandle> sets;
  unsigned int i;
  for (i = 0; i < numfaces; i++)
  {
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
  for (i = 0; i < param_coords.size(); i++)
  {
    CartVect intx = P + param_coords[i] * V;
    for (int j = 0; j < 3; j++)
      intersect_coords.push_back(intx[j]);

  }
  if (_smooth)
  {
    // correct the intersection point and the distance for smooth surfaces
    for (i = 0; i < intersect_entity_handles.size(); i++)
    {
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
  for (int i = 0; i < 4; i++)
  {
    if (_my_gsets[i].find(from) != _my_gsets[i].end())
    {
      this_dim = i;
      break;
    }
  }

  // check target dimension
  if (-1 == this_dim)
  {
    //ProcessError(iBase_FAILURE, "Entity not a geometry entity.");
    return MB_FAILURE;
  }
  else if (0 > to_dim || 3 < to_dim)
  {
    //ProcessError(iBase_FAILURE, "To dimension must be between 0 and 3.");
    return MB_FAILURE;
  }
  else if (to_dim == this_dim)
  {
    //ProcessError(iBase_FAILURE,
    //      "To dimension must be different from entity dimension.");
    return MB_FAILURE;
  }

  ErrorCode rval;
  adjs.clear();
  if (to_dim > this_dim)
  {
    int number;
    rval = MBI->num_parent_meshsets(from, &number, 0);
    rval = MBI->get_parent_meshsets(from, adjs);
    adjs.clear();
    rval = MBI->get_parent_meshsets(from, adjs, to_dim - this_dim);
  }
  else
  {
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

ErrorCode FBEngine::createTag(const char* tag_name, int tag_num_type_values,
    int tag_type, Tag & tag_handle_out)
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

  for (unsigned int i = 0; i < faces.size(); i++)
  {
    if (faces[i] == mbface)
    {
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
  for (int i = 0; i < entities_size; i++)
  {
    measures[i] = 0.;

    int type;
    EntityHandle gset=moab_entities[i];
    rval = getEntType(gset , &type);
    if (MB_SUCCESS != rval) return rval;
    if (type == 1)
    { // edge: get all edges part of the edge set
      Range entities;
      rval = MBI->get_entities_by_type(gset, MBEDGE,  entities);
      if (MB_SUCCESS != rval) return rval;



      for (Range::iterator it = entities.begin(); it!=entities.end(); it++)
      {
        EntityHandle edge=*it;
        CartVect vv[2];
        const EntityHandle *conn2=NULL;
        int  num_nodes;
        rval = MBI->get_connectivity(edge, conn2, num_nodes);
        if (MB_SUCCESS != rval || num_nodes!=2)
          return MB_FAILURE;
        rval = MBI->get_coords(conn2, 2, (double *) &(vv[0][0]) );
        if (MB_SUCCESS != rval) return rval;

        vv[0] = vv[1]-vv[0];
        measures[i] += vv[0].length();
      }
    }
    if (type == 2)
    { // surface
      // get triangles in surface; TODO: quads!
      Range entities;
      rval = MBI->get_entities_by_type(gset, MBTRI,  entities);
      if (MB_SUCCESS != rval) return rval;

      for (Range::iterator it = entities.begin(); it!=entities.end(); it++)
      {
        EntityHandle tri=*it;
        CartVect vv[3];
        const EntityHandle *conn3=NULL;
        int  num_nodes;
        rval = MBI->get_connectivity(tri, conn3, num_nodes);
        if (MB_SUCCESS != rval || num_nodes!=3)
          return MB_FAILURE;
        rval = MBI->get_coords(conn3, 3, (double *) &(vv[0][0]) );
        if (MB_SUCCESS != rval) return rval;

        vv[1] = vv[1]-vv[0];
        vv[2] = vv[2]-vv[0];
        vv[0] = vv[1] * vv[2];
        measures[i] += vv[0].length()/2;// area of triangle
      }

    }
  }
  return MB_SUCCESS;
}

ErrorCode FBEngine::getEntNrmlSense( EntityHandle face, EntityHandle region,
      int& sense )
{
  return MB_NOT_IMPLEMENTED; // not implemented
}

ErrorCode FBEngine::getEgEvalXYZ( EntityHandle edge,
                                 double x, double y, double z,
                                 double& on_x, double& on_y, double& on_z,
                                 double& tngt_i, double& tngt_j, double& tngt_k,
                                 double& cvtr_i, double& cvtr_j, double& cvtr_k )
{
  return MB_NOT_IMPLEMENTED; // not implemented
}
ErrorCode FBEngine::getFcEvalXYZ( EntityHandle face,
                                 double x, double y, double z,
                                 double& on_x, double& on_y, double& on_z,
                                 double& nrml_i, double& nrml_j, double& nrml_k,
                                 double& cvtr1_i, double& cvtr1_j, double& cvtr1_k,
                                 double& cvtr2_i, double& cvtr2_j, double& cvtr2_k )
{
  return MB_NOT_IMPLEMENTED; // not implemented
}

ErrorCode FBEngine::getEgVtxSense( EntityHandle edge, EntityHandle vtx1,
    EntityHandle vtx2, int& sense )
{
  // need to decide first or second vertex
  // important for moab
  int type;

  EntityHandle v1, v2;
  ErrorCode rval = getEntType(vtx1 , &type);
  if (MB_SUCCESS != rval || type != 0) return MB_FAILURE;
  // edge: get one vertex as part of the vertex set
  Range entities;
  rval = MBI->get_entities_by_type(vtx1, MBVERTEX,  entities);
  if (MB_SUCCESS != rval) return rval;
  if (entities.size() <1)
    return MB_FAILURE;
  v1 = entities[0]; // the first vertex
  entities.clear();
  rval = getEntType(vtx2 , &type);
  if (MB_SUCCESS != rval || type != 0) return MB_FAILURE;
  rval = MBI->get_entities_by_type(vtx2, MBVERTEX,  entities);
  if (MB_SUCCESS != rval) return rval;
  if (entities.size() <1)
    return MB_FAILURE;
  v2 = entities[0]; // the first vertex
  entities.clear();
  // now get the edges, and get the first node and the last node in sequence of edges
  // the order is important...
  // these are ordered sets !!
  std::vector<EntityHandle> ents;
  rval = MBI->get_entities_by_type(edge, MBEDGE,ents);
  if (MB_SUCCESS != rval) return rval;
  if (ents.size() <1)
    return MB_FAILURE;

  const EntityHandle* conn=NULL;
  int len;
  EntityHandle startNode, endNode;
  rval = MBI->get_connectivity(ents[0], conn, len);
  if (MB_SUCCESS != rval) return rval;
  startNode = conn[0];
  rval = MBI->get_connectivity(ents[ents.size()-1], conn, len);
  if (MB_SUCCESS != rval) return rval;

  endNode = conn[1];
  sense = 1; //
  if ( (startNode == endNode) && (v1==startNode))
  {
    sense = 0; // periodic
  }
  if ( (startNode == v1)  &&  (endNode ==v2) )
  {
    sense = 1; // forward
  }
  if ( (startNode == v2) && (endNode == v1) )
  {
    sense = -1; // reverse
  }
  return MB_SUCCESS;
}

ErrorCode FBEngine::getEntURange( EntityHandle edge,
                                 double& u_min, double& u_max )
{
  SmoothCurve * smoothCurve = _edges[edge];// this is a map
  // now, call smoothCurve methods
  smoothCurve -> get_param_range(u_min, u_max);
  return MB_SUCCESS;
}

ErrorCode FBEngine::getEntUtoXYZ( EntityHandle edge, double u,
                                 double& x, double& y, double& z )
{
  SmoothCurve * smoothCurve = _edges[edge];// this is a map
  // now, call smoothCurve methods
  smoothCurve -> position_from_u(u, x, y, z);
  return MB_SUCCESS;
}
ErrorCode FBEngine::isEntAdj( EntityHandle entity1, EntityHandle entity2,
      bool& adjacent_out )
{
  int type1, type2;
  ErrorCode rval = getEntType(entity1, &type1);
  if (MB_SUCCESS != rval) return rval;
  rval = getEntType(entity2, &type2);
  if (MB_SUCCESS != rval) return rval;

  Range adjs;
  if (type1 < type2) {
     rval = MBI->get_parent_meshsets(entity1, adjs, type2
           - type1);
     if (MB_SUCCESS != rval) return rval;// MBERRORR("Failed to get parent meshsets in iGeom_isEntAdj.");
  } else {
     rval = MBI->get_child_meshsets(entity2, adjs, type2
           - type1);
     if (MB_SUCCESS != rval) return rval;//MBERRORR("Failed to get child meshsets in iGeom_isEntAdj.");
  }

  adjacent_out = adjs.find(entity2) != _my_gsets[type2].end();

  return MB_SUCCESS;
}

} // namespace moab


