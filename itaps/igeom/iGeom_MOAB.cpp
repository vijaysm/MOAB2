#include <iostream>
#include <map>
#include "iGeom_MOAB.hpp"
#include "iMesh.h"
#include "moab/Interface.hpp"
#include "moab/GeomTopoTool.hpp"
#include "moab/OrientedBoxTreeTool.hpp"
#include "moab/CartVect.hpp"
#include <stdlib.h>
#include <cstring>
#include <map>
#include "assert.h"

using namespace moab;

iBase_Error iGeom_LAST_ERROR;
bool i_created = false; // if interface is created
bool t_created = false;
Range _my_gsets[4];
GeomTopoTool* _my_geomTopoTool = NULL;

bool debug_igeom = false;
bool Debug_surf_eval = false;

// smooth stuff
bool _smooth;
#include "SmoothFaceEval.hpp"
#include "SmoothCurveEval.hpp"

// these smooth faces and edges will be initialized after reading the file
// the maps keep the link between EH in moab (geom sets) and
//   their corresponding smooth counterparts
std::map<EntityHandle, SmoothFaceEval*> _faces;
std::map<EntityHandle, SmoothCurveEval*> _edges;
SmoothFaceEval ** _smthFace = NULL;
SmoothCurveEval ** _smthCurve = NULL;

#define COPY_RANGE(r, vec) {                      \
    EntityHandle *tmp_ptr = reinterpret_cast<EntityHandle*>(vec);	\
    std::copy(r.begin(), r.end(), tmp_ptr);}

static inline void
iGeom_processError(iBase_ErrorType code, const char* desc);

static void
iGeom_get_adjacent_entities(iGeom_Instance instance, const EntityHandle from,
      const int to_dim, Range &adj_ents, int* err);

double get_edge_length(double* p1, double* p2);

// copied from CamalPaveDriver from MeshKit;
// it is used to initialize the smoothing procedure

bool initializeSmoothing(iGeom_Instance instance) {
   //
   // first of all, we need to retrieve all the surfaces from the (root) set
   // in icesheet_test we use iGeom, but maybe that is a stretch
   // get directly the sets with geom dim 2, and from there create the SmoothFaceEval
   Tag geom_tag, gid_tag;
   ErrorCode rval = MBI->tag_get_handle(GEOM_DIMENSION_TAG_NAME, geom_tag);
   rval = MBI->tag_get_handle(GLOBAL_ID_TAG_NAME, gid_tag);

#if 0
   // traverse the model, dimension 2, 1, 0
   Range csets, fsets, vsets;
   //std::vector<moab::EntityHandle> sense_ents;
   //std::vector<int> senses, pgids;
   int dim=2;
   void *dim_ptr = &dim;
   //bool sense;

   moab::GeomTopoTool * my_geomTopoTool = new moab::GeomTopoTool(MBI);

   rval = my_geomTopoTool->construct_obb_trees();
   assert(MB_SUCCESS==rval);

   fsets.clear();
   rval = MBI->get_entities_by_type_and_tag(0, MBENTITYSET,
         &geom_tag, &dim_ptr, 1,
         fsets, 1, false);
   int numSurfaces = fsets.size();
#endif
   int numSurfaces = _my_gsets[2].size();
   //SmoothFaceEval ** smthFace = new SmoothFaceEval *[numSurfaces];
   _smthFace = new SmoothFaceEval *[numSurfaces];
   if (!t_created) {
      GETGTT(instance);
      rval = _my_geomTopoTool->construct_obb_trees();
      if (rval != MB_SUCCESS)
      { 
         std::cout << "Failed to construct obb tree.\n";
         return false;
      }
      t_created = true;
   }
   // there should also be a map from surfaces to evaluators
   //std::map<MBEntityHandle, SmoothFaceEval*> mapSurfaces;

   int i = 0;
   Range::iterator it;
   for (it = _my_gsets[2].begin(); it != _my_gsets[2].end(); it++, i++) {
      EntityHandle face = *it;
      _smthFace[i] = new SmoothFaceEval(MBI, face, _my_geomTopoTool);// geom topo tool will be used for searching,
      // among other things; also for senses in edge sets...
      _faces[face] = _smthFace[i];
   }
#if 0
   csets.clear();
   dim = 1; // get now the curves
   rval = MBI->get_entities_by_type_and_tag(0, MBENTITYSET,
         &geom_tag, &dim_ptr, 1,
         csets, 1, false);
#endif
   int numCurves = _my_gsets[1].size();//csets.size();
   //SmoothCurveEval ** smthCurve = new SmoothCurveEval *[numCurves];
   _smthCurve = new SmoothCurveEval *[numCurves];
   // there should also be a map from surfaces to evaluators
   //std::map<MBEntityHandle, SmoothCurveEval*> mapCurves;

   i = 0;
   for (it = _my_gsets[1].begin(); it != _my_gsets[1].end(); it++, i++) {
      EntityHandle curve = *it;
      _smthCurve[i] = new SmoothCurveEval(MBI, curve);
      _edges[curve] = _smthCurve[i];
   }
#if 0
   // create another mapping for vertices (sets of dimension 0)
   vsets.clear();
   dim = 0; // get now the vertice sets (dimension 0)
   rval = MBI->get_entities_by_type_and_tag(0, MBENTITYSET,
         &geom_tag, &dim_ptr, 1,
         vsets, 1, false);
   int numGeoVertices = vsets.size();
   //SmoothVertex ** smthVertex = new SmoothVertex *[numGeoVertices];
   _smthVertex = new SmoothVertex *[numGeoVertices];
   // there should also be a map from original vertex sets to new vertex sets (or just new vertex??)
   //std::map<MBEntityHandle, SmoothVertex*> mapVertices;

   i=0;
   for ( it = vsets.begin(); it!=vsets.end(); it++, i++)
   {
      MBEntityHandle vertex= *it;
      _smthVertex[i]= new SmoothVertex(MBI, vertex, _mbo);
      _mapVertices[vertex] = _smthVertex[i];
   }
#endif
   // _mb, mapSurfaces, mapCurves, mapVertices are characterizing the geometry/topology of the initial mesh

   //SmoothFaceEval geom_eval(_mb, _set);
   // initialize the smooth mesh evaluator
   //
   //geom_eval.Initialize(): it is decomposed in initializing first the gradients
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
   // the SmoothFaceEval (surface evaluator), that has everything, including the normals!!!
   assert(rval==MB_SUCCESS);

   // create the tag also for control points on the edges
   double defCtrlPoints[9] = { 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   Tag edgeCtrlTag;
   rval = MBI->tag_create("CONTROLEDGE", 9 * sizeof(double), MB_TAG_DENSE,
         edgeCtrlTag, &defCtrlPoints);
   if (MB_SUCCESS != rval)
      return false;

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
      _smthFace[i]->compute_control_points_on_edges(min_dot, edgeCtrlTag,
            markTag);
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

   return true;
}

void iGeom_getDescription(iGeom_Instance instance, char* descr, int* err,
      int descr_len) {
   unsigned int
         len =
               MIN(strlen(iGeom_LAST_ERROR.description), ((unsigned int) descr_len));
   strncpy(descr, iGeom_LAST_ERROR.description, len);
   descr[len] = '\0';
   RETURN(iBase_SUCCESS);
}

void iGeom_getErrorType(iGeom_Instance instance,
/*out*/int *error_type, int *err) {
   *error_type = iGeom_LAST_ERROR.error_type;
   RETURN(iBase_SUCCESS);
}

void iGeom_newGeom(char const* options, iGeom_Instance* instance_out, int* err,
      int options_len) {
   if (*instance_out && !(reinterpret_cast<iMesh_Instance*> (instance_out))) {
      *err = iBase_INVALID_ENTITY_TYPE;
      ERRORR("Passed in instance must be an iMesh_Instance*.");
   }

   // make a new imesh instance
   iMesh_newMesh(options, reinterpret_cast<iMesh_Instance*> (instance_out),
         err, options_len);
   ERRORR("Failure to create instance.");

   i_created = true;

   RETURN(iBase_SUCCESS);
}

void iGeom_dtor(iGeom_Instance instance, int* err) {
   if (i_created) {

      // this was supposed to help, but it crashes, still
      // we have some memory leaks if we do not delete the topo tool,
      //   the smooth evaluators, etc

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
         _smooth = false;
      }

      for (int j = 0; j < 4; j++)
         _my_gsets[j].clear();
      delete _my_geomTopoTool;
      _my_geomTopoTool = NULL;

      iMesh_dtor(IMESH_INSTANCE(instance), err);
      ERRORR("Failed to destruct instance.");
   }

   RETURN(iBase_SUCCESS);
}

void iGeom_load(iGeom_Instance instance, char const* name, char const* options,
      int* err, int name_len, int options_len) {
   // first remove option for smooth facetting
   char * smth = "SMOOTH;";
   const char * res = NULL;

   char * reducedOptions = NULL;
   int reduced_len = options_len;
   if (options)
      res = strstr(options, smth);
   if (res) {
      // extract that option, will not be recognized by our moab/imesh
      reducedOptions = new char[options_len - 6];
      int preLen = (int) (res - options);
      strncpy(reducedOptions, options, preLen);
      int postLen = options_len - 7 - preLen;

      char * tmp = reducedOptions + preLen;

      strncpy(tmp, res + 7, postLen);
      reducedOptions[options_len - 7] = 0;
      reduced_len = options_len - 7;
      std::cout << reducedOptions << std::endl;
      _smooth = true;

   } else {
      reducedOptions = const_cast<char *> (options);
   }
   // load mesh-based geometry
   iMesh_load(IMESH_INSTANCE(instance), NULL, name, reducedOptions, err,
         name_len, reduced_len);
   ERRORR("Failure to load geometry.");

   // keep mesh-based geometries in Range
   GETGTT(instance);
   ErrorCode rval = _my_geomTopoTool->find_geomsets(_my_gsets);
   MBERRORR("Failure to keep geometry list.");

   if (debug_igeom) {
      iBase_EntityHandle *entities;
      int entities_alloc, entities_size;
      for (int i = 0; i < iMesh_ALL_TOPOLOGIES; i++) {
         entities = NULL;
         entities_alloc = 0;
         iMesh_getEntities(IMESH_INSTANCE(instance), NULL, iBase_ALL_TYPES, i,
               &entities, &entities_alloc, &entities_size, err);
         ERRORR("Failed to get entities\n");
         std::cout << "type_geom=" << i << ", number=" << entities_size
               << std::endl;
      }

      iBase_EntitySetHandle *esets = NULL;
      int esets_alloc = 0, esets_size;
      iMesh_getEntSets(IMESH_INSTANCE(instance), NULL, 1, &esets, &esets_alloc,
            &esets_size, err);
      ERRORR("Failed to get entity sets\n");
      std::cout << "entity_geom set number=" << esets_size << std::endl;

      entities = NULL;
      entities_alloc = 0;
      iMesh_getEntities(IMESH_INSTANCE(instance),
            reinterpret_cast<iBase_EntitySetHandle> (*(_my_gsets[0].begin())),
            iBase_ALL_TYPES, 0, &entities, &entities_alloc, &entities_size, err);
      ERRORR("Failed to get entities\n");

      double x, y, z;
      iMesh_getVtxCoord(IMESH_INSTANCE(instance), entities[0], &x, &y, &z, err);
      std::cout << "vertex coords=" << x << ", " << y << ", " << z << std::endl;
   }

   if (_smooth)
      initializeSmoothing(instance);

   RETURN(iBase_SUCCESS);
}

void iGeom_save(iGeom_Instance instance, char const* name, char const* options,
      int* err, int name_len, int options_len) {
   iMesh_save(IMESH_INSTANCE(instance), NULL, name, options, err, name_len,
         options_len);
   ERRORR("Failed get save geomtry instance.");

   RETURN(iBase_SUCCESS);
}

void iGeom_getRootSet(iGeom_Instance instance, iBase_EntitySetHandle* root_set,
      int* err) {
   iMesh_getRootSet(IMESH_INSTANCE(instance), root_set, err);
   ERRORR("Failed get root set.");

   RETURN(iBase_SUCCESS);
}

void iGeom_getBoundBox(iGeom_Instance, double* min_x, double* min_y,
      double* min_z, double* max_x, double* max_y, double* max_z, int* err) {
   RETURN(iBase_NOT_SUPPORTED);
}

void iGeom_getEntities(iGeom_Instance instance,
      iBase_EntitySetHandle set_handle, int entity_type,
      iBase_EntityHandle** entity_handles, int* entity_handles_allocated,
      int* entity_handles_size, int* err) {
   int i;
   if (0 > entity_type || 4 < entity_type) {
      *err = iBase_INVALID_ENTITY_TYPE;
      ERRORR("Bad entity type.");
   } else if (entity_type < 4) {
      *entity_handles_size = _my_gsets[entity_type].size();
      CHECK_SIZE(*entity_handles, *entity_handles_allocated,
            *entity_handles_size, iBase_EntityHandle, NULL);
      COPY_RANGE(_my_gsets[entity_type], *entity_handles);
   } else {
      *entity_handles_size = 0;
      Range total_range;
      for (i = 0; i < 4; i++) {
         total_range.merge(_my_gsets[i]);
      }
      *entity_handles_size = total_range.size();
      CHECK_SIZE(*entity_handles, *entity_handles_allocated,
            *entity_handles_size, iBase_EntityHandle, NULL);
      COPY_RANGE(total_range, *entity_handles);
   }

   RETURN(iBase_SUCCESS);
}

void iGeom_getNumOfType(iGeom_Instance instance,
      iBase_EntitySetHandle set_handle, int entity_type, int* num_out, int* err) {
   if (0 > entity_type || 3 < entity_type) {
      *err = iBase_INVALID_ENTITY_TYPE;
      ERRORR("Bad entity type.");
   }
   *num_out = _my_gsets[entity_type].size();

   RETURN(iBase_SUCCESS);
}

void iGeom_getEntType(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, int* type, int* err) {
   for (int i = 0; i < 4; i++) {
      if (_my_gsets[i].find(MBH_cast(entity_handle)) != _my_gsets[i].end()) {
         *type = i;
         RETURN(iBase_SUCCESS);
      }
   }

   *err = iBase_INVALID_ENTITY_TYPE;
   ERRORR("Entity not a geometry entity.");
}

void iGeom_getArrType(iGeom_Instance instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int** type, int* type_allocated, int* type_size, int* err) {
   CHECK_SIZE(*type, *type_allocated, *type_size, int, NULL);

   int tmp_err;
   *err = iBase_SUCCESS;

   for (int i = 0; i < entity_handles_size; i++) {
      iGeom_getEntType(instance, entity_handles[i], *type + i, &tmp_err);
      if (iBase_SUCCESS != tmp_err) {
         *err = tmp_err;
         ERRORR("Failed to get entity type in iGeom_getArrType.");
      }
   }

   RETURN(iBase_SUCCESS);
}

void iGeom_getEntAdj(iGeom_Instance instance, iBase_EntityHandle entity_handle,
      int to_dimension, iBase_EntityHandle** adj_entities,
      int* adj_entities_allocated, int* adj_entities_size, int* err) {
   Range adjs;
   EntityHandle this_ent = MBH_cast(entity_handle);

   // get adjacent
   iGeom_get_adjacent_entities(instance, this_ent, to_dimension, adjs, err);
   ERRORR("Failed to get adjacent entities in iGeom_getEntAdj.");

   // copy adjacent entities
   *adj_entities_size = adjs.size();
   CHECK_SIZE(*adj_entities, *adj_entities_allocated,
         *adj_entities_size, iBase_EntityHandle, NULL);
   COPY_RANGE(adjs, *adj_entities);

   RETURN(iBase_SUCCESS);
}

void iGeom_getArrAdj(iGeom_Instance instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int requested_entity_type, iBase_EntityHandle** adj_entity_handles,
      int* adj_entity_handles_allocated, int* adj_entity_handles_size,
      int** offset, int* offset_allocated, int* offset_size, int* err) {
   // check offset array size
   Range temp_range, total_range;
   CHECK_SIZE(*offset, *offset_allocated, entity_handles_size + 1, int, NULL);
   *offset_size = entity_handles_size + 1;

   // get adjacent entities
   for (int i = 0; i < entity_handles_size; ++i) {
      (*offset)[i] = total_range.size();
      temp_range.clear();
      iGeom_get_adjacent_entities(instance, MBH_cast(entity_handles[i]),
            requested_entity_type, temp_range, err);
      ERRORR("Failed to get adjacent entities in iGeom_getEntAdj.");

      total_range.merge(temp_range);
   }
   int nTot = total_range.size();
   (*offset)[entity_handles_size] = nTot;

   // copy adjacent entities
   CHECK_SIZE(*adj_entity_handles, *adj_entity_handles_allocated,
         nTot, iBase_EntityHandle, NULL);
   COPY_RANGE(total_range, *adj_entity_handles);
   *adj_entity_handles_size = nTot;

   RETURN(iBase_SUCCESS);
}

void iGeom_getEnt2ndAdj(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, int bridge_dimension, int to_dimension,
      iBase_EntityHandle** adjacent_entities, int* adjacent_entities_allocated,
      int* adjacent_entities_size, int* err) {
   Range to_ents, bridge_ents, tmp_ents;
   iGeom_get_adjacent_entities(instance, MBH_cast(entity_handle),
         bridge_dimension, bridge_ents, err);
   ERRORR("Failed to get adjacent entities in iGeom_getEnt2ndAdj.");

   Range::iterator iter, jter, kter, end_jter;
   Range::iterator end_iter = bridge_ents.end();
   for (iter = bridge_ents.begin(); iter != end_iter; iter++) {
      iGeom_get_adjacent_entities(instance, *iter, to_dimension, tmp_ents, err);
      ERRORR("Failed to get adjacent entities in iGeom_getEnt2ndAdj.");

      for (jter = tmp_ents.begin(); jter != end_jter; jter++) {
         if (to_ents.find(*jter) == to_ents.end()) {
            to_ents.insert(*jter);
         }
      }
      tmp_ents.clear();
   }

   *adjacent_entities_size = to_ents.size();
   CHECK_SIZE(*adjacent_entities, *adjacent_entities_allocated,
         *adjacent_entities_size, iBase_EntityHandle, NULL);
   COPY_RANGE(to_ents, *adjacent_entities);

   RETURN(iBase_SUCCESS);
}

void iGeom_getArr2ndAdj(iGeom_Instance instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int order_adjacent_key, int requested_entity_type,
      iBase_EntityHandle** adj_entity_handles,
      int* adj_entity_handles_allocated, int* adj_entity_handles_size,
      int** offset, int* offset_allocated, int* offset_size, int* err) {
   CHECK_SIZE(*offset, *offset_allocated, entity_handles_size + 1, int, NULL);
   Range bridge_range, temp_range, entity_range, total_range;

   for (int i = 0; i < entity_handles_size; ++i) {
      bridge_range.clear();
      entity_range.clear();
      iGeom_get_adjacent_entities(instance, MBH_cast(entity_handles[i]),
            order_adjacent_key, bridge_range, err);
      ERRORR("Failed to get adjacent entities in iGeom_getArr2ndAdj.");

      Range::iterator iter, jter, end_jter;
      Range::iterator end_iter = bridge_range.end();
      for (iter = bridge_range.begin(); iter != end_iter; iter++) {
         temp_range.clear();
         iGeom_get_adjacent_entities(instance, *iter, requested_entity_type,
               temp_range, err);
         ERRORR("Failed to get adjacent entities in iGeom_getArr2ndAdj.");

         for (jter = temp_range.begin(); jter != end_jter; jter++) {
            if (entity_range.find(*jter) == entity_range.end()) {
               entity_range.insert(*jter);
            }
         }
      }

      (*offset)[i] = total_range.size();
      total_range.merge(entity_range);
   }
   *adj_entity_handles_size = total_range.size();
   (*offset)[entity_handles_size] = *adj_entity_handles_size;

   CHECK_SIZE(*adj_entity_handles, *adj_entity_handles_allocated,
         *adj_entity_handles_size, iBase_EntityHandle, NULL);
   COPY_RANGE(total_range, *adj_entity_handles);

   RETURN(iBase_SUCCESS);
}

void iGeom_isEntAdj(iGeom_Instance instance, iBase_EntityHandle entity_handle1,
      iBase_EntityHandle entity_handle2, int* are_adjacent, int* err) {
   int type1, type2;
   iGeom_getEntType(instance, entity_handle1, &type1, err);
   ERRORR("Failed to get entity type in iGeom_isEntAdj.");
   iGeom_getEntType(instance, entity_handle2, &type2, err);
   ERRORR("Failed to get entity type in iGeom_isEntAdj.");

   ErrorCode rval;
   Range adjs;
   if (type1 < type2) {
      rval = MBI->get_parent_meshsets(MBH_cast(entity_handle1), adjs, type2
            - type1);
      MBERRORR("Failed to get parent meshsets in iGeom_isEntAdj.");
   } else {
      rval = MBI->get_child_meshsets(MBH_cast(entity_handle1), adjs, type2
            - type1);
      MBERRORR("Failed to get child meshsets in iGeom_isEntAdj.");
   }

   *are_adjacent = adjs.find(MBH_cast(entity_handle2))
         != _my_gsets[type2].end();

   RETURN(iBase_SUCCESS);
}

void iGeom_isArrAdj(iGeom_Instance instance,
      iBase_EntityHandle const* entity_handles_1, int entity_handles_1_size,
      iBase_EntityHandle const* entity_handles_2, int entity_handles_2_size,
      int** is_adjacent_info, int* is_adjacent_info_allocated,
      int* is_adjacent_info_size, int* err) {
   int index1 = 0;
   int index2 = 0;
   size_t index1_step, index2_step;
   int count;

   // If either list contains only 1 entry, compare that entry with
   // every entry in the other list.
   if (entity_handles_1_size == entity_handles_2_size) {
      index1_step = index2_step = 1;
      count = entity_handles_1_size;
   } else if (entity_handles_1_size == 1) {
      index1_step = 0;
      index2_step = 1;
      count = entity_handles_2_size;
   } else if (entity_handles_2_size == 1) {
      index1_step = 1;
      index2_step = 0;
      count = entity_handles_1_size;
   } else {
      RETURN(iBase_INVALID_ENTITY_COUNT);
   }

   CHECK_SIZE(*is_adjacent_info, *is_adjacent_info_allocated,
         count, int, NULL);

   for (int i = 0; i < count; ++i) {
      iGeom_isEntAdj(instance, entity_handles_1[index1],
            entity_handles_2[index2], &((*is_adjacent_info)[i]), err);
      ERRORR("Failed to check if entities are adjacent in iGeom_isArrAdj.");

      index1 += index1_step;
      index2 += index2_step;
   }

   RETURN(iBase_SUCCESS);
}

void iGeom_getEntClosestPt(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, double near_x, double near_y,
      double near_z, double* on_x, double* on_y, double* on_z, int* err) {
   ErrorCode rval;
   int type;
   iGeom_getEntType(instance, entity_handle, &type, err);
   ERRORR("Failed to get entity type.");

   if (type == 0) {
      iGeom_getVtxCoord(instance, entity_handle, on_x, on_y, on_z, err);
      ERRORR("Failed to get vertex coordinates.");
   } else if (type == 1) {
      // just copy over the coordinates
      // should be modified
      *on_x = near_x;
      *on_y = near_y;
      *on_z = near_z;
   } else if (type == 2 || type == 3) {
      if (!t_created) {
         GETGTT(instance);
         rval = _my_geomTopoTool->construct_obb_trees();
         MBERRORR("Failed to construct obb tree.");
         t_created = true;
      }

      double point[3] = { near_x, near_y, near_z };
      double point_out[3];
      EntityHandle root, facet_out;
      if (_smooth && 2 == type) {
         EntityHandle geoSet = MBH_cast(entity_handle);
         SmoothFaceEval* smthFace = _faces[geoSet];
         *on_x = near_x;
         *on_y = near_y;
         *on_z = near_z;
         smthFace->move_to_surface(*on_x, *on_y, *on_z);

      } else {
         rval = _my_geomTopoTool->get_root(MBH_cast(entity_handle), root);
         MBERRORR("Failed to get tree root in iGeom_getEntClosestPt.");
         rval = _my_geomTopoTool->obb_tree()->closest_to_location(point, root,
               point_out, facet_out);
         MBERRORR("Failed to get closest point in iGeom_getEntClosestPt.");

         *on_x = point_out[0];
         *on_y = point_out[1];
         *on_z = point_out[2];
      }
   } else
      RETURN(iBase_INVALID_ENTITY_TYPE);

   RETURN(iBase_SUCCESS);
}

void iGeom_getArrClosestPt(iGeom_Instance instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int storage_order, double const* near_coordinates,
      int near_coordinates_size, double** on_coordinates,
      int* on_coordinates_allocated, int* on_coordinates_size, int* err) {
   CHECK_SIZE(*on_coordinates, *on_coordinates_allocated,
         near_coordinates_size, double, NULL);
   for (int i = 0; i < entity_handles_size; i++) {
      if (storage_order == iBase_INTERLEAVED) {
         iGeom_getEntClosestPt(instance, entity_handles[i], near_coordinates[3
               * i], near_coordinates[3 * i + 1], near_coordinates[3 * i + 2],
               on_coordinates[3 * i], on_coordinates[3 * i + 1],
               on_coordinates[3 * i + 2], err);
      } else if (storage_order == iBase_BLOCKED) {
         iGeom_getEntClosestPt(instance, entity_handles[i],
               near_coordinates[i], near_coordinates[i + entity_handles_size],
               near_coordinates[i + 2 * entity_handles_size],
               on_coordinates[i], on_coordinates[i + entity_handles_size],
               on_coordinates[i + 2 * entity_handles_size], err);
      }
      ERRORR("Failed to get closest point.");
   }

   RETURN(iBase_SUCCESS);
}

void iGeom_getEntNrmlXYZ(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, double x, double y, double z,
      double* nrml_i, double* nrml_j, double* nrml_k, int* err) {
   // just do for surface and volume
   int type;
   iGeom_getEntType(instance, entity_handle, &type, err);
   ERRORR("Failed to get entity type in iGeom_getEntNrmlXYZ.");

   if (type != 2 && type != 3) {
      *err = iBase_INVALID_ENTITY_TYPE;
      ERRORR("Entities passed into gentityNormal must be face or volume.");
   }

   if (_smooth && 2 == type) {
      EntityHandle geoSet = MBH_cast(entity_handle);
      SmoothFaceEval* smthFace = _faces[geoSet];
      //*on_x = near_x; *on_y = near_y; *on_z = near_z;
      smthFace-> normal_at(x, y, z, *nrml_i, *nrml_j, *nrml_k);

   } else {
      // get closest location and facet
      double point[3] = { x, y, z };
      double point_out[3];
      EntityHandle root, facet_out;
      _my_geomTopoTool->get_root(MBH_cast(entity_handle), root);
      ErrorCode rval = _my_geomTopoTool->obb_tree()->closest_to_location(point,
            root, point_out, facet_out);
      MBERRORR("Failed to get closest location in iGeom_getEntNrmlXYZ.");

      // get facet normal
      const EntityHandle* conn;
      int len;
      CartVect coords[3], normal;
      rval = MBI->get_connectivity(facet_out, conn, len);
      MBERRORR("Failed to get triangle connectivity in iGeom_getEntNrmlXYZ.");
      if (len != 3)
         RETURN(iBase_FAILURE);

      rval = MBI->get_coords(conn, len, coords[0].array());
      MBERRORR("Failed to get triangle coordinates in iGeom_getEntNrmlXYZ.");

      coords[1] -= coords[0];
      coords[2] -= coords[0];
      normal = coords[1] * coords[2];
      normal.normalize();
      *nrml_i = normal[0];
      *nrml_j = normal[1];
      *nrml_k = normal[2];
   }
   *err = 0;
   RETURN(iBase_SUCCESS);
}

void iGeom_getArrNrmlXYZ(iGeom_Instance instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int storage_order, double const* coordinates, int coordinates_size,
      double** normals, int* normals_allocated, int* normals_size, int* err) {
   // set up iteration according to storage order.
   // allow either gentity_handles or near_coordinates to contain
   // only one value, where that single value is applied for every
   // entry in the other list.
   size_t index = 0;
   size_t coord_step, norm_step = 1, ent_step;
   int count;
   if (3 * entity_handles_size == coordinates_size) {
      coord_step = ent_step = 1;
      count = entity_handles_size;
   } else if (coordinates_size == 3) {
      coord_step = 0;
      ent_step = 1;
      count = entity_handles_size;
   } else if (entity_handles_size == 1) {
      coord_step = 1;
      ent_step = 0;
      count = coordinates_size / 3;
   } else {
      *err = iBase_INVALID_ENTITY_COUNT;
      ERRORR("Mismatched array sizes");
   }

   // check or pre-allocate the coordinate arrays
   CHECK_SIZE(*normals, *normals_allocated, 3*count, double, NULL);

   const double *coord_x, *coord_y, *coord_z;
   double *norm_x, *norm_y, *norm_z;
   if (storage_order == iBase_BLOCKED) {
      coord_x = coordinates;
      coord_y = coord_x + coordinates_size / 3;
      coord_z = coord_y + coordinates_size / 3;
      norm_x = *normals;
      norm_y = norm_x + count;
      norm_z = norm_y + count;
      norm_step = 1;
   } else {
      storage_order = iBase_INTERLEAVED; /* set if unspecified */
      coord_x = coordinates;
      coord_y = coord_x + 1;
      coord_z = coord_x + 2;
      norm_x = *normals;
      norm_y = norm_x + 1;
      norm_z = norm_x + 2;
      coord_step *= 3;
      norm_step = 3;
   }

   for (int i = 0; i < count; ++i) {
      iGeom_getEntNrmlXYZ(instance, entity_handles[index], *coord_x, *coord_y,
            *coord_z, norm_x, norm_y, norm_z, err);
      ERRORR("Failed to get entity normal of point.");

      index += ent_step;
      coord_x += coord_step;
      coord_y += coord_step;
      coord_z += coord_step;
      norm_x += norm_step;
      norm_y += norm_step;
      norm_z += norm_step;
   }

   RETURN(iBase_SUCCESS);
}

void iGeom_getEntNrmlPlXYZ(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, double x, double y, double z,
      double* pt_x, double* pt_y, double* pt_z, double* nrml_i, double* nrml_j,
      double* nrml_k, int* err) {
   // just do for surface and volume
   int type;
   iGeom_getEntType(instance, entity_handle, &type, err);
   ERRORR("Failed to get entity type in iGeom_getEntNrmlPlXYZ.");

   if (type != 2 && type != 3) {
      *err = iBase_INVALID_ENTITY_TYPE;
      ERRORR("Entities passed into gentityNormal must be face or volume.");
   }

   // get closest location and facet
   double point[3] = { x, y, z };
   double point_out[3];
   EntityHandle root, facet_out;
   _my_geomTopoTool->get_root(MBH_cast(entity_handle), root);
   ErrorCode rval = _my_geomTopoTool->obb_tree()->closest_to_location(point,
         root, point_out, facet_out);
   MBERRORR("Failed to get closest location in iGeom_getEntNrmlPlXYZ.");

   // get closest point
   *pt_x = point_out[0];
   *pt_y = point_out[1];
   *pt_z = point_out[2];

   // get facet normal
   const EntityHandle* conn;
   int len;
   CartVect coords[3], normal;
   rval = MBI->get_connectivity(facet_out, conn, len);
   MBERRORR("Failed to get triangle connectivity in iGeom_getEntNrmlPlXYZ.");
   if (len != 3)
      RETURN(iBase_FAILURE);

   rval = MBI->get_coords(conn, len, coords[0].array());
   MBERRORR("Failed to get triangle coordinates in iGeom_getEntNrmlPlXYZ.");

   coords[1] -= coords[0];
   coords[2] -= coords[0];
   normal = coords[1] * coords[2];
   normal.normalize();
   *nrml_i = normal[0];
   *nrml_j = normal[1];
   *nrml_k = normal[2];

   RETURN(iBase_SUCCESS);
}

void iGeom_getArrNrmlPlXYZ(iGeom_Instance instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int storage_order, double const* near_coordinates,
      int near_coordinates_size, double** on_coordinates,
      int* on_coordinates_allocated, int* on_coordinates_size,
      double** normals, int* normals_allocated, int* normals_size, int* err) {
   // set up iteration according to storage order.
   // allow either gentity_handles or near_coordinates to contain
   // only one value, where that single value is applied for every
   // entry in the other list.
   size_t index = 0;
   size_t near_step, on_step = 1, ent_step;
   int count;
   if (3 * entity_handles_size == near_coordinates_size) {
      near_step = ent_step = 1;
      count = entity_handles_size;
   } else if (near_coordinates_size == 3) {
      near_step = 0;
      ent_step = 1;
      count = entity_handles_size;
   } else if (entity_handles_size == 1) {
      near_step = 1;
      ent_step = 0;
      count = near_coordinates_size / 3;
   } else {
      *err = iBase_INVALID_ENTITY_COUNT;
      ERRORR("Mismatched array sizes");
   }

   // check or pre-allocate the coordinate arrays
   CHECK_SIZE(*on_coordinates, *on_coordinates_allocated, 3*count, double, NULL);
   CHECK_SIZE(*normals, *normals_allocated, 3*count, double, NULL);

   const double *near_x, *near_y, *near_z;
   double *on_x, *on_y, *on_z;
   double *norm_x, *norm_y, *norm_z;
   if (storage_order == iBase_BLOCKED) {
      near_x = near_coordinates;
      near_y = near_x + near_coordinates_size / 3;
      near_z = near_y + near_coordinates_size / 3;
      on_x = *on_coordinates;
      on_y = on_x + count;
      on_z = on_y + count;
      norm_x = *normals;
      norm_y = norm_x + count;
      norm_z = norm_y + count;
      on_step = 1;
   } else {
      storage_order = iBase_INTERLEAVED; /* set if unspecified */
      near_x = near_coordinates;
      near_y = near_x + 1;
      near_z = near_x + 2;
      on_x = *on_coordinates;
      on_y = on_x + 1;
      on_z = on_x + 2;
      norm_x = *normals;
      norm_y = norm_x + 1;
      norm_z = norm_x + 2;
      near_step *= 3;
      on_step = 3;
   }

   for (int i = 0; i < count; ++i) {
      iGeom_getEntNrmlPlXYZ(instance, entity_handles[index], *near_x, *near_y,
            *near_z, on_x, on_y, on_z, norm_x, norm_y, norm_z, err);
      ERRORR("Failed to get entity normal of point.");

      //entities += ent_step;
      index += ent_step;
      near_x += near_step;
      near_y += near_step;
      near_z += near_step;
      on_x += on_step;
      on_y += on_step;
      on_z += on_step;
      norm_x += on_step;
      norm_y += on_step;
      norm_z += on_step;
   }

   RETURN(iBase_SUCCESS);
}

void iGeom_getEntTgntXYZ(iGeom_Instance, iBase_EntityHandle entity_handle,
      double x, double y, double z, double* tgnt_i, double* tgnt_j,
      double* tgnt_k, int* err) {
   RETURN(iBase_NOT_SUPPORTED);
}

void iGeom_getArrTgntXYZ(iGeom_Instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int storage_order, double const* coordinates, int coordinates_size,
      double** tangents, int* tangents_allocated, int* tangents_size, int* err) {
   RETURN(iBase_NOT_SUPPORTED);
}

void iGeom_getEntBoundBox(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, double* min_x, double* min_y,
      double* min_z, double* max_x, double* max_y, double* max_z, int* err) {
   ErrorCode rval;
   int type;
   iGeom_getEntType(instance, entity_handle, &type, err);
   ERRORR("Failed to get entity type.");

   if (type == 0) {
      iGeom_getVtxCoord(instance, entity_handle, min_x, min_y, min_z, err);
      ERRORR("Failed to get vertex coordinates.");
      max_x = min_x;
      max_y = min_y;
      max_z = min_z;
   } else if (type == 1) {
      *err = iBase_NOT_SUPPORTED;
      ERRORR("iGeom_getEntBoundBox is not supported for Edge entity type.");
   } else if (type == 2 || type == 3) {
      if (!t_created) {
         GETGTT(instance);
         rval = _my_geomTopoTool->construct_obb_trees();
         MBERRORR("Failed to construct obb tree.");
         t_created = true;
      }

      EntityHandle root;
      CartVect center, axis[3];
      rval = _my_geomTopoTool->get_root(MBH_cast(entity_handle), root);
      MBERRORR("Failed to get tree root in iGeom_getEntBoundBox.");
      rval = _my_geomTopoTool->obb_tree()->box(root, center.array(),
            axis[0].array(), axis[1].array(), axis[2].array());
      MBERRORR("Failed to get closest point in iGeom_getEntBoundBox.");

      CartVect absv[3];
      for (int i=0; i<3; i++)
      {
         absv[i]= CartVect( fabs(axis[i][0]), fabs(axis[i][1]), fabs(axis[i][2]) );
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
      RETURN(iBase_INVALID_ENTITY_TYPE);

   RETURN(iBase_SUCCESS);
}

void iGeom_getArrBoundBox(iGeom_Instance instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int storage_order, double** min_corner, int* min_corner_allocated,
      int* min_corner_size, double** max_corner, int* max_corner_allocated,
      int* max_corner_size, int* err) {
   // check or pre-allocate the coordinate arrays
   CHECK_SIZE(*min_corner, *min_corner_allocated, 3*entity_handles_size, double, NULL);
   CHECK_SIZE(*max_corner, *max_corner_allocated, 3*entity_handles_size, double, NULL);

   size_t step, init;
   if (storage_order == iBase_BLOCKED) {
      step = 1;
      init = entity_handles_size;
   } else {
      step = 3;
      init = 1;
   }
   double *min_x, *min_y, *min_z, *max_x, *max_y, *max_z;
   min_x = *min_corner;
   max_x = *max_corner;
   min_y = min_x + init;
   max_y = max_x + init;
   min_z = min_y + init;
   max_z = max_y + init;

   for (int i = 0; i < entity_handles_size; ++i) {
      iGeom_getEntBoundBox(instance, entity_handles[i], min_x, min_y, min_z,
            max_x, max_y, max_z, err);
      ERRORR("Failed to get entity bounding box in iGeom_getArrBoundBox.");

      min_x += step;
      max_x += step;
      min_y += step;
      max_y += step;
      min_z += step;
      max_z += step;
   }

   RETURN(iBase_SUCCESS);
}

void iGeom_getVtxCoord(iGeom_Instance instance,
      iBase_EntityHandle vertex_handle, double* x, double* y, double* z,
      int* err) {
   int type;
   iGeom_getEntType(instance, vertex_handle, &type, err);
   ERRORR("Failed to get entity type in iGeom_getVtxCoord.");

   if (type != 0) {
      *err = iBase_INVALID_ENTITY_TYPE;
      ERRORR("Entity is not a vertex type.");
   }

   iBase_EntityHandle *verts = NULL;
   int verts_alloc = 0, verts_size;
   iMesh_getEntities(IMESH_INSTANCE(instance),
         reinterpret_cast<iBase_EntitySetHandle> (vertex_handle),
         iBase_ALL_TYPES, iMesh_POINT, &verts, &verts_alloc, &verts_size, err);
   ERRORR("Failed to get vertices.");

   if (verts_size != 1) {
      *err = iBase_FAILURE;
      ERRORR("Vertex has multiple points.");
   }

   iMesh_getVtxCoord(IMESH_INSTANCE(instance), verts[0], x, y, z, err);
   ERRORR("Failed to get vertex coordinate.");

   RETURN(iBase_SUCCESS);
}

void iGeom_getVtxArrCoords(iGeom_Instance instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int storage_order, double** coordinates, int* coordinates_allocated,
      int* coordinates_size, int* err) {
   // check or pre-allocate the coordinate arrays
   CHECK_SIZE(*coordinates, *coordinates_allocated, 3*entity_handles_size, double, NULL);

   double *x, *y, *z;
   size_t step;
   if (storage_order == iBase_BLOCKED) {
      x = *coordinates;
      y = x + entity_handles_size;
      z = y + entity_handles_size;
      step = 1;
   } else {
      x = *coordinates;
      y = x + 1;
      z = x + 2;
      step = 3;
   }

   for (int i = 0; i < entity_handles_size; i++) {
      iGeom_getVtxCoord(instance, entity_handles[i], x, y, z, err);
      x += step;
      y += step;
      z += step;
   }

   RETURN(iBase_SUCCESS);
}

void iGeom_getPntRayIntsct(iGeom_Instance instance, double x, double y, double z,
      double dir_x, double dir_y, double dir_z,
      iBase_EntityHandle** intersect_entity_handles,
      int* intersect_entity_handles_allocated,
      int* intersect_entity_handles_size, int storage_order,
      double** intersect_coords, int* intersect_coords_allocated,
      int* intersect_coords_size, double** param_coords,
      int* param_coords_allocated, int* param_coords_size, int* err) {
   // this is pretty cool
   // we will return only surfaces
   //
   ErrorCode rval;
   if (!t_created) {
      GETGTT(instance);
      ErrorCode rval = _my_geomTopoTool->construct_obb_trees();
      MBERRORR("Failed to construct obb tree.");
      t_created = true;
   }
   // OrientedBoxTreeTool
   // rval = _my_geomTopoTool->obb_tree()->closest_to_location??
 /*  ray_intersect_sets( std::vector<double>& distances_out,
   std::vector<EntityHandle>& sets_out,
   std::vector<EntityHandle>& facets_out,
   EntityHandle root_set,
   double tolerance,
   unsigned min_tolerace_intersections,
   const double ray_point[3],
   const double unit_ray_dir[3],
   const double* ray_length = 0,
   TrvStats* accum = 0 );*/
   //ErrorCode rval =  _my_geomTopoTool->obb_tree()->

   // loop through all faces
   // some tolerances need to be set in advance...
   int numfaces = _my_gsets[2].size();
   // do ray fire
   const double point[] = {x, y, z};
   const double dir[] = {dir_x, dir_y, dir_z};
   CartVect P(point);
   CartVect V(dir);
   unsigned min_tolerace_intersections = 1000;
   double tolerance =0.01; // TODO: how is this used ????
   std::vector<double> distances;
   std::vector<EntityHandle> facets;
   std::vector<EntityHandle> sets;
   int i;
   for (i=0; i<numfaces; i++)
   {
      EntityHandle face = _my_gsets[2][i];
      EntityHandle rootForFace ;
      rval = _my_geomTopoTool->get_root(face, rootForFace) ;
      MBERRORR("Failed to get root of face.");
      std::vector<double> distances_out;
      std::vector<EntityHandle> sets_out;
      std::vector<EntityHandle> facets_out;
      rval = _my_geomTopoTool->obb_tree()-> ray_intersect_sets(
          distances_out,
          sets_out,
          facets_out,
          rootForFace,
          tolerance,
          min_tolerace_intersections,
          point,
          dir);
      unsigned int j;
      for (j=0; j<distances_out.size(); j++)
         distances.push_back(distances_out[j]);
      for (j=0; j<sets_out.size(); j++)
         sets.push_back(sets_out[j]);
      for (j=0; j<facets_out.size(); j++)
         facets.push_back(facets_out[j]);

      MBERRORR("Failed to get ray intersections.");
   }
   // we will say that those lists are the same size, always
   *param_coords_allocated = *param_coords_size = distances.size();
   *param_coords = (double*)malloc((*param_coords_size)*sizeof(double));

   *intersect_coords_allocated = *intersect_coords_size = 3 * distances.size();
   *intersect_coords = (double*)malloc((*intersect_coords_size)*sizeof(double));

   *intersect_entity_handles_allocated = *intersect_entity_handles_size = sets.size();
   *intersect_entity_handles = (iBase_EntityHandle*)malloc((*intersect_entity_handles_size)*sizeof(iBase_EntityHandle));

   // facets.size == distances.size()!!
   for (i=0; i<*param_coords_allocated; i++ )
   {
      (*param_coords)[i] = distances[i];
      CartVect intx = P+distances[i]*V;
      for (int j=0; j<3; j++)
         (*intersect_coords)[3*i+j] = intx[j];
      (*intersect_entity_handles)[i] = (iBase_EntityHandle)sets[i];
   }
   if (_smooth)
   {
      // correct the intersection point and the distance for smooth surfaces
      for (i=0; i<sets.size(); i++)
      {
         //EntityHandle geoSet = MBH_cast(sets[i]);
         SmoothFaceEval* sFace = _faces[sets[i]];
         // correct coordinates and distance from point
         /*moab::ErrorCode ray_intersection_correct(moab::EntityHandle facet, // (IN) the facet where the patch is defined
                  moab::CartVect &pt, // (IN) shoot from
                  moab::CartVect &ray, // (IN) ray direction
                  moab::CartVect &eval_pt, // (INOUT) The intersection point
                  double & distance, // (IN OUT) the new distance
                  bool &outside);*/
         CartVect pos(&(*intersect_coords)[3*i]);
         double dist = (*param_coords)[i];
         bool outside = false;
         rval = sFace->ray_intersection_correct(facets[i],
               P, V, pos, dist, outside);
         MBERRORR("Failed to get better point on ray.");
         (*param_coords)[i] = dist;

         for (int j=0; j<3; j++)
            (*intersect_coords)[3*i+j] = pos[j];
      }
   }
   RETURN(iBase_SUCCESS);;
}

void iGeom_getPntArrRayIntsct(iGeom_Instance, int storage_order,
      const double* coords, int coords_size, const double* directions,
      int directions_size, iBase_EntityHandle** intersect_entity_handles,
      int* intersect_entity_handles_allocated,
      int* intersect_entity_hangles_size, int** offset, int* offset_allocated,
      int* offset_size, double** intersect_coords,
      int* intersect_coords_allocated, int* intersect_coords_size,
      double** param_coords, int* param_coords_allocated,
      int* param_coords_size, int* err) {
}
void iGeom_getEntNrmlSense(iGeom_Instance, iBase_EntityHandle face,
      iBase_EntityHandle region, int* sense_out, int* err) {
}
void iGeom_getArrNrmlSense(iGeom_Instance,
      iBase_EntityHandle const* face_handles, int face_handles_size,
      iBase_EntityHandle const* region_handles, int region_handles_size,
      int** sense, int* sense_allocated, int* sense_size, int* err) {
}

/**\brief Get the sense of an edge with respect to a face
 * Get the sense of an edge with respect to a face.  Sense returned is -1, 0, or 1,
 * representing "reversed", "both", or "forward".  "both" sense indicates that edge bounds
 * the face once with each sense.
 * \param edge Edge being queried
 * \param face Face being queried
 * \param sense_out Sense of edge with respect to face
 */

void iGeom_getEgFcSense(iGeom_Instance instance, iBase_EntityHandle edge,
      iBase_EntityHandle face, int* sense_out, int* err) {
   // this one is important, for establishing the orientation of the edges in faces
   // bummer, I "thought" it is already implemented
   // use senses
   GETGTT(instance);
   std::vector<EntityHandle> faces;
   std::vector<int> senses; // 0 is forward and 1 is backward
   ErrorCode rval = _my_geomTopoTool->get_senses( MBH_cast(edge), faces, senses);
   MBERRORR("Failed to get edge senses in iGeom_getEgFcSense.");
   //
   int index = -1;
   EntityHandle mbfaceSet = MBH_cast(face);
   bool sense_forward = false;
   bool sense_reverse = false;
   for (unsigned int i = 0; i<faces.size(); i++)
   {
      if (faces[i] == mbfaceSet)
      {
         index = i;
         if (senses[i]==0)
               sense_forward =true;
            else
               sense_reverse = true;
      }
   }
   if (index == -1)
   {
      *err = iBase_FAILURE;
      return;
   }

   // 0 is not possible for us, but maybe we should consider this?
   if (sense_forward && sense_reverse)
      *sense_out = 0; // is it really possible for a nice geometry ?
   else
   {
      if (sense_forward) // only sense forward
         *sense_out = 1;
      else if (sense_reverse)
         *sense_out = -1;
   }
   RETURN(iBase_SUCCESS);

}
void iGeom_getEgFcArrSense(iGeom_Instance,
      iBase_EntityHandle const* edge_handles, int edge_handles_size,
      iBase_EntityHandle const* face_handles, int face_handles_size,
      int** sense, int* sense_allocated, int* sense_size, int* err) {
}
void iGeom_getEgVtxSense(iGeom_Instance instance, iBase_EntityHandle edge,
      iBase_EntityHandle vertex1, iBase_EntityHandle vertex2, int* sense_out,
      int* err) {
   // need to decide first or second vertex
   // important for moab
   iBase_EntityHandle *verts1 = NULL;
   int verts_alloc1 = 0, verts_size1;
   iMesh_getEntities(IMESH_INSTANCE(instance),
         reinterpret_cast<iBase_EntitySetHandle> (vertex1),
         iBase_ALL_TYPES, iMesh_POINT, &verts1, &verts_alloc1, &verts_size1, err);
   ERRORR("Failed to get vertices.");

   if (verts_size1 != 1) {
      *err = iBase_FAILURE;
      ERRORR("Vertex has multiple points.");
   }
   iBase_EntityHandle *verts2 = NULL;
   int verts_alloc2 = 0, verts_size2;
   iMesh_getEntities(IMESH_INSTANCE(instance),
         reinterpret_cast<iBase_EntitySetHandle> (vertex2),
         iBase_ALL_TYPES, iMesh_POINT, &verts2, &verts_alloc2, &verts_size2, err);
   ERRORR("Failed to get vertices.");

   if (verts_size2 != 1) {
      *err = iBase_FAILURE;
      ERRORR("Vertex has multiple points.");
   }
   // now get the edges, and get the first node and the last node in sequence of edges
   // the order is important...
   iBase_EntityHandle *medges = NULL;
   int medges_alloc = 0, medges_size;
   iMesh_getEntities(IMESH_INSTANCE(instance),
         reinterpret_cast<iBase_EntitySetHandle> (edge),
         iBase_ALL_TYPES, iMesh_LINE_SEGMENT, &medges, &medges_alloc, &medges_size, err);
   ERRORR("Failed to get mesh edges.");
   // get the first node and the last node, then compare with verts1[0] and verts2[0]

   const EntityHandle* conn;
   int len;
   EntityHandle startNode, endNode;
   ErrorCode rval = MBI->get_connectivity(MBH_cast(medges[0]), conn, len);
   MBERRORR("Failed to get edge connectivity in iGeom_getEgVtxSense.");
   startNode = conn[0];
   rval = MBI->get_connectivity(MBH_cast(medges[medges_size-1]), conn, len);
   MBERRORR("Failed to get edge connectivity in iGeom_getEgVtxSense.");
   endNode = conn[1];
   if (startNode == endNode && MBH_cast(verts1[0])==startNode)
   {
      * sense_out = 0; // periodic
   }
   if (startNode == MBH_cast(verts1[0]) && endNode == MBH_cast(verts2[0]) )
   {
      * sense_out = 1; // forward
   }
   if (startNode == MBH_cast(verts2[0]) && endNode == MBH_cast(verts1[0]) )
   {
      * sense_out = -1; // reverse
   }

   free (medges);
   free(verts1);
   free(verts2);

   RETURN(iBase_SUCCESS);
}
void iGeom_getEgVtxArrSense(iGeom_Instance,
      iBase_EntityHandle const* edge_handles, int edge_handles_size,
      iBase_EntityHandle const* vertex_handles_1, int veretx_handles_1_size,
      iBase_EntityHandle const* vertex_handles_2, int vertex_handles_2_size,
      int** sense, int* sense_allocated, int* sense_size, int* err) {
}
void iGeom_measure(iGeom_Instance instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      double** measures, int* measures_allocated, int* measures_size, int* err) {
   CHECK_SIZE(*measures, *measures_allocated, entity_handles_size, double, NULL);
   for (int i = 0; i < entity_handles_size; i++) {
      (*measures)[i] = 0.;

      int type;
      iGeom_getEntType(instance, entity_handles[i], &type, err);
      ERRORR("Failed to get entity type in iGeom_measure.");

      if (type == 1) { // edge
         iBase_EntityHandle *edges = NULL;
         int edges_alloc = 0, edges_size;
         iMesh_getEntities(IMESH_INSTANCE(instance),
               reinterpret_cast<iBase_EntitySetHandle> (entity_handles[i]),
               iBase_ALL_TYPES, iMesh_LINE_SEGMENT, &edges, &edges_alloc,
               &edges_size, err);
         ERRORR("Failed to get edges.");

         iBase_EntityHandle *adj = NULL;
         int adj_alloc = 0, adj_size;
         int* offset = NULL;
         int offset_alloc, offset_size;
         iMesh_getEntArrAdj(IMESH_INSTANCE(instance), edges, edges_size,
               iBase_VERTEX, &adj, &adj_alloc, &adj_size, &offset,
               &offset_alloc, &offset_size, err);
         ERRORR("Failed to get entity adjacencies in iGeom_measure.");

         for (int j = 0; j < edges_size; j++) {
            double p1[3], p2[3];
            iMesh_getVtxCoord(IMESH_INSTANCE(instance), adj[offset[j]], &p1[0],
                  &p1[1], &p1[2], err);
            ERRORR("Failed to get vertex coordinates in iGeom_measure.");
            iMesh_getVtxCoord(IMESH_INSTANCE(instance), adj[offset[j] + 1],
                  &p2[0], &p2[1], &p2[2], err);
            ERRORR("Failed to get vertex coordinates in iGeom_measure.");
            (*measures)[i] += get_edge_length(p1, p2);
         }
      }
      if (type == 2) { // surface
         // get triangles in surface
         iBase_EntityHandle *tris = NULL;
         int tris_alloc = 0, tris_size;
         iMesh_getEntities(IMESH_INSTANCE(instance),
               reinterpret_cast<iBase_EntitySetHandle> (entity_handles[i]),
               iBase_ALL_TYPES, iMesh_TRIANGLE, &tris, &tris_alloc, &tris_size,
               err);
         ERRORR("Failed to get triangles in iGeom_measure.");

         iBase_EntityHandle *adj = NULL;
         int adj_alloc = 0, adj_size;
         int* offset = NULL;
         int offset_alloc, offset_size;
         iMesh_getEntArrAdj(IMESH_INSTANCE(instance), tris, tris_size,
               iBase_VERTEX, &adj, &adj_alloc, &adj_size, &offset,
               &offset_alloc, &offset_size, err);
         ERRORR("Failed to get triangle adjacencies in iGeom_measure.");

         // calculate sum of area of triangles
         double* p;
         CartVect coords[3];
         for (int j = 0; j < tris_size; j++) {
            for (int k = 0; k < 3; k++) {
               p = coords[k].array();
               iMesh_getVtxCoord(IMESH_INSTANCE(instance), adj[offset[j] + k],
                     p, p + 1, p + 2, err);
               ERRORR("Failed to get vertex coordinates in iGeom_measure.");
            }
            coords[1] -= coords[0];
            coords[2] -= coords[0];
            coords[0] = coords[1] * coords[2];
            (*measures)[i] += coords[0].length();
         }
         (*measures)[i] *= .5;
      }
   }

   RETURN(iBase_SUCCESS);
}

void iGeom_getFaceType(iGeom_Instance, iBase_EntityHandle face_handle,
      char* face_type, int* err, int* face_type_length) {
}
void iGeom_getParametric(iGeom_Instance, int* is_parametric, int* err) {
}
void iGeom_isEntParametric(iGeom_Instance, iBase_EntityHandle entity_handle,
      int* parametric, int* err) {
}
void iGeom_isArrParametric(iGeom_Instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int** is_parametric, int* is_parametric_allocated,
      int* is_parametric_size, int* err) {
}
void iGeom_getEntUVtoXYZ(iGeom_Instance, iBase_EntityHandle entity_handle,
      double u, double v, double* x, double* y, double* z, int* err) {
}
void iGeom_getArrUVtoXYZ(iGeom_Instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int storage_order, double const* uv, int uv_size, double** coordinates,
      int* coordinates_allocated, int* coordinates_size, int* err) {
}

void iGeom_getEntUtoXYZ(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, double u, double* x, double* y,
      double* z, int* err) {
   int type, i, j;
   double tot_length = 0., old_length;
   iGeom_getEntType(instance, entity_handle, &type, err);
   ERRORR("Failed to get entity type in iGeom_getEntUtoXYZ.");

   if (type == 1) { // edge
      if (_smooth)
      {
         // first, find the edge
         EntityHandle geoSet = MBH_cast(entity_handle);
         SmoothCurveEval* smthEdge = _edges[geoSet];
         smthEdge ->position_from_u( u, *x, *y, *z );
         RETURN(iBase_SUCCESS);
      }
      else
      {
         // get edges and verticies of this geometry
         iBase_EntityHandle *edges = NULL;
         int edges_alloc = 0, edges_size;
         iMesh_getEntities(IMESH_INSTANCE(instance),
               reinterpret_cast<iBase_EntitySetHandle> (entity_handle),
               iBase_ALL_TYPES, iMesh_LINE_SEGMENT, &edges, &edges_alloc,
               &edges_size, err);
         ERRORR("Failed to get edges in iGeom_getEntUtoXYZ.");

         iBase_EntityHandle *verts = NULL;
         int verts_alloc = 0, verts_size;
         int* offset = NULL;
         int offset_alloc, offset_size;
         iMesh_getEntArrAdj(IMESH_INSTANCE(instance), edges, edges_size,
               iBase_VERTEX, &verts, &verts_alloc, &verts_size, &offset,
               &offset_alloc, &offset_size, err);
         ERRORR("Failed to get entity adjacencies in iGeom_getEntUtoXYZ.");

         // make vertex loop
         std::vector<iBase_EntityHandle> loop_verts;
         std::map<iBase_EntityHandle, int> edge_map;
         std::map<iBase_EntityHandle, int>::iterator iter;
         for (i = 1; i < edges_size; i++) {
            edge_map[edges[i]] = i;
         }
         iBase_EntityHandle start_vertex = verts[offset[0]];
         iBase_EntityHandle end_vertex = verts[offset[0] + 1];
         for (i = 0; i < edges_size; i++) {
            loop_verts.push_back(start_vertex);
            //edge_map.erase(edges[i]);

            iBase_EntityHandle *adj_edges = NULL;
            int adj_edges_alloc = 0, adj_edges_size;
            iMesh_getEntAdj(IMESH_INSTANCE(instance), end_vertex, iBase_EDGE,
                  &adj_edges, &adj_edges_alloc, &adj_edges_size, err);
            ERRORR("Failed to get entity adjacencies in iGeom_getEntUtoXYZ.");

            for (j = 0; j < adj_edges_size; j++) {
               iter = edge_map.find(adj_edges[j]);
               if (iter != edge_map.end()) {
                  if (end_vertex == verts[offset[iter->second]]) {
                     start_vertex = end_vertex;
                     end_vertex = verts[offset[iter->second] + 1];
                     edge_map.erase(iter->first);
                     break;
                  } else if (end_vertex == verts[offset[iter->second] + 1]) {
                     start_vertex = end_vertex;
                     end_vertex = verts[offset[iter->second]];
                     edge_map.erase(iter->first);
                     break;
                  }
               }
            }
         }

         if (debug_igeom) {
            iBase_EntitySetHandle set;
            iMesh_createEntSet(IMESH_INSTANCE(instance), true, &set, err);
            ERRORR("Problem creating geometry entityset.\n");
            iMesh_addEntArrToSet(IMESH_INSTANCE(instance), &loop_verts[0],
                  loop_verts.size(), set, err);
            ERRORR("Failed to add vertex in entity set\n");

            iBase_EntityHandle *ver = NULL;
            int ver_alloc = 0, ver_size;
            iMesh_getEntities(IMESH_INSTANCE(instance), set, iBase_ALL_TYPES,
                  iMesh_POINT, &ver, &ver_alloc, &ver_size, err);
            ERRORR("Failed to get vertex.");
            std::cout << "ver_size=" << ver_size << std::endl;
            iMesh_save(IMESH_INSTANCE(instance), set, "loop.vtk", 0, err, 8, 0);
         }

         // find proper point in vertex loop with u
         int n_verts = loop_verts.size();
         if (n_verts == edges_size) {
            for (i = 0; i < n_verts; i++) {
               double p1[3], p2[3];
               iMesh_getVtxCoord(IMESH_INSTANCE(instance),
                     loop_verts[i % n_verts], &p1[0], &p1[1], &p1[2], err);
               ERRORR("Failed to get vertex coordinates in iGeom_getEntUtoXYZ.");
               iMesh_getVtxCoord(IMESH_INSTANCE(instance), loop_verts[(i + 1)
                     % n_verts], &p2[0], &p2[1], &p2[2], err);
               ERRORR("Failed to get vertex coordinates in iGeom_getEntUtoXYZ.");
               old_length = tot_length;
               tot_length += get_edge_length(p1, p2);
               if (tot_length > u) {
                  double portion = (u - old_length) / (tot_length - old_length);
                  *x = p1[0] + portion * (p2[0] - p1[0]);
                  *y = p1[1] + portion * (p2[1] - p1[1]);
                  *z = p1[2] + portion * (p2[2] - p1[2]);
                  RETURN(iBase_SUCCESS);
               }
            }
         }
         *err = iBase_FAILURE;
         ERRORR("Failed to xyz point for u.");
      }
   } else
      RETURN(iBase_NOT_SUPPORTED);
}

void iGeom_getArrUtoXYZ(iGeom_Instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      double const* u, int u_size, int storage_order, double** on_coords,
      int* on_coords_allocated, int* on_coords_size, int* err) {
}
void iGeom_getEntXYZtoUV(iGeom_Instance, iBase_EntityHandle entity_handle,
      double x, double y, double z, double* u, double* v, int* err) {
}
void iGeom_getEntXYZtoU(iGeom_Instance, iBase_EntityHandle entity_handle,
      double x, double y, double z, double* u, int* err) {
}
void iGeom_getArrXYZtoUV(iGeom_Instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int storage_order, double const* coordinates, int coordinates_size,
      double** uv, int* uv_allocated, int* uv_size, int* err) {
}
void iGeom_getArrXYZtoU(iGeom_Instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int storage_order, double const* coordinates, int coordinates_size,
      double** u, int* u_allocated, int* u_size, int* err) {
}
void iGeom_getEntXYZtoUVHint(iGeom_Instance, iBase_EntityHandle entity_handle,
      double x, double y, double z, double* u, double* v, int* err) {
}
void iGeom_getArrXYZtoUVHint(iGeom_Instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int storage_order, double const* coords, int coords_size, double** uv,
      int* uv_allocated, int* uv_size, int* err) {
}
void iGeom_getEntUVRange(iGeom_Instance, iBase_EntityHandle entity_handle,
      double* u_min, double* v_min, double* u_max, double* v_max, int* err) {
}

void iGeom_getEntURange(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, double* u_min, double* u_max, int* err) {
   int type;
   iGeom_getEntType(instance, entity_handle, &type, err);
   ERRORR("Failed to get entity type in iGeom_getEntURange.");

   if (type == 1) { // edge
      if (_smooth)
      {
         // first, find the edge
         EntityHandle geoSet = MBH_cast(entity_handle);
         SmoothCurveEval* smthEdge = _edges[geoSet];
         smthEdge ->get_param_range(*u_min, *u_max);
         RETURN(iBase_SUCCESS);
      }
      else
      {
         iBase_EntityHandle *edges = NULL;
         int edges_alloc = 0, edges_size;
         iMesh_getEntities(IMESH_INSTANCE(instance),
               reinterpret_cast<iBase_EntitySetHandle> (entity_handle),
               iBase_ALL_TYPES, iMesh_LINE_SEGMENT, &edges, &edges_alloc,
               &edges_size, err);
         ERRORR("Failed to get edges in iGeom_getEntURange.");

         iBase_EntityHandle *adj = NULL;
         int adj_alloc = 0, adj_size;
         int* offset = NULL;
         int offset_alloc, offset_size;
         iMesh_getEntArrAdj(IMESH_INSTANCE(instance), edges, edges_size,
               iBase_VERTEX, &adj, &adj_alloc, &adj_size, &offset, &offset_alloc,
               &offset_size, err);
         ERRORR("Failed to get entity adjacencies in iGeom_getEntURange.");

         *u_min = *u_max = 0.;
         for (int j = 0; j < edges_size; j++) {
            double p1[3], p2[3];
            iMesh_getVtxCoord(IMESH_INSTANCE(instance), adj[offset[j]], &p1[0],
                  &p1[1], &p1[2], err);
            ERRORR("Failed to get vertex coordinates in iGeom_getEntURange.");
            iMesh_getVtxCoord(IMESH_INSTANCE(instance), adj[offset[j] + 1],
                  &p2[0], &p2[1], &p2[2], err);
            ERRORR("Failed to get vertex coordinates in iGeom_getEntURange.");
            *u_max += get_edge_length(p1, p2);
         }
      }
   } else
      RETURN(iBase_NOT_SUPPORTED);

   RETURN(iBase_SUCCESS);
}

void iGeom_getArrUVRange(iGeom_Instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int storage_order, double** uv_min, int* uv_min_allocated,
      int* uv_min_size, double** uv_max, int* uv_max_allocated,
      int* uv_max_size, int* err) {
}
void iGeom_getArrURange(iGeom_Instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      double** u_min, int* u_min_allocated, int* u_min_size, double** u_max,
      int* u_max_allocated, int* u_max_size, int* err) {
}
void iGeom_getEntUtoUV(iGeom_Instance, iBase_EntityHandle edge_handle,
      iBase_EntityHandle face_handle, double in_u, double* u, double* v,
      int* err) {
}
void iGeom_getVtxToUV(iGeom_Instance, iBase_EntityHandle vertex_handle,
      iBase_EntityHandle face_handle, double* u, double* v, int* err) {
}
void iGeom_getVtxToU(iGeom_Instance, iBase_EntityHandle vertex_handle,
      iBase_EntityHandle edge_handle, double* u, int* err) {
}
void iGeom_getArrUtoUV(iGeom_Instance, iBase_EntityHandle const* edge_handles,
      int edge_handles_size, iBase_EntityHandle const* face_handles,
      int face_handles_size, double const* u_in, int u_in_size,
      int storage_order, double** uv, int* uv_allocated, int* uv_size, int* err) {
}
void iGeom_getVtxArrToUV(iGeom_Instance,
      iBase_EntityHandle const* vertex_handles, int vertex_handles_size,
      iBase_EntityHandle const* face_handles, int face_handles_size,
      int storage_order, double** uv, int* uv_allocated, int* uv_size, int* err) {
}
void iGeom_getVtxArrToU(iGeom_Instance,
      iBase_EntityHandle const* vertex_handles, int vertex_handles_size,
      iBase_EntityHandle const* edge_handles, int edge_handles_size,
      double** u, int* u_allocated, int* u_size, int* err) {
}
void iGeom_getEntNrmlUV(iGeom_Instance, iBase_EntityHandle entity_handle,
      double u, double v, double* nrml_i, double* nrml_j, double* nrml_k,
      int* err) {
}
void iGeom_getArrNrmlUV(iGeom_Instance, iBase_EntityHandle const* face_handles,
      int face_handles_size, int storage_order, double const* parameters,
      int parameters_size, double** normals, int* normals_allocated,
      int* normals_size, int* err) {
}
void iGeom_getEntTgntU(iGeom_Instance, iBase_EntityHandle entity_handle,
      double u, double* tgnt_i, double* tgnt_j, double* tgnt_k, int* err) {
}
void iGeom_getArrTgntU(iGeom_Instance, iBase_EntityHandle const* edge_handles,
      int edge_handles_size, int storage_order, double const* parameters,
      int parameters_size, double** tangents, int* tangents_allocated,
      int* tangents_size, int* err) {
}
void iGeom_getEnt1stDrvt(iGeom_Instance, iBase_EntityHandle entity_handle,
      double u, double v, double** drvt_u, int* drvt_u_allocated,
      int* drvt_u_size, double** drvt_v, int* dvrt_v_allocated,
      int* dvrt_v_size, int* err) {
}
void iGeom_getArr1stDrvt(iGeom_Instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int storage_order, double const* uv, int uv_size, double** dvtr_u,
      int* dvrt_u_allocated, int* dvrt_u_size, int** u_offset,
      int* u_offset_allocated, int* u_offset_size, double** dvrt_v,
      int* dvrt_v_allocated, int* dvrt_v_size, int** v_offset,
      int* v_offset_allocated, int* v_offset_size, int* err) {
}
void iGeom_getEnt2ndDrvt(iGeom_Instance, iBase_EntityHandle entity_handle,
      double u, double v, double** drvt_uu, int* drvt_uu_allocated,
      int* drvt_uu_size, double** drvt_vv, int* dvrt_vv_allocated,
      int* dvrt_vv_size, double** drvt_uv, int* dvrt_uv_allocated,
      int* dvrt_uv_size, int* err) {
}
void iGeom_getArr2ndDrvt(iGeom_Instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int storage_order, double const* uv, int uv_size, double** dvtr_uu,
      int* dvrt_uu_allocated, int* dvrt_uu_size, int** uu_offset,
      int* uu_offset_allocated, int* uu_offset_size, double** dvtr_vv,
      int* dvrt_vv_allocated, int* dvrt_vv_size, int** vv_offset,
      int* vv_offset_allocated, int* vv_offset_size, double** dvrt_uv,
      int* dvrt_uv_allocated, int* dvrt_uv_size, int** uv_offset,
      int* uv_offset_allocated, int* uv_offset_size, int* err) {
}
void iGeom_getFcCvtrUV(iGeom_Instance, iBase_EntityHandle entity_handle,
      double u, double v, double* cvtr1_i, double* cvtr1_j, double* cvtr1_k,
      double* cvtr2_i, double* cvtr2_j, double* cvtr2_k, int* err) {
}
void iGeom_getFcArrCvtrUV(iGeom_Instance,
      iBase_EntityHandle const* face_handles, int face_handles_size,
      int storage_order, double const* uv, int uv_size, double** cvtr_1,
      int* cvtr_1_allocated, int* cvtr_1_size, double** cvtr_2,
      int* cvtr_2_allocated, int* cvtr_2_size, int* err) {
}
void iGeom_isEntPeriodic(iGeom_Instance, iBase_EntityHandle entity_handle,
      int* in_u, int* in_v, int* err) {
}
void iGeom_isArrPeriodic(iGeom_Instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int** in_uv, int* in_uv_allocated, int* in_uv_size, int* err) {
}
void iGeom_isFcDegenerate(iGeom_Instance, iBase_EntityHandle face_handle,
      int* is_degenerate, int* err) {
}
void iGeom_isFcArrDegenerate(iGeom_Instance,
      iBase_EntityHandle const* face_handles, int face_handles_size,
      int** degenerate, int* degenerate_allocated, int* degenerate_size,
      int* err) {
}

void iGeom_getArrTolerance(iGeom_Instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      double** tolerances, int* tolerances_allocated, int* tolerances_size,
      int* err) {
}

void iGeom_initEntIter(iGeom_Instance, iBase_EntitySetHandle entity_set_handle,
      int entity_dimension, iBase_EntityIterator* entity_iterator, int* err) {
}

void iGeom_initEntArrIter(iGeom_Instance,
      iBase_EntitySetHandle entity_set_handle, int entity_dimension,
      int requested_array_size, iBase_EntityArrIterator* entArr_iterator,
      int* err) {
}

void iGeom_getNextEntIter(iGeom_Instance, iBase_EntityIterator,
      iBase_EntityHandle* entity_handle, int* has_data, int* err) {
}

void iGeom_getNextEntArrIter(iGeom_Instance, iBase_EntityArrIterator,
      iBase_EntityHandle** entity_handles, int* entity_handles_allocated,
      int* entity_handles_size, int* has_data, int* err) {
}

void iGeom_resetEntIter(iGeom_Instance, iBase_EntityIterator, int* err) {
}

void iGeom_resetEntArrIter(iGeom_Instance, iBase_EntityArrIterator, int* err) {
}

void iGeom_endEntIter(iGeom_Instance, iBase_EntityIterator, int* err) {
}

void iGeom_endEntArrIter(iGeom_Instance, iBase_EntityArrIterator, int* err) {
}

void iGeom_copyEnt(iGeom_Instance, iBase_EntityHandle source,
      iBase_EntityHandle* copy, int* err) {
}

void iGeom_sweepEntAboutAxis(iGeom_Instance, iBase_EntityHandle geom_entity,
      double angle, double axis_normal_x, double axis_normal_y,
      double axis_normal_z, iBase_EntityHandle* geom_entity2, int* err) {
}

void iGeom_deleteAll(iGeom_Instance instance, int* err) {
   ErrorCode rval;
   for (int i = 0; i < 4; i++) {
      rval = MBI->delete_entities(_my_gsets[i]);
      MBERRORR("Failed to delete entities in iGeom_deleteAll.");
      _my_gsets[i].clear();
   }

   RETURN(iBase_SUCCESS);
}

void iGeom_deleteEnt(iGeom_Instance instance, iBase_EntityHandle entity_handle,
      int* err) {
   int type;
   iGeom_getEntType(instance, entity_handle, &type, err);
   ERRORR("Failed to get entity type in iGeom_deleteEnt.");

   Range::iterator iter = _my_gsets[type].find(MBH_cast(entity_handle));
   if (iter == _my_gsets[type].end()) {
      RETURN(iBase_INVALID_ENTITY_HANDLE);
   }
   _my_gsets[type].erase(iter);

   EntityHandle this_entity = MBH_cast(entity_handle);
   ErrorCode rval = MBI->delete_entities(&this_entity, 1);
   MBERRORR("Failed to delete entity.");
}

void iGeom_createSphere(iGeom_Instance, double radius,
      iBase_EntityHandle* sphere_handle_out, int* err) {
}

void iGeom_createPrism(iGeom_Instance, double height, int n_sides,
      double major_rad, double minor_rad, iBase_EntityHandle* prism_handle_out,
      int* err) {
}

void iGeom_createBrick(iGeom_Instance, double x, double y, double z,
      iBase_EntityHandle* geom_entity, int* err) {
}

void iGeom_createCylinder(iGeom_Instance, double height, double major_rad,
      double minor_rad, iBase_EntityHandle* geom_entity, int* err) {
}

void iGeom_createCone(iGeom_Instance, double height, double major_rad_base,
      double minor_rad_base, double rad_top, iBase_EntityHandle* geom_entity,
      int* err) {
}

void iGeom_createTorus(iGeom_Instance, double major_rad, double minor_rad,
      iBase_EntityHandle* geom_entity, int* err) {
}

void iGeom_moveEnt(iGeom_Instance, iBase_EntityHandle geom_entity, double x,
      double y, double z, int* err) {
}

void iGeom_rotateEnt(iGeom_Instance, iBase_EntityHandle geom_entity,
      double angle, double axis_normal_x, double axis_normal_y,
      double axis_normal_z, int* err) {
}

void iGeom_reflectEnt(iGeom_Instance, iBase_EntityHandle geom_entity,
      double plane_normal_x, double plane_normal_y, double plane_normal_z,
      int* err) {
}

void iGeom_scaleEnt(iGeom_Instance, iBase_EntityHandle geom_entity,
      double scale_x, double scale_y, double scale_z, int* err) {
}

void iGeom_uniteEnts(iGeom_Instance, iBase_EntityHandle const* geom_entities,
      int geom_entities_size, iBase_EntityHandle* geom_entity, int* err) {
}

void iGeom_subtractEnts(iGeom_Instance, iBase_EntityHandle blank,
      iBase_EntityHandle tool, iBase_EntityHandle* geom_entity, int* err) {
}

void iGeom_intersectEnts(iGeom_Instance, iBase_EntityHandle entity2,
      iBase_EntityHandle entity1, iBase_EntityHandle* geom_entity, int* err) {
}

void iGeom_sectionEnt(iGeom_Instance, iBase_EntityHandle geom_entity,
      double plane_normal_x, double plane_normal_y, double plane_normal_z,
      double offset, int reverse, iBase_EntityHandle* geom_entity2, int* err) {
}

void iGeom_imprintEnts(iGeom_Instance, iBase_EntityHandle const* geom_entities,
      int geom_entities_size, int* err) {
}

void iGeom_mergeEnts(iGeom_Instance, iBase_EntityHandle const* geom_entities,
      int geom_entities_size, double tolerance, int* err) {
}

void iGeom_createEntSet(iGeom_Instance instance, int isList,
      iBase_EntitySetHandle* entity_set_created, int *err) {
   iMesh_createEntSet(IMESH_INSTANCE(instance), isList, entity_set_created, err);
   ERRORR("Failed to get create entity set.");
}

void iGeom_destroyEntSet(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, int *err) {
}

void iGeom_isList(iGeom_Instance instance, iBase_EntitySetHandle entity_set,
      int *is_list, int *err) {
}

void iGeom_getNumEntSets(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set_handle, int num_hops, int *num_sets,
      int *err) {
   iMesh_getNumEntSets(IMESH_INSTANCE(instance), entity_set_handle, num_hops,
         num_sets, err);
   ERRORR("Failed to get number of entity sets.");
}

void iGeom_getEntSets(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set_handle, int num_hops,
      iBase_EntitySetHandle** contained_set_handles,
      int* contained_set_handles_allocated, int* contained_set_handles_size,
      int *err) {
   iMesh_getEntSets(IMESH_INSTANCE(instance), entity_set_handle, num_hops,
         contained_set_handles, contained_set_handles_allocated,
         contained_set_handles_size, err);
   ERRORR("Failed to get entity sets.");
}

void iGeom_addEntToSet(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, iBase_EntitySetHandle entity_set,
      int *err) {
   iMesh_addEntToSet(IMESH_INSTANCE(instance), entity_handle, entity_set, err);
   ERRORR("Failed to add entity to set.");
}

void iGeom_rmvEntFromSet(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, iBase_EntitySetHandle entity_set,
      int *err) {
   iMesh_rmvEntFromSet(IMESH_INSTANCE(instance), entity_handle, entity_set, err);
   ERRORR("Failed to remove entity from set.");
}

void iGeom_addEntArrToSet(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_EntitySetHandle entity_set, int *err) {
   iMesh_addEntArrToSet(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, entity_set, err);
   ERRORR("Failed to add entities to set.");
}

void iGeom_rmvEntArrFromSet(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_EntitySetHandle entity_set, int *err) {
   iMesh_rmvEntArrFromSet(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, entity_set, err);
   ERRORR("Failed to remove entities from set.");
}

void iGeom_addEntSet(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set_to_add,
      iBase_EntitySetHandle entity_set_handle, int *err) {
   iMesh_addEntSet(IMESH_INSTANCE(instance), entity_set_to_add,
         entity_set_handle, err);
   ERRORR("Failed to add entity set to entity set.");
}

void iGeom_rmvEntSet(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set_to_remove,
      iBase_EntitySetHandle entity_set_handle, int *err) {
   iMesh_rmvEntSet(IMESH_INSTANCE(instance), entity_set_to_remove,
         entity_set_handle, err);
   ERRORR("Failed to remove entity set from entity set.");
}

void iGeom_isEntContained(iGeom_Instance instance,
      iBase_EntitySetHandle containing_entity_set,
      iBase_EntityHandle contained_entity, int *is_contained, int *err) {
   iMesh_isEntContained(IMESH_INSTANCE(instance), containing_entity_set,
         contained_entity, is_contained, err);
   ERRORR("Failed to check if entity is contained to entity set.");
}

void iGeom_isEntArrContained(iGeom_Instance instance,
      iBase_EntitySetHandle containing_set,
      const iBase_EntityHandle* entity_handles, int num_entity_handles,
      int** is_contained, int* is_contained_allocated, int* is_contained_size,
      int* err) {
   iMesh_isEntArrContained(IMESH_INSTANCE(instance), containing_set,
         entity_handles, num_entity_handles, is_contained,
         is_contained_allocated, is_contained_size, err);
   ERRORR("Failed to check if entities are contained to entity sets.");
}

void iGeom_isEntSetContained(iGeom_Instance instance,
      iBase_EntitySetHandle containing_entity_set,
      iBase_EntitySetHandle contained_entity_set, int *is_contained, int *err) {
   iMesh_isEntSetContained(IMESH_INSTANCE(instance), containing_entity_set,
         contained_entity_set, is_contained, err);
   ERRORR("Failed to check if entity set is contained to entity set.");
}

void iGeom_addPrntChld(iGeom_Instance instance,
      iBase_EntitySetHandle parent_entity_set,
      iBase_EntitySetHandle child_entity_set, int *err) {
   iMesh_addPrntChld(IMESH_INSTANCE(instance), parent_entity_set,
         child_entity_set, err);
   ERRORR("Failed to add parent and child relation of geometry sets.");
}

void iGeom_rmvPrntChld(iGeom_Instance instance,
      iBase_EntitySetHandle parent_entity_set,
      iBase_EntitySetHandle child_entity_set, int *err) {
   iMesh_rmvPrntChld(IMESH_INSTANCE(instance), parent_entity_set,
         child_entity_set, err);
   ERRORR("Failed to remove parent and child relation of geometry sets.");
}

void iGeom_isChildOf(iGeom_Instance instance,
      iBase_EntitySetHandle parent_entity_set,
      iBase_EntitySetHandle child_entity_set, int *is_child, int *err) {
   iMesh_isChildOf(IMESH_INSTANCE(instance), parent_entity_set,
         child_entity_set, is_child, err);
   ERRORR("Failed to check if there is a parent/child relation.");
}

void iGeom_getNumChld(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, int num_hops, int *num_child, int *err) {
   iMesh_getNumChld(IMESH_INSTANCE(instance), entity_set, num_hops, num_child,
         err);
   ERRORR("Failed to get number of children.");
}

void iGeom_getNumPrnt(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, int num_hops, int *num_parent, int *err) {
   iMesh_getNumPrnt(IMESH_INSTANCE(instance), entity_set, num_hops, num_parent,
         err);
   ERRORR("Failed to get number of parents.");
}

void iGeom_getChldn(iGeom_Instance instance,
      iBase_EntitySetHandle from_entity_set, int num_hops,
      iBase_EntitySetHandle** entity_set_handles,
      int* entity_set_handles_allocated, int* entity_set_handles_size, int *err) {
   iMesh_getChldn(IMESH_INSTANCE(instance), from_entity_set, num_hops,
         entity_set_handles, entity_set_handles_allocated,
         entity_set_handles_size, err);
   ERRORR("Failed to get children.");
}

void iGeom_getPrnts(iGeom_Instance instance,
      iBase_EntitySetHandle from_entity_set, int num_hops,
      iBase_EntitySetHandle** entity_set_handles,
      int* entity_set_handles_allocated, int* entity_set_handles_size, int *err) {
   iMesh_getPrnts(IMESH_INSTANCE(instance), from_entity_set, num_hops,
         entity_set_handles, entity_set_handles_allocated,
         entity_set_handles_size, err);
   ERRORR("Failed to get parents.");
}

void iGeom_createTag(iGeom_Instance instance, const char* tag_name,
      int tag_size, int tag_type, iBase_TagHandle* tag_handle, int *err,
      int tag_name_len) {

   iMesh_createTag(IMESH_INSTANCE(instance), tag_name, tag_size, tag_type,
         tag_handle, err, tag_name_len);
   ERRORR("Failure to create tag.");
}

void iGeom_destroyTag(iGeom_Instance instance, iBase_TagHandle tag_handle,
      int forced, int *err) {
   iMesh_destroyTag(IMESH_INSTANCE(instance), tag_handle, forced, err);
   ERRORR("Failure to destroy tag.");
}

void iGeom_getTagName(iGeom_Instance instance, iBase_TagHandle tag_handle,
      char *name, int* err, int name_len) {
   iMesh_getTagName(IMESH_INSTANCE(instance), tag_handle, name, err, name_len);
   ERRORR("Failure to get tag name.");
}

void iGeom_getTagSizeValues(iGeom_Instance instance,
      iBase_TagHandle tag_handle, int *tag_size, int *err) {
   iMesh_getTagSizeValues(IMESH_INSTANCE(instance), tag_handle, tag_size, err);
   ERRORR("Failure to get numbers tag size.");
}

void iGeom_getTagSizeBytes(iGeom_Instance instance, iBase_TagHandle tag_handle,
      int *tag_size, int *err) {
   iMesh_getTagSizeBytes(IMESH_INSTANCE(instance), tag_handle, tag_size, err);
   ERRORR("Failure to get byte tag size.");
}

void iGeom_getTagHandle(iGeom_Instance instance, const char* tag_name,
      iBase_TagHandle *tag_handle, int *err, int tag_name_len) {
   iMesh_getTagHandle(IMESH_INSTANCE(instance), tag_name, tag_handle, err,
         tag_name_len);
   ERRORR("Failure to get tag name.");
}

void iGeom_getTagType(iGeom_Instance instance, iBase_TagHandle tag_handle,
      int *tag_type, int *err) {
   iMesh_getTagType(IMESH_INSTANCE(instance), tag_handle, tag_type, err);
   ERRORR("Failure to get tag type.");
}

void iGeom_setEntSetData(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set_handle, iBase_TagHandle tag_handle,
      const char* tag_value, int tag_value_size, int *err) {
   iMesh_setEntSetData(IMESH_INSTANCE(instance), entity_set_handle, tag_handle,
         tag_value, tag_value_size, err);
   ERRORR("Failure to set arbitrary data to entity set.");
}

void iGeom_setEntSetIntData(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, iBase_TagHandle tag_handle,
      int tag_value, int *err) {
   iMesh_setEntSetIntData(IMESH_INSTANCE(instance), entity_set, tag_handle,
         tag_value, err);
   ERRORR("Failure to set integer data to entity set.");
}

void iGeom_setEntSetDblData(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, iBase_TagHandle tag_handle,
      double tag_value, int *err) {
   iMesh_setEntSetDblData(IMESH_INSTANCE(instance), entity_set, tag_handle,
         tag_value, err);
   ERRORR("Failure to set double data to entity set.");
}

void iGeom_setEntSetEHData(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, iBase_TagHandle tag_handle,
      iBase_EntityHandle tag_value, int *err) {
   iMesh_setEntSetEHData(IMESH_INSTANCE(instance), entity_set, tag_handle,
         tag_value, err);
   ERRORR("Failure to set entity handle data to entity set.");
}

void iGeom_getEntSetData(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set_handle, iBase_TagHandle tag_handle,
      char** tag_value, int* tag_value_allocated, int* tag_value_size, int *err) {
   iMesh_getEntSetData(IMESH_INSTANCE(instance), entity_set_handle, tag_handle,
         tag_value, tag_value_allocated, tag_value_size, err);
   ERRORR("Failure to get arbitrary data from entity set.");
}

void iGeom_getEntSetIntData(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, iBase_TagHandle tag_handle,
      int *out_data, int *err) {
   iMesh_getEntSetIntData(IMESH_INSTANCE(instance), entity_set, tag_handle,
         out_data, err);
   ERRORR("Failure to get integer data from entity set.");
}

void iGeom_getEntSetDblData(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, iBase_TagHandle tag_handle,
      double *out_data, int *err) {
   iMesh_getEntSetDblData(IMESH_INSTANCE(instance), entity_set, tag_handle,
         out_data, err);
   ERRORR("Failure to get double data from entity set.");
}

void iGeom_getEntSetEHData(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, iBase_TagHandle tag_handle,
      iBase_EntityHandle *out_data, int *err) {
   iMesh_getEntSetEHData(IMESH_INSTANCE(instance), entity_set, tag_handle,
         out_data, err);
   ERRORR("Failure to get double data from entity set.");
}

void iGeom_getAllEntSetTags(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set_handle, iBase_TagHandle** tag_handles,
      int* tag_handles_allocated, int* tag_handles_size, int *err) {
   iMesh_getAllEntSetTags(IMESH_INSTANCE(instance), entity_set_handle,
         tag_handles, tag_handles_allocated, tag_handles_size, err);
   ERRORR("Failure to get double data from entity set.");
}

void iGeom_rmvEntSetTag(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set_handle, iBase_TagHandle tag_handle,
      int *err) {
   iMesh_rmvEntSetTag(IMESH_INSTANCE(instance), entity_set_handle, tag_handle,
         err);
   ERRORR("Failure to remove entity set tag.");
}

void iGeom_getArrData(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, char** tag_values, int* tag_values_allocated,
      int* tag_values_size, int *err) {
   iMesh_getArrData(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, tag_values, tag_values_allocated,
         tag_values_size, err);
   ERRORR("Failure to get tag values of arbitrary type on an array of entities.");
}

void iGeom_getIntArrData(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, int** tag_values, int* tag_values_allocated,
      int* tag_values_size, int *err) {
   iMesh_getIntArrData(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, tag_values, tag_values_allocated,
         tag_values_size, err);
   ERRORR("Failure to get integer tag values on an array of entities.");
}

void iGeom_getDblArrData(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, double** tag_values,
      int* tag_values_allocated, int* tag_values_size, int *err) {
   iMesh_getDblArrData(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, tag_values, tag_values_allocated,
         tag_values_size, err);
   ERRORR("Failure to get double tag values on an array of entities.");
}

void iGeom_getEHArrData(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, iBase_EntityHandle** tag_value,
      int* tag_value_allocated, int* tag_value_size, int *err) {
   iMesh_getEHArrData(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, tag_value, tag_value_allocated,
         tag_value_size, err);
   ERRORR("Failure to get entity handle tag values on an array of entities.");
}

void iGeom_setArrData(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, const char* tag_values, int tag_values_size,
      int *err) {
   iMesh_setArrData(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, tag_values, tag_values_size, err);
   ERRORR("Failure to set tag values of arbitrary type on an array of entities.");
}

void iGeom_setIntArrData(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, const int* tag_values, int tag_values_size,
      int *err) {
   iMesh_setIntArrData(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, tag_values, tag_values_size, err);
   ERRORR("Failure to set interger tag values on an array of entities.");
}

void iGeom_setDblArrData(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, const double* tag_values,
      const int tag_values_size, int *err) {
   iMesh_setDblArrData(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, tag_values, tag_values_size, err);
   ERRORR("Failure to set double tag values on an array of entities.");
}

void iGeom_setEHArrData(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, const iBase_EntityHandle* tag_values,
      int tag_values_size, int *err) {
   iMesh_setEHArrData(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, tag_values, tag_values_size, err);
   ERRORR("Failure to set entity handle tag values on an array of entities.");
}

void iGeom_rmvArrTag(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, int *err) {
   iMesh_rmvArrTag(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, err);
   ERRORR("Failure to remove tag values on an array of entities.");
}

void iGeom_getData(iGeom_Instance instance, iBase_EntityHandle entity_handle,
      iBase_TagHandle tag_handle, char** tag_value, int *tag_value_allocated,
      int *tag_value_size, int *err) {
   iMesh_getData(IMESH_INSTANCE(instance), entity_handle, tag_handle,
         tag_value, tag_value_allocated, tag_value_size, err);
   ERRORR("Failure to get tag values of an entity.");
}

void iGeom_getIntData(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, iBase_TagHandle tag_handle,
      int *out_data, int *err) {
   iMesh_getIntData(IMESH_INSTANCE(instance), entity_handle, tag_handle,
         out_data, err);
   ERRORR("Failure to get integer tag values of an entity.");
}

void iGeom_getDblData(iGeom_Instance instance,
      const iBase_EntityHandle entity_handle, const iBase_TagHandle tag_handle,
      double *out_data, int *err) {
   iMesh_getDblData(IMESH_INSTANCE(instance), entity_handle, tag_handle,
         out_data, err);
   ERRORR("Failure to get double tag values of an entity.");
}

void iGeom_getEHData(iGeom_Instance instance, iBase_EntityHandle entity_handle,
      iBase_TagHandle tag_handle, iBase_EntityHandle *out_data, int *err) {
   iMesh_getEHData(IMESH_INSTANCE(instance), entity_handle, tag_handle,
         out_data, err);
   ERRORR("Failure to get entity handle tag values of an entity.");
}

void iGeom_setData(iGeom_Instance instance, iBase_EntityHandle entity_handle,
      iBase_TagHandle tag_handle, const char* tag_value, int tag_value_size,
      int *err) {
   iMesh_setData(IMESH_INSTANCE(instance), entity_handle, tag_handle,
         tag_value, tag_value_size, err);
   ERRORR("Failure to set tag values of an entity.");
}

void iGeom_setIntData(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, iBase_TagHandle tag_handle,
      int tag_value, int *err) {
   iMesh_setIntData(IMESH_INSTANCE(instance), entity_handle, tag_handle,
         tag_value, err);
   ERRORR("Failure to set integer tag values of an entity.");
}

void iGeom_setDblData(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, iBase_TagHandle tag_handle,
      double tag_value, int *err) {
   iMesh_setDblData(IMESH_INSTANCE(instance), entity_handle, tag_handle,
         tag_value, err);
   ERRORR("Failure to set double tag values of an entity.");
}

void iGeom_setEHData(iGeom_Instance instance, iBase_EntityHandle entity_handle,
      iBase_TagHandle tag_handle, iBase_EntityHandle tag_value, int *err) {
   iMesh_setEHData(IMESH_INSTANCE(instance), entity_handle, tag_handle,
         tag_value, err);
   ERRORR("Failure to set entity handle tag values of an entity.");
}

void iGeom_getAllTags(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, iBase_TagHandle** tag_handles,
      int* tag_handles_allocated, int* tag_handles_size, int *err) {
   iMesh_getAllTags(IMESH_INSTANCE(instance), entity_handle, tag_handles,
         tag_handles_allocated, tag_handles_size, err);
   ERRORR("Failure to get all tags.");
}

void iGeom_rmvTag(iGeom_Instance instance, iBase_EntityHandle entity_handle,
      iBase_TagHandle tag_handle, int *err) {
   iMesh_rmvTag(IMESH_INSTANCE(instance), entity_handle, tag_handle, err);
   ERRORR("Failure to remove tag.");
}

void iGeom_subtract(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set_1, iBase_EntitySetHandle entity_set_2,
      iBase_EntitySetHandle* result_entity_set, int *err) {
}

void iGeom_intersect(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set_1, iBase_EntitySetHandle entity_set_2,
      iBase_EntitySetHandle* result_entity_set, int *err) {
}

void iGeom_unite(iGeom_Instance instance, iBase_EntitySetHandle entity_set_1,
      iBase_EntitySetHandle entity_set_2,
      iBase_EntitySetHandle* result_entity_set, int *err) {
}

static inline void iGeom_processError(iBase_ErrorType code, const char* desc) {
   std::strncpy(iGeom_LAST_ERROR.description, desc,
         sizeof(iGeom_LAST_ERROR.description));
   iGeom_LAST_ERROR.error_type = code;
}

static void iGeom_get_adjacent_entities(iGeom_Instance instance,
      const EntityHandle from, const int to_dim, Range &adjs,
      //std::vector<EntityHandle>& adjs,
      int* err) {
   int this_dim = -1;
   for (int i = 0; i < 4; i++) {
      if (_my_gsets[i].find(from) != _my_gsets[i].end()) {
         this_dim = i;
         break;
      }
   }

   // check target dimension
   if (-1 == this_dim) {
      iGeom_processError(iBase_FAILURE, "Entity not a geometry entity.");
      RETURN(iBase_FAILURE);
   } else if (0 > to_dim || 3 < to_dim) {
      iGeom_processError(iBase_FAILURE, "To dimension must be between 0 and 3.");
      RETURN(iBase_FAILURE);
   } else if (to_dim == this_dim) {
      iGeom_processError(iBase_FAILURE,
            "To dimension must be different from entity dimension.");
      RETURN(iBase_FAILURE);
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

   RETURN(iBase_SUCCESS);
}

double get_edge_length(double* p1, double* p2) {
   return std::sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1])
         * (p1[1] - p2[1]) + (p1[2] - p2[2]) * (p1[2] - p2[2]));
}

