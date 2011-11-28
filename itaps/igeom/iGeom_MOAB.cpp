#include <iostream>
#include <map>
#include "iGeom_MOAB.hpp"
#include "moab/GeomTopoTool.hpp"
#include "moab/OrientedBoxTreeTool.hpp"
#include "moab/CartVect.hpp"
#include "FileOptions.hpp"
#include <stdlib.h>
#include <cstring>
#include <map>
#include "assert.h"

using namespace moab;

static int compare_no_case1(const char *str1, const char *str2, size_t n) {
   for (size_t i = 1; i != n && *str1 && toupper(*str1) == toupper(*str2);
        ++i, ++str1, ++str2);
   return toupper(*str2) - toupper(*str1);
}
// Filter out non-MOAB options and remove the "moab:" prefix
static std::string filter_options1(const char *begin, const char *end)
{
  const char *opt_begin = begin;
  const char *opt_end   = begin;

  std::string filtered;
  bool first = true;

  while (opt_end != end) {
    opt_end = std::find(opt_begin, end, ' ');

    if (opt_end-opt_begin >= 5 && compare_no_case1(opt_begin, "moab:", 5) == 0) {
      if (!first)
        filtered.push_back(';');
      first = false;
      filtered.append(opt_begin+5, opt_end);
    }

    opt_begin = opt_end+1;
  }
  return filtered;
}

bool debug_igeom = false;
bool Debug_surf_eval = false;

#define COPY_RANGE(r, vec) {                      \
    EntityHandle *tmp_ptr = reinterpret_cast<EntityHandle*>(vec);	\
    std::copy(r.begin(), r.end(), tmp_ptr);}

#define TAG_HANDLE(tagh) reinterpret_cast<Tag>(tagh)

#define COPY_DOUBLEVEC(r, vec) {                      \
  double *tmp_ptr = reinterpret_cast<double*>(vec); \
  std::copy(r.begin(), r.end(), tmp_ptr);}

void iGeom_getDescription(iGeom_Instance instance, char* descr, int descr_len) {
  iMesh_getDescription( IMESH_INSTANCE(instance), descr, descr_len);
}

void iGeom_getErrorType(iGeom_Instance instance, /*out*/int *error_type) {
  iMesh_getErrorType( IMESH_INSTANCE(instance), /*out*/error_type) ;
}

void iGeom_newGeom(char const* options, iGeom_Instance* instance_out, int* err,
      int options_len) {

  std::string tmp_options = filter_options1(options, options+options_len);
  FileOptions opts(tmp_options.c_str());
  // process some options?

  MBiGeom **mbigeom = reinterpret_cast<MBiGeom**> (instance_out);
  *mbigeom = NULL;
  *mbigeom = new MBiGeom();
  *err = iBase_SUCCESS;
}

void iGeom_dtor(iGeom_Instance instance, int* err) {
  delete FBE_cast(instance);
  *err = iBase_SUCCESS;
}

void iGeom_load(iGeom_Instance instance, char const* name, char const* options,
    int* err, int name_len, int options_len) {
  // first remove option for smooth facetting

  const char smth[] = "SMOOTH;";
  bool smooth = false;
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
    smooth = true;

  } else {
    reducedOptions = const_cast<char *> (options);
  }
  // load mesh-based geometry
  const EntityHandle* file_set = 0;
  ErrorCode rval = MBI->load_file(name, file_set, reducedOptions);
  CHKERR(rval, "can't load mesh file");

  FBEngine * fbe = FBE_cast(instance);
  if (fbe == NULL) {
    *err = iBase_FAILURE;
    return;
  }
  GeomTopoTool * gtt = GETGTT(instance);
  if (gtt == NULL) {
    *err = iBase_FAILURE;
    return;
  }
  // keep mesh-based geometries in Range
  rval = gtt->find_geomsets();
  CHKERR(rval, "Failure to find geometry lists.");

  if (smooth) fbe->set_smooth();// assumes that initialization did not happen yet

  fbe->Init();// major computation

  RETURN(iBase_SUCCESS);
}

void iGeom_save(iGeom_Instance instance, char const* name, char const* options,
      int* err, int name_len, int options_len) {
   iMesh_save(IMESH_INSTANCE(instance), NULL, name, options, err, name_len,
         options_len);
}

void iGeom_getRootSet(iGeom_Instance instance, iBase_EntitySetHandle* root_set,
      int* err) {
  EntityHandle modelSet;
  ErrorCode rval = FBE_cast(instance)->getRootSet(&modelSet);
  CHKERR(rval,"can't get root set ");
  *root_set =  (iBase_EntitySetHandle)modelSet;
  RETURN(iBase_SUCCESS);
}

void iGeom_getBoundBox(iGeom_Instance instance, double* min_x, double* min_y,
      double* min_z, double* max_x, double* max_y, double* max_z, int* err) {
   RETURN(iBase_NOT_SUPPORTED);
}

void iGeom_getEntities(iGeom_Instance instance,
      iBase_EntitySetHandle set_handle, int entity_type,
      iBase_EntityHandle** entity_handles, int* entity_handles_allocated,
      int* entity_handles_size, int* err) {

   if (0 > entity_type || 4 < entity_type) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "Bad entity type.");
   } else/* 0<= entity_type <= 4) */ {
     Range gentities;
     ErrorCode rval= FBE_cast(instance)->getEntities((EntityHandle)set_handle, entity_type, gentities);
     CHKERR(rval,"can't get entities ");
     *entity_handles_size = gentities.size();

     CHECK_SIZE(*entity_handles, *entity_handles_allocated,
          *entity_handles_size, iBase_EntityHandle, NULL);
     COPY_RANGE(gentities, *entity_handles);
   }

   RETURN(iBase_SUCCESS);
}

void iGeom_getNumOfType(iGeom_Instance instance,
      iBase_EntitySetHandle set_handle, int entity_type, int* num_out, int* err) {
   if (0 > entity_type || 3 < entity_type) {
      ERROR(iBase_INVALID_ENTITY_TYPE, "Bad entity type.");
   }
   ErrorCode rval = FBE_cast(instance)->getNumOfType((EntityHandle)set_handle, entity_type, num_out);
   CHKERR(rval,"can't get number of type ");

   RETURN(iBase_SUCCESS);
}


void iGeom_getEntType(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, int* type, int* err) {

  ErrorCode rval = FBE_cast(instance)->getEntType((EntityHandle)entity_handle, type);
  CHKERR(rval,"can't get entity type ");

  RETURN(iBase_SUCCESS);
}

void iGeom_getArrType(iGeom_Instance instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int** type, int* type_allocated, int* type_size, int* err) {
   CHECK_SIZE(*type, *type_allocated, *type_size, int, NULL);

   int tmp_err;

   for (int i = 0; i < entity_handles_size; i++) {
      iGeom_getEntType(instance, entity_handles[i], *type + i, &tmp_err);
      if (iBase_SUCCESS != tmp_err) {
         ERROR(tmp_err, "Failed to get entity type in iGeom_getArrType.");
      }
   }

   RETURN(iBase_SUCCESS);
}

void iGeom_getEntAdj(iGeom_Instance instance, iBase_EntityHandle entity_handle,
      int to_dimension, iBase_EntityHandle** adj_entities,
      int* adj_entities_allocated, int* adj_entities_size, int* err) {
   Range adjs;
   EntityHandle this_ent = MBH_cast(entity_handle);

   ErrorCode rval = FBE_cast(instance)->getEntAdj(this_ent, to_dimension,
       adjs);

   CHKERR(rval, "Failed to get adjacent entities in iGeom_getEntAdj.");

   // copy adjacent entities
   *adj_entities_size = adjs.size();
   CHECK_SIZE(*adj_entities, *adj_entities_allocated,
         *adj_entities_size, iBase_EntityHandle, NULL);
   COPY_RANGE(adjs, *adj_entities);

   RETURN(iBase_SUCCESS);
}

// I suspect this is wrong
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
      ErrorCode rval = FBE_cast(instance)->getEntAdj( MBH_cast(entity_handles[i]),
          requested_entity_type,
          temp_range);
      CHKERR(rval, "Failed to get adjacent entities in iGeom_getArrAdj.");
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
   ErrorCode rval = FBE_cast(instance)->getEntAdj(MBH_cast(entity_handle), bridge_dimension,
       bridge_ents);

   CHKERR(rval, "Failed to get adjacent entities in iGeom_getEnt2ndAdj.");

   Range::iterator iter, jter, kter, end_jter;
   Range::iterator end_iter = bridge_ents.end();
   for (iter = bridge_ents.begin(); iter != end_iter; iter++) {
     rval = FBE_cast(instance)->getEntAdj(*iter, to_dimension,
         tmp_ents);

      CHKERR(rval, "Failed to get adjacent entities in iGeom_getEnt2ndAdj.");

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
   // not implemented
// who would need this monster, anyway?
   RETURN(iBase_FAILURE);
}

void iGeom_isEntAdj(iGeom_Instance instance, iBase_EntityHandle entity_handle1,
      iBase_EntityHandle entity_handle2, int* are_adjacent, int* err) {

  bool adjacent_out;
  ErrorCode rval = FBE_cast(instance)->isEntAdj(MBH_cast(entity_handle1), MBH_cast(entity_handle2),
     adjacent_out);
  CHKERR(rval, "Failed to get adjacent info");
  *are_adjacent = (int)adjacent_out; // 0 or 1, really

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
      FWDERR();

      index1 += index1_step;
      index2 += index2_step;
   }

   RETURN(iBase_SUCCESS);
}

void iGeom_getEntClosestPt(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, double near_x, double near_y,
      double near_z, double* on_x, double* on_y, double* on_z, int* err) {

   ErrorCode rval = FBE_cast(instance)->getEntClosestPt(MBH_cast(entity_handle), near_x,
        near_y, near_z, on_x, on_y, on_z);
   CHKERR(rval, "Failed to get closest point");

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
      FWDERR();
   }

   RETURN(iBase_SUCCESS);
}

void iGeom_getEntNrmlXYZ(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, double x, double y, double z,
      double* nrml_i, double* nrml_j, double* nrml_k, int* err) {

  ErrorCode rval = FBE_cast(instance)->getEntNrmlXYZ(MBH_cast(entity_handle),  x,
       y,  z,   nrml_i,  nrml_j,  nrml_k);
  CHKERR(rval, "Failed to get normal");
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
      ERROR(iBase_INVALID_ENTITY_COUNT, "Mismatched array sizes");
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
      FWDERR();

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
   FWDERR();

   if (type != 2 && type != 3) {
      ERROR(iBase_INVALID_ENTITY_TYPE,
             "Entities passed into gentityNormal must be face or volume.");
   }

   // do 2 searches, so it is not fast enough
   iGeom_getEntClosestPt(instance,
          entity_handle, x, y, z,  pt_x, pt_y, pt_z,  err);

   FWDERR();
   iGeom_getEntNrmlXYZ(instance,
         entity_handle, *pt_x, *pt_y, *pt_z,
          nrml_i,   nrml_j,  nrml_k,   err);
   FWDERR();

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
      ERROR(iBase_INVALID_ENTITY_COUNT, "Mismatched array sizes");
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
      FWDERR();

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

void iGeom_getEntTgntXYZ(iGeom_Instance instance,
                         iBase_EntityHandle entity_handle,
                         double x, double y, double z,
                         double* tgnt_i, double* tgnt_j, double* tgnt_k,
                         int* err) {
   RETURN(iBase_NOT_SUPPORTED);
}

void iGeom_getArrTgntXYZ(iGeom_Instance instance,
                         iBase_EntityHandle const* entity_handles,
                         int entity_handles_size,
                         int storage_order, double const* coordinates,
                         int coordinates_size,
                         double** tangents, int* tangents_allocated,
                         int* tangents_size, int* err) {
   RETURN(iBase_NOT_SUPPORTED);
}

void iGeom_getEntBoundBox(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, double* min_x, double* min_y,
      double* min_z, double* max_x, double* max_y, double* max_z, int* err) {
   ErrorCode rval;
   int type;
   iGeom_getEntType(instance, entity_handle, &type, err);
   FWDERR();

   if (type == 0) {
      iGeom_getVtxCoord(instance, entity_handle, min_x, min_y, min_z, err);
      FWDERR();
      max_x = min_x;
      max_y = min_y;
      max_z = min_z;
   } else if (type == 1) {
     // it could be relatively easy to support
      *err = iBase_NOT_SUPPORTED;
      FWDERR();
   } else if (type == 2 || type == 3) {

      EntityHandle root;
      CartVect center, axis[3];
      GeomTopoTool * gtt = GETGTT(instance);
      if (!gtt)
        ERROR(iBase_FAILURE, "Can't get geom topo tool.");
      rval = gtt->get_root(MBH_cast(entity_handle), root);
      CHKERR(rval, "Failed to get tree root in iGeom_getEntBoundBox.");
      rval = gtt->obb_tree()->box(root, center.array(),
            axis[0].array(), axis[1].array(), axis[2].array());
      CHKERR(rval, "Failed to get box from obb tree.");

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
      FWDERR();

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
   ErrorCode rval = FBE_cast(instance)->getVtxCoord(MBH_cast(vertex_handle), x, y, z);
   CHKERR(rval, "Failed to vertex position");
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
  // storage order is ignored
  std::vector<EntityHandle> intersect_handles;
  std::vector<double> coords;
  std::vector<double> params;
  ErrorCode rval = FBE_cast(instance)->getPntRayIntsct(x, y, z, dir_x,
       dir_y, dir_z,intersect_handles, coords, params);
  CHKERR(rval,"can't get ray intersections ");
  *intersect_entity_handles_size = (int)intersect_handles.size();

  CHECK_SIZE(*intersect_entity_handles, *intersect_entity_handles_allocated,
      *intersect_entity_handles_size, iBase_EntityHandle, NULL);
  *intersect_coords_size = 3*(int)intersect_handles.size();
  CHECK_SIZE(*intersect_coords, *intersect_coords_allocated,
      *intersect_coords_size, double, NULL);
  *param_coords_size=(int)intersect_handles.size();
  CHECK_SIZE(*param_coords, *param_coords_allocated,
      *param_coords_size, double, NULL);

  COPY_RANGE(intersect_handles, *intersect_entity_handles);

  COPY_DOUBLEVEC(params, *param_coords);
  if (storage_order == iBase_BLOCKED) {
    int sz=(int)intersect_handles.size();
    for (int i=0; i<sz; i++)
    {
      *intersect_coords[i]=coords[3*i];
      *intersect_coords[sz+i]=coords[3*i+1];
      *intersect_coords[2*sz+i]=coords[3*i+2];
    }
  } else {
    COPY_DOUBLEVEC(coords, *intersect_coords);
  }

   RETURN(iBase_SUCCESS);
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
  // not implemented
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
  ErrorCode rval = FBE_cast(instance)->getEgFcSense(MBH_cast(edge), MBH_cast (face),
     *sense_out);

  CHKERR(rval, "Failed to get edge senses in iGeom_getEgFcSense.");
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

  ErrorCode rval = FBE_cast(instance)->getEgVtxSense(MBH_cast(edge), MBH_cast(vertex1),
      MBH_cast(vertex2), *sense_out);
  CHKERR(rval, "Failed to get vertex sense wrt edge in iGeom_getEgVtxSense");
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
  ErrorCode rval = FBE_cast(instance)->measure((EntityHandle *) (entity_handles) ,
      entity_handles_size,  *measures);
  CHKERR(rval, "Failed to get measures");
  RETURN(iBase_SUCCESS);
}

void iGeom_getFaceType(iGeom_Instance, iBase_EntityHandle face_handle,
      char* face_type, int* err, int* face_type_length) {
}
void iGeom_getParametric(iGeom_Instance instance, int* is_parametric, int* err) {
  *is_parametric = 0; //(false)
  RETURN(iBase_SUCCESS);
}
void iGeom_isEntParametric(iGeom_Instance instance, iBase_EntityHandle entity_handle,
      int* parametric, int* err) {
  int type = -1;
  iGeom_getEntType(instance, entity_handle, &type, err);
  if (type==1)
    *parametric = 1;// true
  else
    *parametric = 0; // false
  RETURN(iBase_SUCCESS);
}
void iGeom_isArrParametric(iGeom_Instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int** is_parametric, int* is_parametric_allocated,
      int* is_parametric_size, int* err) {
  // not implemented
}
void iGeom_getEntUVtoXYZ(iGeom_Instance instance, iBase_EntityHandle entity_handle,
      double u, double v, double* x, double* y, double* z, int* err) {
  RETURN(iBase_NOT_SUPPORTED);
}
void iGeom_getArrUVtoXYZ(iGeom_Instance ,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      int storage_order, double const* uv, int uv_size, double** coordinates,
      int* coordinates_allocated, int* coordinates_size, int* err) {
}


void iGeom_getEntUtoXYZ(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, double u, double* x, double* y,
      double* z, int* err) {
   int type ;
   iGeom_getEntType(instance, entity_handle, &type, err);
   FWDERR();

   if (type != 1)  // not edge
     RETURN(iBase_NOT_SUPPORTED);

   ErrorCode rval = FBE_cast(instance)->getEntUtoXYZ(
       (EntityHandle) entity_handle, u, *x, *y, *z );
   CHKERR(rval, "Failed to get position from parameter");
   RETURN(iBase_SUCCESS);
}

void iGeom_getArrUtoXYZ(iGeom_Instance,
      iBase_EntityHandle const* entity_handles, int entity_handles_size,
      double const* u, int u_size, int storage_order, double** on_coords,
      int* on_coords_allocated, int* on_coords_size, int* err) {
  // not implemented
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
  ErrorCode rval = FBE_cast(instance)->getEntURange((EntityHandle) entity_handle,
                 *u_min,  *u_max );
  CHKERR(rval, "Failed to get range");
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
  // it means deleting some sets from moab db ; is this what we want?
}

void iGeom_deleteEnt(iGeom_Instance instance, iBase_EntityHandle entity_handle,
      int* err) {
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
// start copy old

void iGeom_createEntSet(iGeom_Instance instance, int isList,
      iBase_EntitySetHandle* entity_set_created, int *err) {
   iMesh_createEntSet(IMESH_INSTANCE(instance), isList, entity_set_created, err);
   FWDERR();
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
}

void iGeom_getEntSets(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set_handle, int num_hops,
      iBase_EntitySetHandle** contained_set_handles,
      int* contained_set_handles_allocated, int* contained_set_handles_size,
      int *err) {
   iMesh_getEntSets(IMESH_INSTANCE(instance), entity_set_handle, num_hops,
         contained_set_handles, contained_set_handles_allocated,
         contained_set_handles_size, err);
}

void iGeom_addEntToSet(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, iBase_EntitySetHandle entity_set,
      int *err) {
   iMesh_addEntToSet(IMESH_INSTANCE(instance), entity_handle, entity_set, err);
}

void iGeom_rmvEntFromSet(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, iBase_EntitySetHandle entity_set,
      int *err) {
   iMesh_rmvEntFromSet(IMESH_INSTANCE(instance), entity_handle, entity_set, err);
}

void iGeom_addEntArrToSet(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_EntitySetHandle entity_set, int *err) {
   iMesh_addEntArrToSet(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, entity_set, err);
}

void iGeom_rmvEntArrFromSet(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_EntitySetHandle entity_set, int *err) {
   iMesh_rmvEntArrFromSet(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, entity_set, err);
}

void iGeom_addEntSet(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set_to_add,
      iBase_EntitySetHandle entity_set_handle, int *err) {
   iMesh_addEntSet(IMESH_INSTANCE(instance), entity_set_to_add,
         entity_set_handle, err);
}

void iGeom_rmvEntSet(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set_to_remove,
      iBase_EntitySetHandle entity_set_handle, int *err) {
   iMesh_rmvEntSet(IMESH_INSTANCE(instance), entity_set_to_remove,
         entity_set_handle, err);
}

void iGeom_isEntContained(iGeom_Instance instance,
      iBase_EntitySetHandle containing_entity_set,
      iBase_EntityHandle contained_entity, int *is_contained, int *err) {
   iMesh_isEntContained(IMESH_INSTANCE(instance), containing_entity_set,
         contained_entity, is_contained, err);
}

void iGeom_isEntArrContained(iGeom_Instance instance,
      iBase_EntitySetHandle containing_set,
      const iBase_EntityHandle* entity_handles, int num_entity_handles,
      int** is_contained, int* is_contained_allocated, int* is_contained_size,
      int* err) {
   iMesh_isEntArrContained(IMESH_INSTANCE(instance), containing_set,
         entity_handles, num_entity_handles, is_contained,
         is_contained_allocated, is_contained_size, err);
}

void iGeom_isEntSetContained(iGeom_Instance instance,
      iBase_EntitySetHandle containing_entity_set,
      iBase_EntitySetHandle contained_entity_set, int *is_contained, int *err) {
   iMesh_isEntSetContained(IMESH_INSTANCE(instance), containing_entity_set,
         contained_entity_set, is_contained, err);
}

void iGeom_addPrntChld(iGeom_Instance instance,
      iBase_EntitySetHandle parent_entity_set,
      iBase_EntitySetHandle child_entity_set, int *err) {
   iMesh_addPrntChld(IMESH_INSTANCE(instance), parent_entity_set,
         child_entity_set, err);
}

void iGeom_rmvPrntChld(iGeom_Instance instance,
      iBase_EntitySetHandle parent_entity_set,
      iBase_EntitySetHandle child_entity_set, int *err) {
   iMesh_rmvPrntChld(IMESH_INSTANCE(instance), parent_entity_set,
         child_entity_set, err);
}

void iGeom_isChildOf(iGeom_Instance instance,
      iBase_EntitySetHandle parent_entity_set,
      iBase_EntitySetHandle child_entity_set, int *is_child, int *err) {
   iMesh_isChildOf(IMESH_INSTANCE(instance), parent_entity_set,
         child_entity_set, is_child, err);
}

void iGeom_getNumChld(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, int num_hops, int *num_child, int *err) {
   iMesh_getNumChld(IMESH_INSTANCE(instance), entity_set, num_hops, num_child,
         err);
}

void iGeom_getNumPrnt(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, int num_hops, int *num_parent, int *err) {
   iMesh_getNumPrnt(IMESH_INSTANCE(instance), entity_set, num_hops, num_parent,
         err);
}

void iGeom_getChldn(iGeom_Instance instance,
      iBase_EntitySetHandle from_entity_set, int num_hops,
      iBase_EntitySetHandle** entity_set_handles,
      int* entity_set_handles_allocated, int* entity_set_handles_size, int *err) {
   iMesh_getChldn(IMESH_INSTANCE(instance), from_entity_set, num_hops,
         entity_set_handles, entity_set_handles_allocated,
         entity_set_handles_size, err);
}

void iGeom_getPrnts(iGeom_Instance instance,
      iBase_EntitySetHandle from_entity_set, int num_hops,
      iBase_EntitySetHandle** entity_set_handles,
      int* entity_set_handles_allocated, int* entity_set_handles_size, int *err) {
   iMesh_getPrnts(IMESH_INSTANCE(instance), from_entity_set, num_hops,
         entity_set_handles, entity_set_handles_allocated,
         entity_set_handles_size, err);
}

void iGeom_createTag(iGeom_Instance instance, const char* tag_name,
      int tag_size, int tag_type, iBase_TagHandle* tag_handle, int *err,
      int tag_name_len) {

   iMesh_createTag(IMESH_INSTANCE(instance), tag_name, tag_size, tag_type,
         tag_handle, err, tag_name_len);
}


void iGeom_destroyTag(iGeom_Instance instance, iBase_TagHandle tag_handle,
      int forced, int *err) {
  ErrorCode rval = MBI->tag_delete(TAG_HANDLE(tag_handle));
  CHKERR(rval, "Failed to delete tag");
  RETURN(iBase_SUCCESS);
}

void iGeom_getTagName(iGeom_Instance instance, iBase_TagHandle tag_handle,
      char *name, int* err, int name_len) {
   iMesh_getTagName(IMESH_INSTANCE(instance), tag_handle, name, err, name_len);
}

void iGeom_getTagSizeValues(iGeom_Instance instance,
      iBase_TagHandle tag_handle, int *tag_size, int *err) {
   iMesh_getTagSizeValues(IMESH_INSTANCE(instance), tag_handle, tag_size, err);
}

void iGeom_getTagSizeBytes(iGeom_Instance instance, iBase_TagHandle tag_handle,
      int *tag_size, int *err) {
   iMesh_getTagSizeBytes(IMESH_INSTANCE(instance), tag_handle, tag_size, err);
}

void iGeom_getTagHandle(iGeom_Instance instance, const char* tag_name,
      iBase_TagHandle *tag_handle, int *err, int tag_name_len) {
   iMesh_getTagHandle(IMESH_INSTANCE(instance), tag_name, tag_handle, err,
         tag_name_len);
}

void iGeom_getTagType(iGeom_Instance instance, iBase_TagHandle tag_handle,
      int *tag_type, int *err) {
   iMesh_getTagType(IMESH_INSTANCE(instance), tag_handle, tag_type, err);
}

void iGeom_setEntSetData(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set_handle, iBase_TagHandle tag_handle,
      const void* tag_value, int tag_value_size, int *err) {
   iMesh_setEntSetData(IMESH_INSTANCE(instance), entity_set_handle, tag_handle,
         tag_value, tag_value_size, err);
}

void iGeom_setEntSetIntData(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, iBase_TagHandle tag_handle,
      int tag_value, int *err) {
   iMesh_setEntSetIntData(IMESH_INSTANCE(instance), entity_set, tag_handle,
         tag_value, err);
}

void iGeom_setEntSetDblData(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, iBase_TagHandle tag_handle,
      double tag_value, int *err) {
   iMesh_setEntSetDblData(IMESH_INSTANCE(instance), entity_set, tag_handle,
         tag_value, err);
}

void iGeom_setEntSetEHData(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, iBase_TagHandle tag_handle,
      iBase_EntityHandle tag_value, int *err) {
   iMesh_setEntSetEHData(IMESH_INSTANCE(instance), entity_set, tag_handle,
         tag_value, err);
}

void iGeom_setEntSetESHData(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, iBase_TagHandle tag_handle,
      iBase_EntitySetHandle tag_value, int *err) {
   iMesh_setEntSetESHData(IMESH_INSTANCE(instance), entity_set, tag_handle,
         tag_value, err);
}

void iGeom_getEntSetData(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set_handle, iBase_TagHandle tag_handle,
      void** tag_value, int* tag_value_allocated, int* tag_value_size, int *err) {
   iMesh_getEntSetData(IMESH_INSTANCE(instance), entity_set_handle, tag_handle,
         tag_value, tag_value_allocated, tag_value_size, err);
}

void iGeom_getEntSetIntData(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, iBase_TagHandle tag_handle,
      int *out_data, int *err) {
   iMesh_getEntSetIntData(IMESH_INSTANCE(instance), entity_set, tag_handle,
         out_data, err);
}

void iGeom_getEntSetDblData(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, iBase_TagHandle tag_handle,
      double *out_data, int *err) {
   iMesh_getEntSetDblData(IMESH_INSTANCE(instance), entity_set, tag_handle,
         out_data, err);
}

void iGeom_getEntSetEHData(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, iBase_TagHandle tag_handle,
      iBase_EntityHandle *out_data, int *err) {
   iMesh_getEntSetEHData(IMESH_INSTANCE(instance), entity_set, tag_handle,
         out_data, err);
}

void iGeom_getEntSetESHData(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set, iBase_TagHandle tag_handle,
      iBase_EntitySetHandle *out_data, int *err) {
   iMesh_getEntSetESHData(IMESH_INSTANCE(instance), entity_set, tag_handle,
         out_data, err);
}

void iGeom_getAllEntSetTags(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set_handle, iBase_TagHandle** tag_handles,
      int* tag_handles_allocated, int* tag_handles_size, int *err) {
   iMesh_getAllEntSetTags(IMESH_INSTANCE(instance), entity_set_handle,
         tag_handles, tag_handles_allocated, tag_handles_size, err);
}

void iGeom_rmvEntSetTag(iGeom_Instance instance,
      iBase_EntitySetHandle entity_set_handle, iBase_TagHandle tag_handle,
      int *err) {
   iMesh_rmvEntSetTag(IMESH_INSTANCE(instance), entity_set_handle, tag_handle,
         err);
}

void iGeom_getArrData(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, void** tag_values, int* tag_values_allocated,
      int* tag_values_size, int *err) {
   iMesh_getArrData(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, tag_values, tag_values_allocated,
         tag_values_size, err);
}

void iGeom_getIntArrData(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, int** tag_values, int* tag_values_allocated,
      int* tag_values_size, int *err) {
   iMesh_getIntArrData(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, tag_values, tag_values_allocated,
         tag_values_size, err);
}

void iGeom_getDblArrData(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, double** tag_values,
      int* tag_values_allocated, int* tag_values_size, int *err) {
   iMesh_getDblArrData(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, tag_values, tag_values_allocated,
         tag_values_size, err);
}

void iGeom_getEHArrData(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, iBase_EntityHandle** tag_value,
      int* tag_value_allocated, int* tag_value_size, int *err) {
   iMesh_getEHArrData(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, tag_value, tag_value_allocated,
         tag_value_size, err);
}

void iGeom_getESHArrData(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, iBase_EntitySetHandle** tag_value,
      int* tag_value_allocated, int* tag_value_size, int *err) {
   iMesh_getESHArrData(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, tag_value, tag_value_allocated,
         tag_value_size, err);
}

void iGeom_setArrData(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, const void* tag_values, int tag_values_size,
      int *err) {
   iMesh_setArrData(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, tag_values, tag_values_size, err);
}

void iGeom_setIntArrData(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, const int* tag_values, int tag_values_size,
      int *err) {
   iMesh_setIntArrData(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, tag_values, tag_values_size, err);
}

void iGeom_setDblArrData(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, const double* tag_values,
      const int tag_values_size, int *err) {
   iMesh_setDblArrData(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, tag_values, tag_values_size, err);
}

void iGeom_setEHArrData(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, const iBase_EntityHandle* tag_values,
      int tag_values_size, int *err) {
   iMesh_setEHArrData(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, tag_values, tag_values_size, err);
}

void iGeom_setESHArrData(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, const iBase_EntitySetHandle* tag_values,
      int tag_values_size, int *err) {
   iMesh_setESHArrData(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, tag_values, tag_values_size, err);
}

void iGeom_rmvArrTag(iGeom_Instance instance,
      const iBase_EntityHandle* entity_handles, int entity_handles_size,
      iBase_TagHandle tag_handle, int *err) {
   iMesh_rmvArrTag(IMESH_INSTANCE(instance), entity_handles,
         entity_handles_size, tag_handle, err);
}

void iGeom_getData(iGeom_Instance instance, iBase_EntityHandle entity_handle,
      iBase_TagHandle tag_handle, void** tag_value, int *tag_value_allocated,
      int *tag_value_size, int *err) {
   iMesh_getData(IMESH_INSTANCE(instance), entity_handle, tag_handle,
         tag_value, tag_value_allocated, tag_value_size, err);
}

void iGeom_getIntData(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, iBase_TagHandle tag_handle,
      int *out_data, int *err) {
   iMesh_getIntData(IMESH_INSTANCE(instance), entity_handle, tag_handle,
         out_data, err);
}

void iGeom_getDblData(iGeom_Instance instance,
      const iBase_EntityHandle entity_handle, const iBase_TagHandle tag_handle,
      double *out_data, int *err) {
   iMesh_getDblData(IMESH_INSTANCE(instance), entity_handle, tag_handle,
         out_data, err);
}

void iGeom_getEHData(iGeom_Instance instance, iBase_EntityHandle entity_handle,
      iBase_TagHandle tag_handle, iBase_EntityHandle *out_data, int *err) {
   iMesh_getEHData(IMESH_INSTANCE(instance), entity_handle, tag_handle,
         out_data, err);
}

void iGeom_getESHData(iGeom_Instance instance, iBase_EntityHandle entity_handle,
      iBase_TagHandle tag_handle, iBase_EntitySetHandle *out_data, int *err) {
   iMesh_getESHData(IMESH_INSTANCE(instance), entity_handle, tag_handle,
         out_data, err);
}

void iGeom_setData(iGeom_Instance instance, iBase_EntityHandle entity_handle,
      iBase_TagHandle tag_handle, const void* tag_value, int tag_value_size,
      int *err) {
   iMesh_setData(IMESH_INSTANCE(instance), entity_handle, tag_handle,
         tag_value, tag_value_size, err);
}

void iGeom_setIntData(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, iBase_TagHandle tag_handle,
      int tag_value, int *err) {
   iMesh_setIntData(IMESH_INSTANCE(instance), entity_handle, tag_handle,
         tag_value, err);
}

void iGeom_setDblData(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, iBase_TagHandle tag_handle,
      double tag_value, int *err) {
   iMesh_setDblData(IMESH_INSTANCE(instance), entity_handle, tag_handle,
         tag_value, err);
}

void iGeom_setEHData(iGeom_Instance instance, iBase_EntityHandle entity_handle,
      iBase_TagHandle tag_handle, iBase_EntityHandle tag_value, int *err) {
   iMesh_setEHData(IMESH_INSTANCE(instance), entity_handle, tag_handle,
         tag_value, err);
}

void iGeom_setESHData(iGeom_Instance instance, iBase_EntityHandle entity_handle,
      iBase_TagHandle tag_handle, iBase_EntitySetHandle tag_value, int *err) {
   iMesh_setESHData(IMESH_INSTANCE(instance), entity_handle, tag_handle,
         tag_value, err);
   }

void iGeom_getAllTags(iGeom_Instance instance,
      iBase_EntityHandle entity_handle, iBase_TagHandle** tag_handles,
      int* tag_handles_allocated, int* tag_handles_size, int *err) {
   iMesh_getAllTags(IMESH_INSTANCE(instance), entity_handle, tag_handles,
         tag_handles_allocated, tag_handles_size, err);
   }

void iGeom_rmvTag(iGeom_Instance instance, iBase_EntityHandle entity_handle,
      iBase_TagHandle tag_handle, int *err) {
   iMesh_rmvTag(IMESH_INSTANCE(instance), entity_handle, tag_handle, err);
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


