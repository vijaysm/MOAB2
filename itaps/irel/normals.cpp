/*
 * normals.cpp
 *
 *  Created on: Oct 30, 2015
 *     will take a cubit model and compute normals at all vertices, using exact geometrry info
 */


#include "iGeom.h"
#include "iMesh.h"
#include "iRel.h"

#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "MBiMesh.hpp"
#include "moab/Core.hpp"


#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)
#ifdef MESHDIR
#ifdef HAVE_OCC
#define DEFAULT_GEOM STRINGIFY(MESHDIR/brick.stp)
#define DEFAULT_MESH STRINGIFY(MESHDIR/brick.h5m)
#else
#define DEFAULT_GEOM STRINGIFY(MESHDIR/cyl2.cub)
#define DEFAULT_MESH STRINGIFY(MESHDIR/cyl2.cub)
#endif
#else
#error Specify MESHDIR to compile test
#endif

#define CHECK_SIZE_C(type, array, allocated_size, size)  \
          if (NULL == *array || *allocated_size == 0) {\
            *array = (type *) malloc(sizeof(type) * size); \
            *allocated_size = size;} \
          else if (*allocated_size < size) { \
             printf("   Array passed in is non-zero but too short.\n"); }

typedef void* iRel_EntityHandle;


int print_geom_info(iGeom_Instance geom, iBase_EntityHandle gent)
{
    /* print information about this entity */
  int ent_type;
  int result;
  const char *type_names[] = {"Vertex", "Edge", "Face", "Region"};

  iGeom_getEntType(geom, gent, &ent_type, &result);

  if (iBase_SUCCESS != result) {
    printf("Trouble getting entity adjacencies or types.");
    return 0;
  }

  printf("%s 0x%lx\n", type_names[ent_type], (unsigned long)gent);

  return 1;
}

int print_mesh_info(iMesh_Instance mesh, iBase_EntityHandle ment)
{
    /* print information about this entity */

    /* get adjacencies first; assume not more than 50 */
  iBase_EntityHandle adj_ents[500], *adj_ents_ptr = adj_ents;
   int ent_types[500], *ent_types_ptr = ent_types;
  int adj_ents_alloc = 500, adj_ents_size, ent_types_size,
    ent_types_allocated = 500;
  int result;

  iBase_TagHandle *ment_tags = NULL;
  int ment_tags_size, ment_tags_alloc;

  char **tag_names;
  int i;

  const char *type_names[] = {"Vertex", "Edge", "Face", "Region"};

  iMesh_getEntAdj(mesh, ment, iBase_ALL_TYPES,
                  &adj_ents_ptr, &adj_ents_alloc, &adj_ents_size,
                  &result);

  if (iBase_SUCCESS != result) return 0;

    /* put this ent on the end, then get types */
  adj_ents[adj_ents_size] = ment;
  iMesh_getEntArrType(mesh, adj_ents, adj_ents_size+1,
                      &ent_types_ptr, &ent_types_allocated,
                      &ent_types_size, &result);
  if (iBase_SUCCESS != result) {
    printf("Trouble getting entity adjacencies or types.");
    return 0;
  }

    /* get tags on ment */
  iMesh_getAllTags(mesh, ment,
                   &ment_tags, &ment_tags_alloc, &ment_tags_size,
                   &result);

  printf("Trouble getting tags on an entity or their names.");

    /* while we're at it, get all the tag names */

  tag_names = (char **) malloc(ment_tags_size * sizeof(char*));

  for (i = 0; i < ment_tags_size; i++) {
    tag_names[i] = (char*)malloc(120*sizeof(char));
    iMesh_getTagName(mesh, ment_tags[i], tag_names[i], &result, 120);
  }

    /* now print the information */
  printf("%s %ld:\n", type_names[ent_types[ent_types_size-1]], (long)ment);
  printf("Adjacencies:");
  for (i = 0; i < adj_ents_size; i++) {
    if (i > 0) printf(", ");
    printf("%s %ld", type_names[ent_types[i]],
           (long)adj_ents[i]);
  }
  printf("\nTags: \n");
  for (i = 0; i < ment_tags_size; i++) {
    int tag_type;

    printf("%s ", tag_names[i]);
    iMesh_getTagType(mesh, ment_tags[i], &tag_type, &result);
    if (iBase_SUCCESS != result)
      printf("(trouble getting type...)\n");
    else {
      char *dum_handle = NULL;
      int dum_handle_alloc = 0, dum_handle_size = 0;
      int int_data;
      double dbl_data;
      iBase_EntityHandle eh_data;

      switch (tag_type) {
        case iBase_INTEGER:
          iMesh_getIntData(mesh, ment, ment_tags[i], &int_data, &result);
          printf("(Int value=%d)", int_data);
          break;
        case iBase_DOUBLE:
          iMesh_getDblData(mesh, ment, ment_tags[i], &dbl_data, &result);
          printf("(Dbl value=%f)", dbl_data);
          break;
        case iBase_ENTITY_HANDLE:
          iMesh_getEHData(mesh, ment, ment_tags[i], &eh_data, &result);
          printf("(EH value=%ld)", (long)eh_data);
          break;
        case iBase_BYTES:
          iMesh_getData(mesh, ment, ment_tags[i],
                        (void**)&dum_handle, &dum_handle_alloc,
                        &dum_handle_size, &result);
          if (NULL != dum_handle && dum_handle_size > 0)
            printf("(Opaque value=%c)", dum_handle[0]);
          break;
      }
    }

    printf("\n");
  }
  printf("(end tags)\n\n");
  free(ment_tags);

  return 1;
}

/*!
  @li Load a geom and a mesh file
*/
int load_geom_mesh(const char *geom_filename,
                        const char *mesh_filename,
                        iGeom_Instance geom,
                        iMesh_Instance mesh)
{
    /* load a geom */
  int result;
  iGeom_load(geom, geom_filename, 0, &result, strlen(geom_filename), 0);
  if (iBase_SUCCESS != result) {
    printf("ERROR : can not load a geometry\n");
    return 0;
  }

    /* load a mesh */
  iMesh_load(mesh, 0, mesh_filename, 0, &result, strlen(mesh_filename), 0);
  if (iBase_SUCCESS != result) {
    printf("ERROR : can not load a mesh\n");
    return 0;
  }

  return 1;
}

/*!

  @li Create relation between geom and mesh
*/
int create_relation(iRel_Instance assoc,
                         iGeom_Instance geom,
                         iMesh_Instance mesh,
                         iRel_PairHandle *pair)
{
  int result;
  iBase_Instance iface1, iface2;
  int type1, type2;
  int ent_or_set1, ent_or_set2;
  int status1, status2;

  iRel_PairHandle tmp_pair;
  iRel_PairHandle *pair_ptr = &tmp_pair;
  int pairs_alloc = 1, pairs_size;

    /* create an relation, entity to set */
  iRel_createPair(assoc,
                  geom, iRel_ENTITY, iRel_IGEOM_IFACE, iRel_ACTIVE,
                  mesh, iRel_SET,    iRel_IMESH_IFACE, iRel_ACTIVE,
                  pair, &result);
  if (iBase_SUCCESS != result) {
    printf("Couldn't create a new relation.\n");
    return 0;
  }

  iRel_getPairInfo(assoc, *pair,
                   &iface1, &ent_or_set1, &type1, &status1,
                   &iface2, &ent_or_set2, &type2, &status2, &result);
  if (iBase_SUCCESS != result) {
    printf("Couldn't retrieve relation info.\n");
    return 0;
  }
  if (iface1 != geom || ent_or_set1 != iRel_ENTITY ||
      type1 != iRel_IGEOM_IFACE || iface2 != mesh || ent_or_set2 != iRel_SET ||
      type2 != iRel_IMESH_IFACE) {
    printf("Unexpected relation info returned.\n");
    return 0;
  }

  iRel_findPairs(assoc, geom, &pair_ptr, &pairs_alloc, &pairs_size, &result);
  if (iBase_SUCCESS != result) {
    printf("Couldn't find relation pair when querying geom.\n");
    return 0;
  }
  if (pairs_size != 1 || tmp_pair != *pair) {
    printf("Unexpected relation pairs returned when querying geom.\n");
    return 0;
  }

  iRel_findPairs(assoc, mesh, &pair_ptr, &pairs_alloc, &pairs_size, &result);
  if (iBase_SUCCESS != result) {
    printf("Couldn't find relation pair when querying mesh.\n");
    return 0;
  }
  if (pairs_size != 1 || tmp_pair != *pair) {
    printf("Unexpected relation pairs returned when querying mesh.\n");
    return 0;
  }

  return 1;
}

/*!
  @li Check relation between geom and mesh
*/
int relate_geom_mesh(iRel_Instance assoc,
                          iGeom_Instance geom,
                          iMesh_Instance mesh,
                          iRel_PairHandle pair)
{
  int result;

  iBase_EntityHandle *gentities = NULL;
  int gentities_size = 0, gentities_alloc = 0;

  iBase_EntitySetHandle *mentity_handles = NULL;
  int mentity_handles_size = 0, mentity_handles_alloc = 0;

  const char *dim_tag_name = "GEOM_DIMENSION";
  iBase_TagHandle dim_tag_mesh;

  iBase_EntitySetHandle *mentities_vec;
  int mentities_vec_size = 0;
  int i;

  iBase_EntitySetHandle *out_mentities = NULL;
  int out_mentities_size = 0, out_mentities_alloc = 0;

  iBase_EntitySetHandle *out_mentities2 = NULL;
  int out_mentities2_size = 0, out_mentities2_alloc = 0;

  iBase_EntityHandle *out_gentities = NULL;
  int out_gentities_size = 0, out_gentities_alloc = 0;

    /* relate geometry entities with coresponding mesh entity sets */
  iGeom_getEntities(geom, NULL,
                    iBase_VERTEX,
                    &gentities,
                    &gentities_alloc,
                    &gentities_size, &result);
  if (iBase_SUCCESS != result) {
    printf("Failed to get gentities by type in relate_geom_mesh.\n");
    return 0;
  }

  iRel_inferEntArrRelations(assoc, pair,
                            gentities, gentities_size, 0,
                            &result);
  if (iBase_SUCCESS != result) {
    printf("Failed to relate geom entities in relate_geom_mesh.\n");
    return 0;
  }

    /* relate coresponding mesh entity sets for geometry entities */
    /* get 1-dimensional mesh entitysets */
  iMesh_getEntSets(mesh, NULL, 1,
                   &mentity_handles, &mentity_handles_alloc,
                   &mentity_handles_size, &result);
  if (iBase_SUCCESS != result) {
    printf("Problem to get all entity sets.\n");
    return 0;
  }

    /* get geom dimension tags for mesh entitysets */
  iMesh_createTag(mesh, dim_tag_name, 1, iBase_INTEGER,
                  &dim_tag_mesh, &result, 15);
  if (iBase_SUCCESS != result && result != iBase_TAG_ALREADY_EXISTS) {
    printf("Couldn't create geom dim tag for mesh entities.\n");
    return 0;
  }

    /* get 1-dimensional mesh entitysets */
  mentities_vec = (iBase_EntitySetHandle*)
    malloc(mentity_handles_size*sizeof(iBase_EntitySetHandle));
  for (i = 0; i < mentity_handles_size; i++) { /* test */
    int dim;
    iMesh_getEntSetIntData(mesh, mentity_handles[i], dim_tag_mesh,
                           &dim, &result);
    if (iBase_SUCCESS != result)
      continue;

    if (dim == 1)
      mentities_vec[mentities_vec_size++] = mentity_handles[i];
  }

  iRel_inferSetArrRelations(assoc, pair,
                            mentities_vec, mentities_vec_size,
                            1, &result);
  if (iBase_SUCCESS != result) {
    printf("Failed to relate mesh entities in relate_geom_mesh.\n");
    return 0;
  }

    /* relate all geometry and mesh entities */
  iRel_inferAllRelations(assoc, pair, &result);
  if (iBase_SUCCESS != result) {
    printf("Failed to relate all geom and mesh entities in relate_geom_mesh.\n");
    return 0;
  }

    /* reset geom entities list and get all geom entities (prev
       only vertices) */
  free(gentities);
  gentities = NULL;
  gentities_alloc = 0;
  iGeom_getEntities(geom, NULL,
                    iBase_ALL_TYPES,
                    &gentities,
                    &gentities_alloc,
                    &gentities_size, &result);
  if (iBase_SUCCESS != result) {
    printf("Failed to get gentities by type in relate_geom_mesh.\n");
    return 0;
  }

    /* get related mesh entity sets for geometry entities */
  iRel_getEntArrSetArrRelation(assoc, pair,
                               gentities, gentities_size, 0,
                               &out_mentities, &out_mentities_alloc,
                               &out_mentities_size,
                               &result);
  if (iBase_SUCCESS != result) {
    printf("Failed to get geom entities in relate_geom_mesh.\n");
    return 0;
  }

  if (out_mentities_size != gentities_size) {
    printf("Number of input geom entities and output mesh entity sets should be same\n");
    return 0;
  }

    /* now try deleting this relation */
  iRel_rmvEntArrRelation(assoc, pair, gentities, gentities_size, 0, &result);
  if (iBase_SUCCESS != result) {
    printf("Failed to remove relation in relate_geom_mesh.\n");
    return 0;
  }

  iRel_getEntArrSetArrRelation(assoc, pair,
                               gentities, gentities_size, 0,
                               &out_mentities2, &out_mentities2_alloc,
                               &out_mentities2_size,
                               &result);
  if (iBase_SUCCESS == result) {
    printf("Shouldn't have gotten mesh sets in relate_geom_mesh.\n");
    return 0;
  }

    /* restore the relation, since we need it later */
  iRel_setEntArrSetArrRelation(assoc, pair, gentities, gentities_size,
                               out_mentities, out_mentities_size, &result);
  if (iBase_SUCCESS != result) {
    printf("Failed to restore relation in relate_geom_mesh.\n");
    return 0;
  }

    /* get related geometry entities for mesh entity sets */
  iRel_getSetArrEntArrRelation(assoc, pair,
                               out_mentities, out_mentities_size, 1,
                               &out_gentities, &out_gentities_alloc,
                               &out_gentities_size,
                               &result);
  if (iBase_SUCCESS != result) {
    printf("Failed to get mesh entities in relate_geom_mesh.\n");
    return 0;
  }

  if (out_mentities_size != out_gentities_size) {
    printf("Number of input mesh entity sets and output geom entities should be same\n");
    return 0;
  }
  free(mentity_handles);
  mentity_handles = NULL;
  free(gentities);
  gentities=NULL;
  free(mentity_handles);
  mentity_handles=NULL;
  free(out_mentities);
  out_mentities=NULL;
  free(mentities_vec);
  mentities_vec=NULL;
  free(out_gentities);
  out_gentities=NULL;
  return 1;
}

/*!
  @li compute normals onto the given geometry
*/
int compute_normals(iRel_Instance assoc,
                         iGeom_Instance geom,
                         iMesh_Instance mesh,
                         iRel_PairHandle pair)
{
  int result;
  int i;

  iBase_EntityHandle *gentities = NULL;
  int gentities_size = 0, gentities_alloc = 0;

  iBase_EntitySetHandle *out_mentities = NULL;
  int out_mentities_size, out_mentities_alloc = 0;

    /* get all the face geo entities, and find relation to some mesh entity */
  iGeom_getEntities(geom, NULL, iBase_FACE,
                    &gentities, &gentities_alloc,
                    &gentities_size, &result);
  if (iBase_SUCCESS != result) {
    printf("Problem getting all geom entities.\n" );
    return 0;
  }
  for (i = 0; i < gentities_size; i++) {
      print_geom_info(geom, gentities[i]);
  }


  iRel_getEntArrSetArrRelation(assoc, pair,
                               gentities, gentities_size, 0,
                               &out_mentities, &out_mentities_alloc,
                               &out_mentities_size,
                               &result);
    /* might not all be */
  if (iBase_SUCCESS != result) {
    printf("Failed to get mesh entities related to geom entities in compute_normals.\n");
    return 0;
  }

    /* check that they're all non-null */
  if (out_mentities_size != gentities_size) {
    printf("Number of mesh & related geom entities don't match.\n");
    return 0;
  }

    /* check to make sure they're mesh sets; how to do that? */
  for (i = 0; i < out_mentities_size; i++) {
    int is_list;
    iMesh_isList(mesh, (iBase_EntitySetHandle)out_mentities[i], &is_list, &result);
    if (iBase_SUCCESS != result) {
      printf("Entity set returned from classification wasn't valid.\n");
      return 0;
    }
  }

  moab::Interface * mb=reinterpret_cast<MBiMesh*>(mesh)->mbImpl;
  if (!mb)
    return 0;
  moab::ErrorCode rval = MB_SUCCESS;

  moab::Tag idtag;
  rval = mb->tag_get_handle("GLOBAL_ID", idtag);

  for (i = 0; i < gentities_size; i++) {

    EntityHandle meshSet=(EntityHandle)out_mentities[i];
    Range faces;
    rval = mb->get_entities_by_dimension(meshSet, 2, faces); //MB_CHK_ERR(rval);

    Range vertices;
    rval = mb->get_connectivity(faces, vertices);

    std::vector<double> coords;
    coords.resize(vertices.size()*3);
    rval = mb->get_coords(vertices, &coords[0]);

    int id;
    rval = mb->tag_get_data(idtag, &meshSet, 1, &id);
    printf("surface %d  has %d faces and %d vertices \n", id, faces.size(),  vertices.size() );

    /**\brief Get the normal vector AND closest point on an entity(ies) at
     *        given position(s)
     *
     * Get the normal vector AND closest point on an entity(ies) at given
     * position(s).  If either the number of entities or number of coordinate
     * triples is unity, then all points or entities are queried for that entity
     * or point, respectively, otherwise each point corresponds to each entity.
     * storage_order should be a value in the iBase_StorageOrder enum.
     * \param instance iGeom instance handle
     * \param entity_handles Entity(ies) being queried
     * \param entity_handles_size Number of entity(ies) being queried
     * \param storage_order Storage order in near_coordinates array
     * \param near_coordinates Starting coordinates
     * \param near_coordinates_size Number of values in near_coordinates array
     * \param on_coordinates Closest point array
     * \param on_coordinates_allocated Allocated size of closest point array
     * \param on_coordinates_size Occupied size of closest point array
     * \param normals Normal array
     * \param normals_allocated Allocated size of normal array
     * \param normals_size Occupied size of normal array
     * \param *err Pointer to error type returned from function
     */

  /*iGeom_getArrNrmlPlXYZ( iGeom_Instance instance,
                              iBase_EntityHandle const* entity_handles,
                              int entity_handles_size,
                              int storage_order,
                              double const* near_coordinates,
                              int near_coordinates_size,
                              double** on_coordinates,
                              int* on_coordinates_allocated,
                              int* on_coordinates_size,
                              double** normals,
                              int* normals_allocated,
                              int* normals_size,
                              int* err );*/

   /* void iGeom_getEntNrmlXYZ( iGeom_Instance instance,
                                iBase_EntityHandle entity_handle,
                                double x,
                                double y,
                                double z,
                                double* nrml_i,
                                double* nrml_j,
                                double* nrml_k,
                                int* err );
    */
    std::vector<double> normals;
    normals.resize(3*vertices.size());
    int sizenor=(int)(3*vertices.size());
    std::vector<double> on_coordinates;
    on_coordinates.resize(3*vertices.size());

    int on_coordinates_allocated= (int)(3*vertices.size());
    int on_coordinates_size;

    int sizeall;
    normals.resize(vertices.size()*3);
    /*iGeom_getArrNrmlPlXYZ( geom, &gentities[i], 1, iBase_INTERLEAVED,
        &coords[0], (int)(vertices.size()*3) , &(&on_coordinates[0]), &on_coordinates_allocated,
                                  &on_coordinates_size,
                                  &(&normals[0]),
                                  &sizenor,
                                  &sizeall,&result);


      if (iBase_SUCCESS != result) {
      printf("Problem getting normals.\n" );
      return 0;
    }*/

    for (int j=0; j<(int) vertices.size(); j++)
    {
      iGeom_getEntNrmlXYZ( geom, gentities[i], coords[3*j], coords[3*j+1],coords[3*j+2],
          &normals[3*j], &normals[3*j+1], &normals[3*j+2], &result );
      if (iBase_SUCCESS != result) {
        printf("Problem getting normals.\n" );
        return 0;
      }
      EntityHandle vertex = vertices[j];
      rval = mb->tag_get_data(idtag, &vertex, 1, &id);
      printf("v: %4d  coor: %12.6f %12.6f %12.6f norm:  %12.6f %12.6f %12.6f \n",
       id, coords[3*j] , coords[3*j+1], coords[3*j+2] ,normals[3*j], normals[3*j+1] ,
       normals[3*j+2]);;
    }


  }

  free(gentities);
  gentities=NULL;
  free(out_mentities);
  out_mentities=NULL;
    /* ok, we're done */
  return 1;
}


int main( int argc, char *argv[] )
{
    /* Check command line arg */
  const char *geom_filename = DEFAULT_GEOM;
  const char *mesh_filename = DEFAULT_MESH;

  int result;

  iGeom_Instance geom;
  iMesh_Instance mesh;
  iRel_Instance assoc;
  iRel_PairHandle pair;

  printf("Usage: %s %s\n", geom_filename, mesh_filename);

  if (argc == 2 && !strcmp(argv[1], "-h")) {
    printf("Usage: %s <geom_filename> <mesh_filename>\n",
           argv[0]);
    return 1;
  }
  else if (argc == 2) {
    geom_filename = argv[1];
    mesh_filename = argv[1];
  }
  else if (argc == 3) {
    geom_filename = argv[1];
    mesh_filename = argv[2];
  }

    /* initialize the Geometry */
  iGeom_newGeom(0, &geom, &result, 0);

    /* initialize the Mesh */
  iMesh_newMesh(0, &mesh, &result, 0);

    /* initialize the Associate */
  iRel_create(0, &assoc, &result, 0);

  printf("   load_geom_mesh: \n");
  result = load_geom_mesh(geom_filename, mesh_filename,
                               geom, mesh);

  printf("   create_relation: \n");
  result = create_relation(assoc, geom, mesh, &pair);

  printf("   relate_geom_mesh: \n");
  result = relate_geom_mesh(assoc, geom, mesh, pair);

    /* compute normals */
  printf("   compute normals: \n");
  result = compute_normals(assoc, geom, mesh, pair);

    /* summary */

  iRel_destroy(assoc, &result);
  iMesh_dtor(mesh, &result);
  iGeom_dtor( geom, &result);


  return 0;
}

