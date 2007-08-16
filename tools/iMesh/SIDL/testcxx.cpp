/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

/** ITAPS Mesh Interface Unit Test
 * 
 * This program tests ITAPS mesh interface functions through the SIDL interface.
 * In a nutshell, the test creates (or accesses) an implementation instance,
 * tells it to load a file (specified on the command line), then calls a
 * number of functions evaluating that mesh.  This test also tests mesh, set
 * and tag creation and evaluation functions.
 *
 * Usage: testcxx <mesh_file_name>
 *
 * Compiling
 * ---------
 * 1. Build your server (your implementation of the ITAPS mesh interface) under 
 *    SIDL/Babel and put it in a library in the Babel-generated server directory, 
 *    ${SERVER_DIR}/libimpl.so.  Any include files needed to instantiate this 
 *    server from applications should also be in that directory.  See note a) below.
 * 2. Change the definition of IMPLEMENTATION_CLASS below to be the 
 *    namespace-qualified name of your implementation class.  This definition is
 *    used in the main to instantiate the server using the SIDL _create function, e.g.
 *    .   ITAPSB::Mesh mesh = IMPLEMENTATION_CLASS._create();
 *    (we put the definition of IMPLEMENTATION_CLASS at the top of this file so that
 *    you don't have to understand the rest of the source code to use this test).  See
 *    note b) below.
 * 3. Include the file(s) needed which declare the implementation's namespace and
 *    class
 * 4. Compile this test file to an object file:
 *        g++ -fpic -I. -I${BABEL}/include -I${SERVER_DIR} -c testcxx.cpp
 * 5. Link the application:
 *        g++ -o testcxx testcxx.o -ltest -L${BABEL_LIBS} -lsidl
 * 
 * Notes 
 * -----
 * a) The files don't absolutely need to be arranged as described above.  For
 *    example, some projects put their custom implementation files somewhere besides
 *    ${SERVER_DIR}, to keep them separate from Babel-generated files.
 * b) This test assumes interface instances are created using a direct call to the
 *    SIDL _create() function for that class.  Another way to do this would be to
 *    call a factory linked into the test.  To change the method used to construct
 *    the interface instance for this test, see the use of IMPLEMENTATION_CLASS towards
 *    the end of this file.
 */

//-------------------------------------------------------
// CHANGE THESE DEFINITIONS TO SUIT YOUR IMPLEMENTATION

// #include here any files with declarations needed to instantiate your
// implementation instance
#include "iMesh_SIDL.hh"

// define IMPLEMENTATION_CLASS to be the namespace::class of your implementation
#define IMPLEMENTATION_CLASS "MOAB"
//-------------------------------------------------------


#include <iostream>
#include "iBase.hh"
#include "iMesh.hh"

typedef void* Entity_Handle;
typedef void* EntitySet_Handle;

using namespace std;

#define ARRAY_PTR(array, type) reinterpret_cast<type*>(array._get_ior()->d_firstElement)
#define ARRAY_SIZE(array) (array._is_nil() ? 0 : array.upper(0) - array.lower(0) + 1)

/*!
prints out a result string based on the value of error_code
*/
void handle_error_code(const bool result,
                       int &number_failed,
                       int &/*number_not_implemented*/,
                       int &number_successful)
{
  if (result) {
    cout << "Success";
    number_successful++;
  }
  else {
    cout << "Failure";    
    number_failed++;
  }
}

#define CAST_ITAPS_INTERFACE(var_in, var_out, iface) \
          iBase::iface var_out; \
          try {var_out = var_in;}\
          catch(iBase::Error) {\
            cerr << "Error: current interface doesn't support iface." << endl; \
            return false;}

#define CAST_IMESH_INTERFACE(var_in, var_out, iface) \
          iMesh::iface var_out; \
          try {var_out = var_in;}\
          catch(iBase::Error) {\
            cerr << "Error: current interface doesn't support iface." << endl; \
            return false;}
/*!
@test 
Load Mesh
@li Load a mesh file
*/
bool load_mesh_test(const std::string filename, iMesh::Mesh &mesh)
{
  // load a mesh
  std::string options;
  try {
    mesh.load(0, filename, options);
  }
  catch (iBase::Error) {
    cerr << "ERROR : can not load a mesh" << endl;
    return false;
  }

  return true;
}

/*!
@test
ITAPS topology dimension Test
@li Check 2d topology dimensions
*/
bool topology_dimension_test(iMesh::Mesh &mesh)
{
  // first get 2D entities
  sidl::array<Entity_Handle> faces;
  int faces_size;
  try {
    mesh.getEntities(NULL, iBase::EntityType_FACE, 
                     iMesh::EntityTopology_ALL_TOPOLOGIES, faces, faces_size);
  } catch (iBase::Error err) {
    cerr << "Failed to get faces in entity_sets_test."
	 << endl;
      return false;
  }

  // get dimensions of faces
  //sidl::array<int32_t> dimensions;
  sidl::array<iBase::EntityType> dimensions;
  
  CAST_IMESH_INTERFACE(mesh, mesh_arr, Arr);

  int dim_size;
  try {
    mesh_arr.getEntArrType(faces, ARRAY_SIZE(faces), dimensions, dim_size);
  } catch (iBase::Error err) {
    cerr << "Failed to get dimensions of faces in topology_test."
	 << endl;
    return false;
  }

  if (dim_size != ARRAY_SIZE(faces)) {
    cerr << "Didn't get the right number of types in topology_test."
         << endl;
    return false;
  }
    
  // make sure all elements are 2D
  int num_faces = ARRAY_SIZE(faces);

  for (int i = 0; i < num_faces; i++) {
    if (dimensions.get(i) != iBase::EntityType_FACE) {
      return false;
    }
  }

  return true;
}

/*!
@test
ITAPS topology adjacency Test
@li Check topology information
@li Check adjacency
@li Get interior and exterior elements
*/
// make each topological entity vectors, check their topology
// types, get interior and exterior faces of hexes
bool topology_adjacency_test(iMesh::Mesh &mesh)
{
  CAST_IMESH_INTERFACE(mesh, mesh_arr, Arr);
    //CAST_IMESH_INTERFACE(mesh, mesh_aesq, AdvancedEntitySetQuery);

  int top = iMesh::EntityTopology_POINT;
  int num_test_top = iMesh::EntityTopology_ALL_TOPOLOGIES;
  std::vector<Entity_Handle> **entity_vectors = 
    new std::vector<Entity_Handle> *[num_test_top];

  // make the array of vectors having each topology entities
  for (int i = 0; i < num_test_top; i++)
    entity_vectors[i] = new std::vector<Entity_Handle>;

  // fill the vectors of each topology entities
  // like lines vector, polygon vector, triangle vector,
  // quadrilateral, polyhedrron, tet, hex, prism, pyramid,
  // septahedron vectors
  for (; top < num_test_top; top++) {
    sidl::array<Entity_Handle> entities;
    int entities_size;
    try {
      mesh.getEntities(NULL, iBase::EntityType_ALL_TYPES,
                       static_cast<iMesh::EntityTopology>(top),
                       entities, entities_size);
    } catch (iBase::Error err) {
      cerr << "Failed to get entities in adjacencies_test." << endl;
      return false;
    }

    int num_entities = ARRAY_SIZE(entities);

    if (num_entities > 0) {
      sidl::array<iMesh::EntityTopology> topologies;
      int topologies_size;
      try {
	mesh_arr.getEntArrTopo(entities, num_entities, 
                                 topologies, topologies_size);
      } catch (iBase::Error err) {
	cerr << "Failed to get topologies in adjacencies_test." << endl;
	return false;
      }  

      if (topologies_size != num_entities) {
        cerr << "Didn't get the right number of topologies in topology_adjacency_test."
             << endl;
        return false;
      }
    
        // put entities into vectors of each topology
      for (int j = 0; j < num_entities; j++) {
        bool found = false;
        for (int top_i = iMesh::EntityTopology_POINT;
             top_i < iMesh::EntityTopology_ALL_TOPOLOGIES; top_i++) {
          if (!found && topologies.get(j) == top_i) {
            entity_vectors[top_i]->push_back(entities[j]);
            found = true;
          }
        }
        if (!found)
          cerr << "Didn't find entity type for this topology." << endl;
      }
    }
  }

  // check number of entities for each topology
  for (int top_j = 0; top_j < num_test_top; top_j++) {
    int num_tops = 0;
    try {
      mesh.getNumOfTopo(NULL, static_cast<iMesh::EntityTopology>(top_j), 
                        num_tops);
    } catch (iBase::Error err) {
      cerr << "Failed to get number of topologies in adjacencies_test." << endl;
      return false;
    }
    
    if (static_cast<int>(entity_vectors[top_j]->size()) != num_tops)
      return false;
  }

  // change if 3d topology entities are added or removed
  for (int region_type = iMesh::EntityTopology_TETRAHEDRON;
       region_type < iMesh::EntityTopology_ALL_TOPOLOGIES; region_type++) {
    // get all adjacent faces of regions 
    std::vector<Entity_Handle> region_vector =
      *entity_vectors[region_type];
    
    int num_region = region_vector.size();
    
    if (num_region > 0) {
      int lower[] = {0};
      int upper[] = {num_region - 1};
      sidl::array<Entity_Handle> regions = 
	sidl::array<Entity_Handle>::createCol(1, lower, upper);
      sidl::array<Entity_Handle> adj_faces;
      sidl::array<int> face_offsets;
      
      for (int k = 0; k < num_region; k++) {
        regions.set(k, region_vector[k]);
      }

      int adj_faces_size = 0, face_offsets_size = 0;
      try {
	mesh_arr.getEntArrAdj(regions, num_region, iBase::EntityType_FACE,
                                adj_faces, adj_faces_size, 
                                face_offsets, face_offsets_size);
      } catch (iBase::Error err) {
	cerr << "Failed to get adjacent faces of regions in adjacencies_test." << endl;
	return false;
      }

      if (num_region+1 != face_offsets_size) {
        cerr << "Number of offsets didn't agree with number of regions in topology_adjacency_test."
             << endl;
        return false;
      }

      bool face_loaded = false;

      // check # of faces really loaded
      if (region_type == iMesh::EntityTopology_TETRAHEDRON) {
	if (adj_faces_size == 4*num_region)
	  face_loaded = true;
      }
      else if (region_type == iMesh::EntityTopology_HEXAHEDRON) {
	if (adj_faces_size == 6*num_region)
	  face_loaded = true;
      }
      else if (region_type == iMesh::EntityTopology_PRISM) {
	if (adj_faces_size == 5*num_region)
	  face_loaded = true;
      }
      else if (region_type == iMesh::EntityTopology_PYRAMID) {
	if (adj_faces_size == 5*num_region)
	  face_loaded = true;
      }
      else if (region_type == iMesh::EntityTopology_SEPTAHEDRON) {
	if (adj_faces_size == 7*num_region)
	  face_loaded = true;
      }
      else
	face_loaded = false;

      // get all adjacent regions of adjacent faces
      sidl::array<Entity_Handle> adj_regions;
      int adj_regions_size = 0;
      sidl::array<int> region_offsets;
      int region_offsets_size = 0;
      try {
        mesh_arr.getEntArrAdj(adj_faces, adj_faces_size, iBase::EntityType_REGION,
                              adj_regions, adj_regions_size, 
                              region_offsets, region_offsets_size);
      } catch (iBase::Error err) {
	cerr << "Failed to get regions from faces in adjacencies_test." << endl;
	return false;
      }

      if (adj_faces_size+1 != region_offsets_size) {
        cerr << "Number of offsets didn't agree with number of faces in topology_adjacency_test."
             << endl;
        return false;
      }

      std::vector<Entity_Handle> interior;
      std::vector<Entity_Handle>::iterator iter;
      interior.clear();
      
      // find the interior faces having two adjacent regions
      for (int i = 0; i < adj_faces_size; i++) {
        int next_offset = 0;
        if (i == adj_faces_size-1) next_offset = adj_regions_size;
        else next_offset = region_offsets.get(i+1);

        if (next_offset - region_offsets.get(i) == 2) {
          iter = std::find(interior.begin(), interior.end(), adj_faces.get(i));
          if(iter == interior.end()){
            interior.push_back(adj_faces.get(i));
          }
        }
      }

      // now remove any interior faces from the previous adjacent faces list
      // and we should be left with exterior faces
      std::vector<Entity_Handle> exterior;
      std::vector<Entity_Handle>::iterator jter;
      
      for (int i = 0; i < adj_faces_size; i++) {
	jter = std::find(interior.begin(), interior.end(), adj_faces.get(i));
	if(jter == interior.end())
	  exterior.push_back(adj_faces.get(i));
      }

      int num_faces_per_region = face_offsets.get(1) - face_offsets.get(0);

      // check # of exterior and interior faces
      // should be #exterior_faces + 2 * #interior_faces = #faces_per_region
      // * #regions
      if (face_loaded)
	if (exterior.size()+2*interior.size() !=
	    num_faces_per_region*region_vector.size())
	  return false;
    }
  }
  
  for (int i = 0; i < num_test_top; i++)
    delete entity_vectors[i];

  delete [] entity_vectors;
  
  return true;
}

/*!
@test
ITAPS EntityType Connectivity Test
@li Get coordinates for all type enities
*/

bool entity_connectivity_test(iMesh::Mesh &mesh)
{
  CAST_IMESH_INTERFACE(mesh, mesh_arr, Arr);

  int type = iBase::EntityType_EDGE;

  for (; type < iBase::EntityType_ALL_TYPES; type++) {
    sidl::array<int> offsets;
    int offsets_size;
    sidl::array<int> indices;
    int indices_size;
    sidl::array<iMesh::EntityTopology> topologies;
    int topologies_size;
    
    try {

      mesh.getVtxCoordIndex(NULL, static_cast<iBase::EntityType>(type), 
                            iMesh::EntityTopology_ALL_TOPOLOGIES,
                            iBase::EntityType_VERTEX,
                            offsets, offsets_size, indices, indices_size,
                            topologies, topologies_size);
    } catch (iBase::Error err) {
      cerr << "Failed to get indices of vertices in connectivity_test." << endl;
      return false;
    }

    if (ARRAY_SIZE(offsets) != offsets_size) {
      cerr << "Number of offsets didn't agree with array size in connectivity_test."
           << endl;
      return false;
    }

    if (ARRAY_SIZE(indices) != indices_size) {
      cerr << "Number of indices didn't agree with array size in connectivity_test."
           << endl;
      return false;
    }

    if (ARRAY_SIZE(topologies) != topologies_size) {
      cerr << "Number of topologies didn't agree with array size in connectivity_test."
           << endl;
      return false;
    }

    sidl::array<int> offsets1;
    int offsets1_size;
    sidl::array<int> indices1;
    sidl::array<int> entity_sets;
    int entity_sets_size;
    sidl::array<Entity_Handle> entities;
    int entities_size;

    try {
      mesh.getAdjEntities(NULL, static_cast<iBase::EntityType>(type), 
                          iMesh::EntityTopology_ALL_TOPOLOGIES, iBase::EntityType_VERTEX, 
                          entities, entities_size,
                          offsets1, offsets1_size,
                          entity_sets, entity_sets_size);
    } catch (iBase::Error err) {
      cerr << "Failed to get indices of adjacent entity vertices in connectivity_test." 
	   << endl;
      return false;
    }

    if (ARRAY_SIZE(entities) != entities_size ||
        ARRAY_SIZE(offsets1) != offsets1_size ||
        ARRAY_SIZE(entity_sets) != entity_sets_size) {
      cerr << "Number of elements didn't agree with array size for an array in connectivity_test."
           << endl;
      return false;
    }
  }

  return true;
}

/*!
@test
ITAPS entity sets sub test
@li Check entity sets
*/

// helper function used to report errors in # sets
bool check_esets(iBase::EntSet &mesh_eset, const int num_sets);

bool entity_sets_subtest(iMesh::Mesh &mesh, bool is_list,
                         int num_iter)
{
  int num_type = iBase::EntityType_ALL_TYPES - iBase::EntityType_VERTEX;
  int num_all_entities_super = 0;
  EntitySet_Handle es_array[num_type];
  int number_array[num_type];
  int ent_type = iBase::EntityType_VERTEX;
  
  CAST_IMESH_INTERFACE(mesh, mesh_arr, Arr);
  CAST_IMESH_INTERFACE(mesh, mesh_ent, Entity);
  CAST_ITAPS_INTERFACE(mesh, mesh_eset, EntSet);
  CAST_ITAPS_INTERFACE(mesh, mesh_esbool, SetBoolOps);
  CAST_ITAPS_INTERFACE(mesh, mesh_esrel, SetRelation);

  // get the number of whole mesh
  int n_whole_mesh = 0;
  try {
    mesh_eset.getNumEntSets(NULL, 1, n_whole_mesh);
  } catch (iBase::Error) {
    cerr << "Problem to get the number of all entity sets in whole mesh."
	 << endl;
    return false;
  }

  // add entities to entitysets by type
  for (; ent_type < num_type; ent_type++) {
    // initialize the entityset
    try {
      mesh_eset.createEntSet(is_list, es_array[ent_type]);
    } catch (iBase::Error) {
      cerr << "Problem creating entityset." << endl;
      return false;
    }

    // get entities by type in total "mesh"
    sidl::array<Entity_Handle> entities;
    int entities_size = 0;
    try {
      mesh.getEntities(NULL, static_cast<iBase::EntityType>(ent_type),
                       iMesh::EntityTopology_ALL_TOPOLOGIES, entities, entities_size);
    } catch (iBase::Error err) {
      cerr << "Failed to get entities by type in entity_sets_test."
	   << endl;
      return false;
    }

    if (ARRAY_SIZE(entities) != entities_size) {
      cerr << "Number of entities didn't agree with array size in entity_sets_subtest."
           << endl;
      return false;
    }

    // add entities into entity set
    if (0 != entities_size) {
      try {
        mesh_eset.addEntArrToSet(entities, entities_size, es_array[ent_type]);
      } catch (iBase::Error) {
        cerr << "Failed to add entities in entity_sets_test."
	   << endl;
        return false;
      }
    }
    
    // Check to make sure entity set really has correct number of entities in it
    try {
      mesh.getNumOfType(es_array[ent_type],
                        static_cast<iBase::EntityType>(ent_type), 
                        number_array[ent_type]);
      
    } catch (iBase::Error) {
      cerr << "Failed to get number of entities by type in entity_sets_test."
	   << endl;
      return false;
    }  

    // compare the number of entities by type
    if (number_array[ent_type] != entities_size)
    {
      cerr << "Number of entities by type is not correct"
		<< std::endl;
      return false;
    }

    // add to number of all entities in super set
    num_all_entities_super += entities_size;
  }

  if (!check_esets(mesh_eset, n_whole_mesh + num_type)) return false;

  // make a super set having all entitysets
  EntitySet_Handle super_set = NULL;
  try {
    mesh_eset.createEntSet(is_list, super_set);
  } catch (iBase::Error) {
    cerr << "Failed to create a super set in entity_sets_test."
	 << endl;
    return false;
  }

  try {
    for (int i = 0; i < num_type; i++)
      mesh_eset.addEntSet(es_array[i], super_set);
  } catch (iBase::Error) {
    cerr << "Failed to add a set to a super set in entity_sets_test."
	 << endl;
    return false;
  }

  if (!check_esets(mesh_eset, n_whole_mesh + num_type + 1)) return false;

    //----------TEST BOOLEAN OPERATIONS----------------//

  EntitySet_Handle temp_es1, temp_es2, temp_es3;
  try {
    mesh_eset.createEntSet(is_list, temp_es1);
  } catch (iBase::Error) {
    cerr << "Failed to create a super set in entity_sets_test."
	 << endl;
    return false;
  }
  
  if (!check_esets(mesh_eset, n_whole_mesh + num_type + 2)) return false;

  // Subtract
  // add all EDGEs and FACEs to temp_es1
  // get all EDGE entities
  sidl::array<Entity_Handle> edges;
  sidl::array<Entity_Handle> faces;
  sidl::array<Entity_Handle> temp_entities1;
  sidl::array<Entity_Handle> temp_entities2;
  int edges_size = 0, faces_size = 0, 
    temp_entities1_size = 0, temp_entities2_size = 0;
  

  try {
    mesh.getEntities(es_array[iBase::EntityType_EDGE], iBase::EntityType_EDGE, iMesh::EntityTopology_ALL_TOPOLOGIES, 
                     edges, edges_size);

  } catch (iBase::Error err) {
    cerr << "Failed to get edge entities in entity_sets_test."
	 << endl;
    return false;
  }

  // add EDGEs to es1
  try {
    mesh_eset.addEntArrToSet(edges, edges_size, temp_es1);
  } catch (iBase::Error) {
    cerr << "Failed to add edge entities in entity_sets_test."
	 << endl;
    return false;
  }

  // get all FACE entities
  try {
    mesh.getEntities(es_array[iBase::EntityType_FACE], iBase::EntityType_FACE, 
                     iMesh::EntityTopology_ALL_TOPOLOGIES, 
                     faces, faces_size);
  } catch (iBase::Error err) {
    cerr << "Failed to get edge entities in entity_sets_test."
	 << endl;
    return false;
  }

  // add FACEs to es1
  try {
    mesh_eset.addEntArrToSet(faces, faces_size, temp_es1);
  } catch (iBase::Error) {
    cerr << "Failed to add face entities in entity_sets_test."
	 << endl;
    return false;
  }

  // subtract EDGEs
  
  try {
    mesh_esbool.subtract(temp_es1, es_array[iBase::EntityType_EDGE], temp_es2);
  } catch (iBase::Error) {
    cerr << "Failed to subtract entitysets in entity_sets_test."
	 << endl;
    return false;
  }

  try {
    mesh.getEntities(temp_es2, iBase::EntityType_FACE, iMesh::EntityTopology_ALL_TOPOLOGIES, 
                     temp_entities1, temp_entities1_size);
  } catch (iBase::Error err) {
    cerr << "Failed to get edge entities in entity_sets_test."
	 <<endl;
    return false;
  }

  if (faces_size != temp_entities1_size) {
    cerr << "not match number of entitysets after subtraction \
             in entity_sets_test." << endl;
    return false;
  }

  // check there's nothing but faces in face_es
  sidl::array<iBase::EntityType> types;
  int types_size;
  
  try {
    mesh_arr.getEntArrType(temp_entities1, temp_entities1_size,
                            types, types_size);
  } catch (iBase::Error err) {
    cerr << "Failed to get types of entities in entity_sets_test." << endl;
    return false;
  }
  for (int i = 0; i < types_size; i++) {
    if (types.get(i) != iBase::EntityType_FACE) {
      
      cerr << "wrong entity type for face test in entity_sets_test." << endl;
      return false;
    }
  }

  try {
    mesh_eset.destroyEntSet(temp_es2);
  } catch (iBase::Error err) {
    cerr << "Failed to destroy temp es2." << endl;
    return false;
  }
    
  if (!check_esets(mesh_eset, n_whole_mesh + num_type + 2)) return false;
  
  //------------Intersect------------
  //

  // clean out the temp_ms1
  try {
    mesh_eset.rmvEntArrFromSet(faces, faces_size, temp_es1);
  } catch (iBase::Error) {
    cerr << "Failed to remove face entities in entity_sets_test."
	 << endl;
    return false;
  }

  // check if it is really cleaned out
  int num_rest;
  try {
    mesh.getNumOfType(temp_es1, iBase::EntityType_FACE, num_rest);
  } catch (iBase::Error) {
    cerr << "Failed to get number of entities by type in entity_sets_test."
	 << endl;
    return false;
  }

  if (num_rest != 0) {
    cerr << "failed to remove correctly." << endl;
    return false;
  }
  
  // add EDGEs to temp es1
  try {
    mesh_eset.addEntArrToSet(edges, edges_size, temp_es1);
  } catch (iBase::Error) {
    cerr << "Failed to add edge entities in entity_sets_test."
	 << endl;
    return false;
  }

  // add FACEs to temp es1
  try {
    mesh_eset.addEntArrToSet(faces, faces_size, temp_es1);
  } catch (iBase::Error) {
    cerr << "Failed to add edge entities in entity_sets_test."
	 << endl;
    return false;
  }

  // intersect temp_es1 with edges meshset 
  // temp_ms1 entityset is altered
  try {
    mesh_esbool.intersect(temp_es1, es_array[iBase::EntityType_EDGE], temp_es2);
  } catch (iBase::Error) {
    cerr << "Failed to intersect in entity_sets_test."
	 << endl;
    return false;
  }

  // try to get FACEs, but there should be nothing but EDGE
  try {
    mesh.getEntities(temp_es2, iBase::EntityType_FACE, iMesh::EntityTopology_ALL_TOPOLOGIES, 
                     temp_entities2, temp_entities2_size);
  } catch (iBase::Error err) {
    cerr << "Failed to get face entities in entity_sets_test."
	 <<endl;
    return false;
  }

  if (temp_entities2_size != 0) {
    cerr << "wrong number of faces." << endl;
    return false;
  }

  try {
    mesh_eset.destroyEntSet(temp_es2);
  } catch (iBase::Error err) {
    cerr << "Failed to destroy temp es2." << endl;
    return false;
  }
    
  if (!check_esets(mesh_eset, n_whole_mesh + num_type + 2)) return false;

  //-------------Unite--------------

  // get all regions
  sidl::array<Entity_Handle> regions;
  int regions_size;

  try {
    mesh_eset.createEntSet(is_list, temp_es2);
  } catch (iBase::Error) {
    cerr << "Failed to create a temp entityset in entity_sets_test."
	 << endl;
    return false;
  }

  try {
    mesh.getEntities(es_array[iBase::EntityType_REGION], 
                     iBase::EntityType_REGION, 
                     iMesh::EntityTopology_ALL_TOPOLOGIES, 
                     regions, regions_size);
  } catch (iBase::Error err) {
    cerr << "Failed to get region entities in entity_sets_test." << endl;
    return false;
  }

  // add REGIONs to temp es2
  try {
    mesh_eset.addEntArrToSet(regions, regions_size, temp_es2);
  } catch (iBase::Error) {
    cerr << "Failed to add region entities in entity_sets_test." << endl;
    return false;
  }

  // unite temp_es1 and temp_es2
  try {
    mesh_esbool.unite(temp_es1, temp_es2, temp_es3);
  } catch (iBase::Error) {
    cerr << "Failed to unite in entity_sets_test." << endl;
    return false;
  }

  // perform the check
  int num_regions;
  try {
    mesh.getNumOfType(temp_es3, iBase::EntityType_REGION, num_regions);
  } catch (iBase::Error) {
    cerr << "Failed to get number of region entities by type in entity_sets_test."
	 << endl;
    return false;
  }
  
  if (num_regions != number_array[iBase::EntityType_REGION]) {
    cerr << "different number of regions in entity_sets_test." << endl;
    return false;
  }

  if (!check_esets(mesh_eset, n_whole_mesh + num_type + 4)) return false;

  //--------Test parent/child stuff in entiysets-----------

  // Add 2 meshsets as children to another
  EntitySet_Handle parent_child;
  try {
    mesh_eset.createEntSet(is_list, parent_child);
  } catch (iBase::Error) {
    cerr << "Problem creating entityset in entity_sets_test."
	 << endl;
    return false;
  }

  try {
    mesh_esrel.addPrntChld(es_array[iBase::EntityType_VERTEX], parent_child);
  } catch (iBase::Error) {
    cerr << "Problem add parent in entity_sets_test."
	 << endl;
    return false;
  }

  // check if parent is really added
  sidl::array<Entity_Handle> parents;
  int parents_size = 0;
  try {
    mesh_esrel.getPrnts(parent_child, 0, parents, parents_size);
  } catch (iBase::Error) {
    cerr << "Problem getting parents in entity_sets_test."
	 << endl;
    return false;
  }

  if (parents_size != 1) {
    cerr << "number of parents is not correct in entity_sets_test."
	 << endl;
    return false;
  }

  // get the number of child entitysets
  int temp_numb;
  try {
    mesh_esrel.getNumChld(es_array[iBase::EntityType_VERTEX], 0,
                          temp_numb);
  } catch (iBase::Error) {
    cerr << "Problem getting number of children in entity_sets_test."
	 << endl;
    return false;
  }

  if (temp_numb != 1) {
    cerr << "number of children is not correct in entity_sets_test."
	 << endl;
    return false;
  }

  // parent_child and es_array[iBase::EntityType_VERTEX] should be related
  try {
    int is_child;
    mesh_esrel.isChildOf(es_array[iBase::EntityType_VERTEX], parent_child, 
                         is_child);
    if (!is_child) {
      cerr << "parent_child and es_array[iBase::EntityType_EDGE] should be related" << endl;
      return false;
    }
  } catch (iBase::Error) {
    cerr << "Problem checking relation in entity_sets_test."
	 << endl;
    return false;
  }
    
  // es_array[iBase::EntityType_FACE] and es_array[iBase::EntityType_REGION] are not related
  try {
    int is_child;
    mesh_esrel.isChildOf(es_array[iBase::EntityType_FACE], 
                         es_array[iBase::EntityType_REGION],
                         is_child);
    if (is_child) {
      cerr << "es_array[iBase::EntityType_REGION] and es_array[iBase::EntityType_FACE] should not be related" << endl;
      return false;
    }
  } catch (iBase::Error) {
    cerr << "Problem checking relation in entity_sets_test."
	 << endl;
    return false;
  }
  
  if (!check_esets(mesh_eset, n_whole_mesh + num_type + 5)) return false;

  //--------test modify and query functions-----------------------------
  
  // get all entity sets in super set
  sidl::array<void*> es_array1;
  int es_array1_size;
  try {
    mesh_eset.getEntSets(super_set, 0, es_array1, es_array1_size);
  } catch (iBase::Error) {
    cerr << "Problem to get entity sets in super set."
	 << endl;
    return false;
  }

  int num_super;

  // get the number of entity sets in super set
  try {
    mesh_eset.getNumEntSets(super_set, 0, num_super);
  } catch (iBase::Error) {
    cerr << "Problem to get the number of all entity sets in super set."
	 << endl;
    return false;
  }

  // the number of entity sets in super set should be same
  if (num_super != es_array1_size) {
    cerr << "the number of entity sets in super set should be same." << endl;
    return false;
  }

  // get all entities in super set
  sidl::array<Entity_Handle> all_entities;
  int all_entities_size;
  try {
    mesh.getEntities(super_set, iBase::EntityType_ALL_TYPES,
                     iMesh::EntityTopology_ALL_TOPOLOGIES, 
                     all_entities, all_entities_size);
  } catch (iBase::Error) {
    cerr << "Problem to get all entities in super set."
	 << endl;
    return false;
  }
  
  // compare the number of all entities in super set
  // NOTE: TEST COMMENTED OUT UNTIL RESOLUTION OF WHETHER GETENTITIES
  // SHOULD GET A NUM_HOPS ARGUMENT
  //  if (num_all_entities_super != all_entities_size) {
  //    cerr << "number of all entities in super set should be same." << endl;
  //    return false;
  //  }

  // test add, remove and get all entitiy sets using super set
  // check GetAllEntitySets works recursively and dosen't return
  // multi sets
  for (int k = 0; k < num_super; k++) {
    // add entity sets of super set to each entity set of super set
    // make multiple child super sets
    EntitySet_Handle es_k = es_array1.get(k);
    try {
      for (int l = 0; l < es_array1_size; l++)
        mesh_eset.addEntSet(es_array1.get(l), es_k);
    } catch (iBase::Error) {
      cerr << "Problem to add entity set to entityset."
	   << endl;
      return false;
    }

    // add super set to each entity set
    try {
      mesh_eset.addEntSet(super_set, es_k);
    } catch (iBase::Error) {
      cerr << "Problem to add super set to entitysets." << endl;
      return false;
    }

    // add one entity sets multiple times
    for (int l = 0; l < 3; l++) {
      try {
        mesh_eset.addEntSet(temp_es1, es_k);
      } catch (iBase::Error) {
	cerr << "Problem to add temp set to entitysets." << endl;
	return false;
      }
    }
  }

  // get adjacent face of hexes
  sidl::array<Entity_Handle> adj_faces;
  int adj_faces_size;
  sidl::array<int> face_offsets;
  sidl::array<int> face_in_sets;
  int face_offsets_size, face_in_sets_size;

  try {
    mesh.getAdjEntities(NULL, iBase::EntityType_ALL_TYPES,
                        iMesh::EntityTopology_HEXAHEDRON, iBase::EntityType_FACE,
                        adj_faces, adj_faces_size,
                        face_offsets, face_offsets_size,
                        face_in_sets, face_in_sets_size);
  } catch (iBase::Error) {
    cerr << "Problem to get adjacent entities in entitysets_test."
	 << endl;
    return false;
  }

  // get all hexes and get faces of that hexes
  sidl::array<Entity_Handle> hexes;
  int hexes_size;

  try {
    mesh.getEntities(NULL, iBase::EntityType_ALL_TYPES,
                     iMesh::EntityTopology_HEXAHEDRON, hexes, hexes_size);
  } catch (iBase::Error err) {
    cerr << "Failed to get hexes in entity_sets_test."
	 << endl;
      return false;
  }
  
  EntitySet_Handle hex_set;

  try {
    mesh_eset.createEntSet(false, hex_set);
  } catch (iBase::Error) {
    cerr << "Problem creating entityset in entity_sets_test."
	 << endl;
    return false;
  }
  
  try {
    mesh_eset.addEntArrToSet(hexes, hexes_size, hex_set);
  } catch (iBase::Error) {
    cerr << "Failed to add hexes in entity_sets_test." << endl;
    return false;
  }
  
  // get adjacent faces of all hexes
  sidl::array<Entity_Handle> adj_faces1;
  int adj_faces1_size;
  sidl::array<int> face_offsets1;
  sidl::array<int> face_in_sets1;
  int face_offsets1_size, face_in_sets1_size;

  try {
    mesh.getAdjEntities(hex_set,
                        iBase::EntityType_ALL_TYPES,
                        iMesh::EntityTopology_HEXAHEDRON, iBase::EntityType_FACE,
                        adj_faces1, adj_faces1_size, 
                        face_offsets1, face_offsets1_size,
                        face_in_sets1, face_in_sets1_size);
  } catch (iBase::Error err) {
    cerr << "Failed to get faces from hexes in entityset_test." << endl;
    return false;
  }

  // compare number of faces
  if (adj_faces_size != adj_faces1_size ||
      face_offsets_size != face_offsets1_size ||
      face_in_sets_size != face_in_sets1_size)
    return false;
  
  if (!check_esets(mesh_eset, n_whole_mesh + num_type + 6)) return false;

  return true;
}

bool check_esets(iBase::EntSet &mesh_eset, const int num_sets) 
{
  int entity_sets_size;
  
  try {
    mesh_eset.getNumEntSets(NULL, 1, entity_sets_size);
  } catch (iBase::Error) {
    cerr << "Problem to get all entity sets in mesh."
	 << endl;
    return false;
  }
  if (entity_sets_size != num_sets) {
    cerr << "the number of entity sets in whole mesh should be "
         << num_sets
         << ", actual number is " << entity_sets_size << "." << endl;
    return false;
  }

  return true;
}

/*!
@test
ITAPS entity sets Test
@li Check entity sets
*/
bool entity_sets_test(iMesh::Mesh &mesh)
{
  int iter_num = 0;
  // check 
  for (int i = 0; i < 2; i++) {
      iter_num++;

      bool result = entity_sets_subtest(mesh, i, iter_num);
      if (!result)
	return result;
  }
  
  return true;
}

/*!
@test 
Vertex Coordinates
@li Get coordinates of vertices by 2 different methods then compare them
*/
bool vertex_coordinates_test(iMesh::Mesh &mesh)
{

  // check storage order
  iBase::StorageOrder this_order = iBase::StorageOrder_UNDETERMINED;
  try {
    mesh.getDfltStorage(this_order);
  }
  catch (iBase::Error) {
    std::cout << "failed to get preferred storage order in vertex_coordinates_test."
	      << std::endl;
    return false;
  }

  // get the coordinates in one array
  sidl::array<double> all_coords;
  int all_coords_size;
  sidl::array<int> in_entity_set;
  int in_entity_set_size;

  try {
    mesh.getAllVtxCoords(0, all_coords, all_coords_size,
                         in_entity_set, in_entity_set_size,
                         this_order);
  }
  catch(iBase::Error) {
    std::cout << "failed to get vertex coordinates in vertex_coordinates_test."
              << std::endl;
    return false;
  }

    // make sure all vertices marked as being in the model
  int not_in_set = 0;
  for (int i = 0; i <= ARRAY_SIZE(in_entity_set)-1; i++)
    if (1 != in_entity_set.get(i)) not_in_set++;

  if (0 != not_in_set) {
    std::cout << "vertex_coordinates_test: " << not_in_set 
              << " vertices not marked as being in whole mesh." << std::endl;
    return false;
  }

  // now get the vertex coordinates from a vertex array
  // need to get the vertices in the model
  sidl::array<void*> verts;
  int verts_size;
  try {
    mesh.getEntities(0, iBase::EntityType_VERTEX, 
                     iMesh::EntityTopology_POINT, verts, verts_size);
  }
  catch (iBase::Error) {
    std::cout << "failed to get entities in vertex_coordinates_test."
	      << std::endl;
    return false;
  }

  // get the coordinates in one array
  sidl::array<double> vert_coords;
  int vert_coords_size;

  try {
    mesh.getVtxArrCoords(verts, verts_size, this_order, 
                         vert_coords, vert_coords_size);
  }
  catch (iBase::Error) {
    std::cout << "failed to get vertex cooridinate of entities in vertex_coordinates_test."
	      << std::endl;
    return false;
  }

  // two lists should be the same length
  if (ARRAY_SIZE(all_coords) != ARRAY_SIZE(vert_coords)) {
    std::cout << "vertex_coordinates_test: vertex arrays different lengths." 
              << std::endl;
    return false;
  }

    // compare the two lists, they should be the same
  int num_bad = 0;

    // can produce an error by rounding off coordinates
  for (int i = 0; i <= ARRAY_SIZE(all_coords)-1; i++)
    if (all_coords.get(i) != vert_coords.get(i)) num_bad++;

  if (0 != num_bad) {
    std::cout << "vertex_coordinates_test: different values in vertex coordinate"
              << " arrays" << std::endl;
    return false;
  }
  
    // if we get here, this test was successful
  return true;
}

/*!
@test
ITAPS Tag Info test
@li Tests tagCreate, tagDelete, tagGetName, tagGetSize, tagGetHandle
*/
bool tag_info_test(iMesh::Mesh &mesh)
{
  CAST_ITAPS_INTERFACE(mesh, mesh_tag, Tag);
  
  // Create a tag
  void *tag_handle = NULL;
  std::string tag_name("int_tag");
  try {mesh_tag.createTag(tag_name, 1, iBase::TagValueType_INTEGER, 
                          tag_handle);}
  catch (iBase::Error) {
    std::cout << "Failed to create tag int_tag in vertex_tag_test." << std::endl;
    return false;
  }

    // check tag info functions
  std::string dum_name;
  try {mesh_tag.getTagName(tag_handle, dum_name);}
  catch (iBase::Error) {
    std::cerr << "Couldn't get name of tag just created." << std::endl;
    return false;
  }
  if (dum_name != "int_tag") {
    std::cerr << "Tag names didn't match." << std::endl;
    return false;
  }
  
  int dum_size;
  try {mesh_tag.getTagSizeBytes(tag_handle, dum_size);}
  catch (iBase::Error) {
    std::cerr << "Couldn't get size of tag just created." << std::endl;
    return false;
  }
  if (dum_size != sizeof(int)) {
    std::cerr << "Tag sizes didn't match." << std::endl;
    return false;
  }

  void *dum_handle;
  try {mesh_tag.getTagHandle(tag_name, dum_handle);}
  catch (iBase::Error) {
    std::cerr << "Couldn't get handle of tag just created." << std::endl;
    return false;
  }
  if (dum_handle != tag_handle) {
    std::cerr << "Tag handles didn't match." << std::endl;
    return false;
  }

    // test non-forced version of tagDelete; forced version tested later,
    // when we have entities around
  try {mesh_tag.destroyTag(tag_handle, false);}
  catch (iBase::Error) {
    std::cerr << "Couldn't delete tag just created." << std::endl;
    return false;
  }

    // look for that tag, to make sure it got deleted
  bool error = false;
  try {mesh_tag.getTagHandle(tag_name, dum_handle);}
  catch (iBase::Error) {
    error = true;
  }
  if (!error) {
    std::cerr << "tagGetHandle was able to find handle for deleted tag." << std::endl;
    return false;
  }

  return true;
}

bool vertex_int_tag_test(iMesh::Mesh &mesh, sidl::array<void*> &verts,
                         void *&int_tag) 
{
  CAST_ITAPS_INTERFACE(mesh, mesh_tag, Tag);
  CAST_ITAPS_INTERFACE(mesh, mesh_etag, EntTag);
  CAST_ITAPS_INTERFACE(mesh, mesh_atag, ArrTag);

  // create a tag
  std::string tag_name("int_tag");
  try {mesh_tag.createTag(tag_name, 1, iBase::TagValueType_INTEGER, 
                          int_tag);}
  catch (iBase::Error) {
    std::cout << "Failed to create tag int_tag in vertex_int_tag_test." << std::endl;
    return false;
  }

  // put a value in the first vertex and retrieve
  sidl::array<void*> dum_array = verts.create1d(1);
  int dum_array_size = 1;
  sidl::array<char> dum_val_array = dum_val_array.create1d(sizeof(int));
  int dum_val_array_size = sizeof(int);
  dum_array.set(0, verts.get(0));
  (ARRAY_PTR(dum_val_array, int))[0] = 11;
  try {
    mesh_atag.setArrData(dum_array, dum_array_size, int_tag, 
                         dum_val_array, dum_val_array_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to set int tag (val=11) in vertex_int_tag_test." << std::endl;
    return false;
  }

  try {mesh_atag.getArrData(dum_array, dum_array_size, int_tag,
                            dum_val_array, dum_val_array_size);}
  catch (iBase::Error) {
    std::cout << "Failed to get int tag (val=11) in vertex_int_tag_test." << std::endl;
    return false;
  }

  if(ARRAY_PTR(dum_val_array, int)[0] != 11) {
    
    std::cout << "Value of vertex tag (val=11) wrong." << std::endl;
    return false;
  }

  // put a value in the last vertex and retrieve
  ARRAY_PTR(dum_val_array, int)[0] = 12;

  try {
    mesh_atag.setArrData(dum_array, dum_array_size, int_tag, 
                         dum_val_array, dum_val_array_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to set int tag (val=12) in vertex_int_tag_test." << std::endl;
    return false;
  }

  try {mesh_atag.getArrData(dum_array, dum_array_size, int_tag,
                            dum_val_array, dum_val_array_size);}
  catch (iBase::Error) {
    std::cout << "Failed to get int tag (val=12) in vertex_int_tag_test." << std::endl;
    return false;
  }

  if(ARRAY_PTR(dum_val_array, int)[0] != 12) {
    
    std::cout << "Value of vertex tag (val=12) wrong." << std::endl;
    return false;
  }

    // ok, the int tag test worked
  return true;
}

bool vertex_double_tag_test(iMesh::Mesh &mesh, sidl::array<void*> &verts,
                            void *&double_tag) 
{
  CAST_ITAPS_INTERFACE(mesh, mesh_atag, ArrTag);

  // create a tag
  std::string tag_name("double_tag");
  try {mesh_atag.createTag(tag_name, 1, iBase::TagValueType_DOUBLE, 
                           double_tag);}
  catch (iBase::Error) {
    std::cout << "Failed to create tag double_tag in vertex_double_tag_test." << std::endl;
    return false;
  }

  // put a value in the first vertex and retrieve
  sidl::array<void*> dum_array = verts.create1d(1);
  int dum_array_size = 1;
  sidl::array<char> dum_val_array = dum_val_array.create1d(sizeof(double));
  int dum_val_array_size = ARRAY_SIZE(dum_val_array);
  dum_array.set(0, verts.get(0));
  (ARRAY_PTR(dum_val_array, double))[0] = 1.0e6;
  try {
    mesh_atag.setArrData(dum_array, dum_array_size, double_tag, 
                         dum_val_array, dum_val_array_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to set double tag (val=1.0e6) in vertex_double_tag_test." << std::endl;
    return false;
  }

  try {mesh_atag.getArrData(dum_array, dum_array_size, double_tag,
                            dum_val_array, dum_val_array_size);}
  catch (iBase::Error) {
    std::cout << "Failed to get double tag (val=1.0e6) in vertex_double_tag_test." << std::endl;
    return false;
  }

  if((ARRAY_PTR(dum_val_array, double))[0] != 1.0e6) {
    
    std::cout << "Value of vertex tag (val=1.0e6) wrong." << std::endl;
    return false;
  }

  // put a value in the last vertex and retrieve
  (ARRAY_PTR(dum_val_array, double))[0] = 2.0e9;

  try {
    mesh_atag.setArrData(dum_array, dum_array_size, double_tag, 
                         dum_val_array, dum_val_array_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to set double tag (val=2.0e9) in vertex_double_tag_test." << std::endl;
    return false;
  }

  try {mesh_atag.getArrData(dum_array, dum_array_size, double_tag,
                            dum_val_array, dum_val_array_size);}
  catch (iBase::Error) {
    std::cout << "Failed to get double tag (val=2.0e9) in vertex_double_tag_test." << std::endl;
    return false;
  }

  if((ARRAY_PTR(dum_val_array, double))[0] != 2.0e9) {
    std::cout << "Value of vertex tag (val=2.0e9) wrong." << std::endl;
    return false;
  }

    // ok, the double tag test worked
  return true;
}

bool vertex_struct_tag_test(iMesh::Mesh &mesh, sidl::array<void*> &verts,
                            void *&struct_tag) 
{
  CAST_ITAPS_INTERFACE(mesh, mesh_atag, ArrTag);
      
  // Add a struct Vertex Tag to the database
  struct TagStruct {
    double test_double;
    int test_int1, test_int2;
  };

  // create a tag
  std::string tag_name("struct_tag");
  try {mesh_atag.createTag(tag_name, sizeof(TagStruct), iBase::TagValueType_BYTES, 
                           struct_tag);}
  catch (iBase::Error) {
    std::cout << "Failed to create tag struct_tag in vertex_struct_tag_test." << std::endl;
    return false;
  }

  // put a value in the first vertex and retrieve
  sidl::array<void*> dum_array = verts.create1d(1);
  int dum_array_size = 1;
  dum_array.set(0, verts.get(0));
    // careful setting the value, since tags are opaque
  sidl::array<char> dum_val_array = dum_val_array.create1d(sizeof(TagStruct));
  int dum_val_array_size = sizeof(TagStruct);
  TagStruct dum_struct = {3.0e12, 2, 3};
  (ARRAY_PTR(dum_val_array, TagStruct))[0] = dum_struct;
  
  try {
    mesh_atag.setArrData(dum_array, dum_array_size, struct_tag, 
                         dum_val_array, dum_val_array_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to set struct tag in vertex_struct_tag_test." << std::endl;
    return false;
  }

  try {mesh_atag.getArrData(dum_array, dum_array_size, struct_tag,
                            dum_val_array, dum_val_array_size);}
  catch (iBase::Error) {
    std::cout << "Failed to get struct tag in vertex_struct_tag_test." << std::endl;
    return false;
  }

  if(dum_val_array_size != sizeof(TagStruct)) {
    std::cout << "Size of vertex struct tag wrong." << std::endl;
    return false;
  }
  if((ARRAY_PTR(dum_val_array, TagStruct))[0].test_int1 != dum_struct.test_int1 ||
     (ARRAY_PTR(dum_val_array, TagStruct))[0].test_int2 != dum_struct.test_int2 ||
     (ARRAY_PTR(dum_val_array, TagStruct))[0].test_double != dum_struct.test_double) {
    std::cout << "Value of vertex struct tag wrong." << std::endl;
    return false;
  }

    // ok, the bool tag test worked
  return true;
}

bool vertex_tag_delete_test(iMesh::Mesh &mesh, sidl::array<void*> &verts) 
{
    // test forced, unforced deletion of tags from entities
  CAST_ITAPS_INTERFACE(mesh, mesh_etag, EntTag);
  CAST_ITAPS_INTERFACE(mesh, mesh_atag, ArrTag);

    // test getAllTagHandles for first entity
  sidl::array<void*> all_tags, dum_entities;
  int all_tags_size, dum_entities_size;
  try {mesh_etag.getAllTags(verts.get(0), all_tags, all_tags_size);}
  catch (iBase::Error) {
    std::cerr << "Couldn't get all tag handles from vertex." << std::endl;
    return false;
  }

    // try removing one
  if (all_tags_size > 0) {
    dum_entities = dum_entities.create1d(1);
    dum_entities_size = 1;
    dum_entities.set(0, verts.get(0));
    try {mesh_atag.rmvArrTag(dum_entities, dum_entities_size, all_tags.get(0));}
    catch (iBase::Error) {
      std::cerr << "Couldn't remove tag from vertex." << std::endl;
      return false;
    }
  
    // try deleting another unforced - should not get error
    if (ARRAY_SIZE(all_tags) > 1) {
      bool delete_err = false;
      try {mesh_atag.destroyTag(all_tags.get(1), false);}
      catch (iBase::Error) {
        delete_err = true;
      }
      if (!delete_err) {
        std::cerr << "Error when unforced-deleting tag in use on a vertex." << std::endl;
        return false;
      }

      // now try doing a forced delete
      try {mesh_atag.destroyTag(all_tags.get(1), true);}
      catch (iBase::Error) {
        std::cerr << "Couldn't force-delete a tag in use on a vertex." << std::endl;
        return false;
      }
    }
  }

    // ok, we're done
  return true;
}

bool vertex_tag_test(iMesh::Mesh &mesh)
{
  bool info_err = tag_info_test(mesh);
  
  // get all the vertices
  sidl::array<void*> verts;
  int verts_size;
  try {mesh.getEntities(NULL, iBase::EntityType_ALL_TYPES, 
                        iMesh::EntityTopology_POINT, verts, verts_size);}
  catch(iBase::Error) {
    std::cout << "entitysetGetEntities failed in vertex_tag_test." << std::endl;
    return false;
  }

  void *int_tag, *double_tag, *struct_tag;
  
    // vertex int tag
  bool int_err = vertex_int_tag_test(mesh, verts, int_tag);
  
    // vertex double tag
  bool double_err = vertex_double_tag_test(mesh, verts, double_tag);
  
    // vertex struct tag
  bool struct_err = vertex_struct_tag_test(mesh, verts, struct_tag);

  bool tag_delete_err = vertex_tag_delete_test(mesh, verts);
  
  return (info_err && int_err && double_err && 
          struct_err && tag_delete_err);
}

bool entityset_int_tag_test(iMesh::Mesh &mesh, sidl::array<void*> &faces,
                            void *&int_tag) 
{
  CAST_ITAPS_INTERFACE(mesh, mesh_etag, EntTag);
  
  // create a tag
  std::string tag_name("set_int_tag");
  try {
    mesh_etag.createTag(tag_name, 1, iBase::TagValueType_INTEGER, 
                        int_tag);
  }
  catch (iBase::Error) {
    std::cout << "Failed to create tag int_tag in entityset_int_tag_test." << std::endl;
    return false;
  }

  // put a value in the first set and retrieve
  sidl::array<char> dum_val;
  dum_val = dum_val.create1d(sizeof(int));
  int *dum_val_ptr = ARRAY_PTR(dum_val, int);
  dum_val_ptr[0] = 11;
  int dum_val_size = ARRAY_SIZE(dum_val);

  try {
    mesh_etag.setData(faces.get(0), int_tag, dum_val, sizeof(int));
  }
  catch (iBase::Error) {
    std::cout << "Failed to set int tag (val=11) in entityset_int_tag_test." << std::endl;
    return false;
  }

  dum_val_ptr[0] = 0;
  try {
    mesh_etag.getData(faces.get(0), int_tag,
                      dum_val, dum_val_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to get int tag (val=11) in entityset_int_tag_test."
	      << std::endl;
    return false;
  }

  if (dum_val_ptr[0] != 11) {
    std::cout << "Value of entityset tag (val=11) wrong." << std::endl;
    return false;
  }

  // put a value in the last faces and retrieve
  dum_val_ptr[0] = 12;
  try {
    mesh_etag.setData(faces.get(ARRAY_SIZE(faces)-1), int_tag,
                      dum_val, dum_val_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to set int tag (val=12) in entityset_int_tag_test."
	      << std::endl;
    return false;
  }

  try {
    mesh_etag.getData(faces.get(ARRAY_SIZE(faces)-1), int_tag,
                      dum_val, dum_val_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to get int tag (val=12) in entityset_int_tag_test."
	      << std::endl;
    return false;
  }
  
  if(dum_val_ptr[0] != 12) {
    std::cout << "Value of entityset tag (val=12) wrong." << std::endl;
    return false;
  }

  // ok, the int tag test worked
  return true;
}

bool entityset_double_tag_test(iMesh::Mesh &mesh, sidl::array<void*> &faces,
			       void *&double_tag) 
{
  CAST_ITAPS_INTERFACE(mesh, mesh_etag, EntTag);

  // create a tag
  std::string tag_name("set_double_tag");
  try {mesh_etag.createTag(tag_name, 1, iBase::TagValueType_DOUBLE, double_tag);}
  catch (iBase::Error) {
    std::cout << "Failed to create tag set_double_tag in entityset_double_tag_test." << std::endl;
    return false;
  }

  sidl::array<char> dum_val;
  int dum_val_size = sizeof(double);
  dum_val = dum_val.create1d(dum_val_size);
  double *dum_val_ptr = ARRAY_PTR(dum_val, double);
  dum_val_ptr[0] = 1.0e6;

  try {
    mesh_etag.setData(faces.get(0), double_tag, dum_val, dum_val_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to set double tag (val=1.0e6) in entityset_double_tag_test."
	      << std::endl;
    return false;
  }

  dum_val_ptr[0] = 0.0;
  try {
    mesh_etag.getData(faces.get(0), double_tag,
                      dum_val, dum_val_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to get double tag (val=1.0e6) in entityset_double_tag_test."
	      << std::endl;
    return false;
  }

  if (dum_val_ptr[0] != 1.0e6) {
    std::cout << "Value of entityset tag (val=1.0e6) wrong." 
	      << std::endl;
    return false;
  }

  // put a value in the last vertex and retrieve
  dum_val_ptr[0] = 2.0e9;
  try {
    mesh_etag.setData(faces.get(ARRAY_SIZE(faces)-1),
                      double_tag, dum_val, dum_val_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to set double tag (val=2.0e9) in entityset_double_tag_test."
	      << std::endl;
    return false;
  }

  dum_val_ptr[0] = 0.0;
  try {
    mesh_etag.getData(faces.get(ARRAY_SIZE(faces)-1), double_tag,
                      dum_val, dum_val_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to get double tag (val=2.0e9) in entityset_double_tag_test."
	      << std::endl;
    return false;
  }

  if (dum_val_ptr[0] != 2.0e9) {
    std::cout << "Value of entityset tag (val=2.0e9) wrong." << std::endl;
    return false;
  }

  // ok, the double tag test worked
  return true;
}

bool entityset_struct_tag_test(iMesh::Mesh &mesh, sidl::array<void*> &faces,
			       void *&struct_tag) 
{
  CAST_ITAPS_INTERFACE(mesh, mesh_etag, EntTag);
      
  // Add a struct entityset Tag to the database
  struct SetTagStruct {
    double test_double;
    int test_int1, test_int2;
  };

  // create a tag
  std::string tag_name("set_struct_tag");
  try {mesh_etag.createTag(tag_name, sizeof(SetTagStruct), iBase::TagValueType_BYTES, struct_tag);}
  catch (iBase::Error) {
    std::cout << "Failed to create entityset tag struct_tag in entityest_struct_tag_test." << std::endl;
    return false;
  }

  // put a value in the first set and retrieve
  sidl::array<char> dum_val;
  int dum_val_size = sizeof(SetTagStruct);
  dum_val = dum_val.create1d(dum_val_size);
  SetTagStruct *dum_val_ptr = ARRAY_PTR(dum_val, SetTagStruct);
  SetTagStruct struct_val = {3.0e12, 2, 3};
  dum_val_ptr[0] = struct_val;

  try {
    mesh_etag.setData(faces.get(0), struct_tag, dum_val, dum_val_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to set struct tag in entityset_struct_tag_test."
	      << std::endl;
    return false;
  }

  SetTagStruct dum_struct = {0.0, 0, 0};
  dum_val_ptr[0] = dum_struct;
  try {
    mesh_etag.getData(faces.get(0), struct_tag,
                      dum_val, dum_val_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to get struct tag in entityset_struct_tag_test."
	      << std::endl;
    return false;
  }

  if (dum_val_size != sizeof(SetTagStruct)) {
    std::cout << "Size of entityest struct tag wrong." << std::endl;
    return false;
  }

  if(struct_val.test_int1 != dum_val_ptr[0].test_int1 ||
     struct_val.test_int2 != dum_val_ptr[0].test_int2 ||
     struct_val.test_double != dum_val_ptr[0].test_double) {
    std::cout << "Value of entityset struct tag wrong." << std::endl;
    return false;
  }
  
  // ok, the struct tag test worked
  return true;
}

bool entityset_tag_delete_test(iMesh::Mesh &mesh, sidl::array<void*> &faces) 
{
    // test forced, unforced deletion of tags from entityset
  CAST_ITAPS_INTERFACE(mesh, mesh_etag, EntTag);
  CAST_ITAPS_INTERFACE(mesh, mesh_atag, ArrTag);

    // test getAllTagHandles for first entity
  sidl::array<void*> all_tags, dum_entities;
  int all_tags_size, dum_entities_size;
  try {
    mesh_etag.getAllTags(faces.get(0), all_tags, all_tags_size);
  }
  catch (iBase::Error) {
    std::cerr << "Couldn't get all tag handles from entityset." << std::endl;
    return false;
  }

    // try removing one
  if (all_tags_size > 0) {
    dum_entities = dum_entities.create1d(1);
    dum_entities_size = 1;
    dum_entities.set(0, faces.get(0));
    try {
      mesh_atag.rmvArrTag(dum_entities, dum_entities_size, all_tags.get(0));
    }
    catch (iBase::Error) {
      std::cerr << "Couldn't remove tag from entityset." << std::endl;
      return false;
    }
  
    // try deleting another unforced - should not get error
    if (all_tags_size > 1) {
      bool delete_err = false;
      try {mesh_atag.destroyTag(all_tags.get(1), false);}
      catch (iBase::Error) {
        delete_err = true;
      }
      if (!delete_err) {
        std::cerr << "Error when unforced-deleting tag in use on an entityset." << std::endl;
        return false;
      }

      // now try doing a forced delete
      try {mesh_atag.destroyTag(all_tags.get(1), true);}
      catch (iBase::Error) {
        std::cerr << "Couldn't force-delete a tag in use on an entityset." << std::endl;
        return false;
      }
    }
  }

    // ok, we're done
  return true;
}

/*!
@test
ITAPS entityset tag Test
@li Check entityset tags
*/
bool entityset_tag_test(iMesh::Mesh &mesh)
{
  CAST_ITAPS_INTERFACE(mesh, mesh_stag, SetTag);
  CAST_ITAPS_INTERFACE(mesh, mesh_eset, EntSet);

  // get the sets
  sidl::array<void*> esets;
  int esets_size;
  try {
    mesh_eset.getEntSets (0, 1, esets, esets_size);
  }
  catch (iBase::Error) {
    std::cout << "entitysetGetEntities failed in entityset_tag_test." 
	      << std::endl;
    return false;
  }

  void *int_tag, *double_tag, *struct_tag;
  
    // entityset int tag
  bool int_success = entityset_int_tag_test(mesh, esets, int_tag);

    // entityeset double tag
  bool double_success = entityset_double_tag_test(mesh, esets, double_tag);

    // entityset struct tag
  bool struct_success = entityset_struct_tag_test(mesh, esets, struct_tag);

  bool tag_delete_success = entityset_tag_delete_test(mesh, esets);

  return (int_success && double_success && struct_success
	  && tag_delete_success);
}

bool mesh_int_tag_test(iMesh::Mesh &mesh, void *&int_tag) 
{
  CAST_ITAPS_INTERFACE(mesh, mesh_stag, SetTag);
  
  // create a tag
  std::string tag_name("mesh_int_tag");
  try {
    mesh_stag.createTag(tag_name, 1, iBase::TagValueType_INTEGER, int_tag);
  }
  catch (iBase::Error) {
    std::cout << "Failed to create tag int_tag in mesh_int_tag_test." << std::endl;
    return false;
  }

  // put a value in the first set and retrieve
  int dum_val_size = sizeof(int);
  sidl::array<char> dum_val = dum_val.create1d(dum_val_size);
  ARRAY_PTR(dum_val, int)[0] = 11;

  try {
    mesh_stag.setEntSetData(0, int_tag, dum_val, dum_val_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to set int tag (val=11) in mesh_int_tag_test." << std::endl;
    return false;
  }

  ARRAY_PTR(dum_val, int)[0] = 0;
  try {
    mesh_stag.getEntSetData(0, int_tag,
                            dum_val, dum_val_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to get int tag (val=11) in mesh_int_tag_test."
	      << std::endl;
    return false;
  }

  if (dum_val_size != sizeof(int) || ARRAY_SIZE(dum_val) != sizeof(int) ||
      ARRAY_PTR(dum_val, int)[0] != 11) {
    std::cout << "Value of mesh tag (val=11) wrong." << std::endl;
    cout << "dum_val=" << dum_val << ",dum_val_ptr=" << ARRAY_PTR(dum_val, int) << endl;
    return false;
  }

  // ok, the int tag test worked
  return true;
}

bool mesh_double_tag_test(iMesh::Mesh &mesh, void *&double_tag) 
{
  CAST_ITAPS_INTERFACE(mesh, mesh_stag, SetTag);

  // create a tag
  std::string tag_name("mesh_double_tag");
  try {mesh_stag.createTag(tag_name, 1, iBase::TagValueType_DOUBLE, double_tag);}
  catch (iBase::Error) {
    std::cout << "Failed to create tag mesh_double_tag in mesh_double_tag_test." << std::endl;
    return false;
  }

  int dum_val_size = sizeof(double);
  sidl::array<char> dum_val = dum_val.create1d(dum_val_size);
  ARRAY_PTR(dum_val, double)[0] = 1.0e6;

  try {
    mesh_stag.setEntSetData(0, double_tag, dum_val, dum_val_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to mesh double tag (val=1.0e6) in mesh_double_tag_test."
	      << std::endl;
    return false;
  }

  ARRAY_PTR(dum_val, double)[0] = 0.0;
  try {
    mesh_stag.getEntSetData(0, double_tag,
                            dum_val, dum_val_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to get double tag (val=1.0e6) in mesh_double_tag_test."
	      << std::endl;
    return false;
  }

  if ((double) ARRAY_PTR(dum_val, double)[0] != 1.0e6) {
    std::cout << "Value of mesh tag (val=1.0e6) wrong." 
	      << std::endl;
    return false;
  }

  // ok, the double tag test worked
  return true;
}

bool mesh_struct_tag_test(iMesh::Mesh &mesh, void *&struct_tag) 
{
  CAST_ITAPS_INTERFACE(mesh, mesh_stag, SetTag);
      
  // Add a struct entityset Tag to the database
  struct SetTagStruct {
    double test_double;
    int test_int1, test_int2;
  };

  // create a tag
  std::string tag_name("mesh_struct_tag");
  try {mesh_stag.createTag(tag_name, sizeof(SetTagStruct), iBase::TagValueType_BYTES, struct_tag);}
  catch (iBase::Error) {
    std::cout << "Failed to create mesh tag struct_tag in mesh_struct_tag_test." << std::endl;
    return false;
  }

  // put a value in the first set and retrieve
  int dum_val_size = sizeof(SetTagStruct);
  sidl::array<char> dum_val = dum_val.create1d(dum_val_size);
  SetTagStruct *dum_val_ptr = ARRAY_PTR(dum_val, SetTagStruct);
  SetTagStruct struct_val = {3.0e12, 2, 3};
  dum_val_ptr[0].test_double = struct_val.test_double;
  dum_val_ptr[0].test_int1 = struct_val.test_int1;
  dum_val_ptr[0].test_int2 = struct_val.test_int2;
  
  try {
    mesh_stag.setEntSetData(0, struct_tag, dum_val, dum_val_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to set struct tag in mesh_struct_tag_test."
	      << std::endl;
    return false;
  }

  SetTagStruct dum_struct = {0.0, 0, 0};
  ARRAY_PTR(dum_val, SetTagStruct)[0] = dum_struct;
  try {
    mesh_stag.getEntSetData(0, struct_tag,
                            dum_val, dum_val_size);
  }
  catch (iBase::Error) {
    std::cout << "Failed to get struct tag in mesh_struct_tag_test."
	      << std::endl;
    return false;
  }

  if (dum_val_size != sizeof(SetTagStruct)) {
    std::cout << "Size of mesh struct tag wrong." << std::endl;
    return false;
  }

  dum_val_ptr = ARRAY_PTR(dum_val, SetTagStruct);
  if(struct_val.test_int1 != dum_val_ptr[0].test_int1 ||
     struct_val.test_int2 != dum_val_ptr[0].test_int2 ||
     struct_val.test_double != dum_val_ptr[0].test_double) {
    std::cout << "Value of mesh struct tag wrong." << std::endl;
    return false;
  }
  
  // ok, the struct tag test worked
  return true;
}

bool mesh_tag_delete_test(iMesh::Mesh &mesh) 
{
    // test forced, unforced deletion of tags from mesh
  CAST_ITAPS_INTERFACE(mesh, mesh_stag, SetTag);

    // test getAllTagHandles of mesh
  sidl::array<void*> all_tags;
  int all_tags_size;
  try {
    mesh_stag.getAllEntSetTags(0, all_tags, all_tags_size);
  }
  catch (iBase::Error) {
    std::cerr << "Couldn't get all tag handles from mesh." << std::endl;
    return false;
  }

    // try removing one
  if (all_tags_size > 0) {
    try {
      mesh_stag.rmvEntSetTag(0, all_tags.get(0));
    }
    catch (iBase::Error) {
      std::cerr << "Couldn't remove tag from mesh." << std::endl;
      return false;
    }
  
    // try deleting another unforced - should not get error
    if (all_tags_size > 1) {
      bool delete_err = false;
      try {mesh_stag.destroyTag(all_tags.get(1), false);}
      catch (iBase::Error) {
        delete_err = true;
      }
      if (!delete_err) {
        std::cerr << "Error when unforced-deleting tag in use on an mesh." << std::endl;
        return false;
      }

      // now try doing a forced delete
      try {mesh_stag.destroyTag(all_tags.get(1), true);}
      catch (iBase::Error) {
        std::cerr << "Couldn't force-delete a tag in use on an mesh." << std::endl;
        return false;
      }
    }
  }

    // ok, we're done
  return true;
}

/*!
@test
ITAPS mesh tag Test
@li Check mesh tags
*/
bool mesh_tag_test(iMesh::Mesh &mesh)
{
  void *int_tag, *double_tag, *struct_tag;
  
  // mesh int tag
  bool int_success = mesh_int_tag_test(mesh, int_tag);
    
  // mesh double tag
  bool double_success = mesh_double_tag_test(mesh, double_tag);
    
  // mesh struct tag
  bool struct_success = mesh_struct_tag_test(mesh, struct_tag);
  
  bool tag_delete_success = mesh_tag_delete_test(mesh);
  
  return (int_success && double_success && struct_success
	  && tag_delete_success);
}
  
int main( int argc, char *argv[] )
{
  // Check command line arg
  std::string filename;

  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " <mesh_filename>" << std::endl;
    return 1;
  }
  filename = argv[1];

  bool result;
  int number_tests = 0;
  int number_tests_successful = 0;
  int number_tests_not_implemented = 0;
  int number_tests_failed = 0;

  // initialize the Mesh
  iMesh::Mesh mesh;
  std::string options;
  try {
//    iMesh_SIDL::MeshSidl::newMesh(options, mesh);
    mesh = iMesh_SIDL::MeshSidl::_create();
  }
  catch (iBase::Error err) {
    cout << "Trouble creating a mesh instance." << endl;
    return 1;
  }

  // Print out Header information
  cout << "\n\nITAPS TEST PROGRAM:\n\n";

  // load_mesh test
  cout << "   load_mesh: ";
  result = load_mesh_test(filename, mesh);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  cout << "\n";

  // topology_adjacency_test
  cout << "   topology_adjacency_test: ";
  result = topology_adjacency_test(mesh);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  cout << "\n";

  // entity connectivity test
  cout << "   entity_connectivity_test: ";
  result = entity_connectivity_test(mesh);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  cout << "\n";
    
  // vertex_coordinates_test
  cout << "   vertex_coordinates_test: ";
  result = vertex_coordinates_test(mesh);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  cout << "\n";
  
  // topology dimension test
  cout << "   topology_dimension_test: ";
  result = topology_dimension_test(mesh);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  cout << "\n";

  // entity sets test
  cout << "   entity_sets_test: ";
  result = entity_sets_test(mesh);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  cout << "\n";

  // vertex tag test
  cout << "   vertex_tag_test: ";
  result = vertex_tag_test(mesh);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  cout << "\n";

  // entityset tag test
  cout << "   entityset_tag_test: ";
  result = entityset_tag_test(mesh);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  cout << "\n";

  // mesh tag test
  cout << "   mesh_tag_test: ";
  result = mesh_tag_test(mesh);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  cout << "\n";

  // summary

  cout << "\nITAPS TEST SUMMARY: \n"
       << "   Number Tests:           " << number_tests << "\n"
       << "   Number Successful:      " << number_tests_successful << "\n"
       << "   Number Not Implemented: " << number_tests_not_implemented << "\n"
       << "   Number Failed:          " << number_tests_failed 
       << "\n\n" << endl;
  
  return 0;
}








