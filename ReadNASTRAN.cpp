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

#include "ReadNASTRAN.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <assert.h>
#include <cmath>

#include "MBInterface.hpp"
#include "MBReadUtilIface.hpp"
#include "MBInternals.hpp" // for MB_START_ID
#include "MBRange.hpp"
#include "FileOptions.hpp"
#include "FileTokenizer.hpp"
#include "MBTagConventions.hpp"
#include "MBCN.hpp"

MBReaderIface* ReadNASTRAN::factory( MBInterface* iface ) { 
  return new ReadNASTRAN( iface );
}

// constructor
ReadNASTRAN::ReadNASTRAN(MBInterface* impl)
  : MBI(impl), fileIDTag(0) {
    assert(NULL != impl);
    void *ptr = 0;
    MBI->query_interface("MBReadUtilIface", &ptr);
    assert(NULL != ptr);
    readMeshIface = reinterpret_cast<MBReadUtilIface*>(ptr);
}

// destructor
ReadNASTRAN::~ReadNASTRAN() {
  if (readMeshIface) {
    MBI->release_interface("MBReadUtilIface", readMeshIface);
    readMeshIface = 0;
  }
}

MBErrorCode ReadNASTRAN::read_tag_values( const char*        file_name,
                                          const char*        tag_name,
                                          const FileOptions& opts,
                                          std::vector<int>&  tag_values_out,
                                          const IDTag*       subset_list,
                                          int                subset_list_length )
{
  return MB_NOT_IMPLEMENTED;
}

// load the file as called by the MBInterface function
MBErrorCode ReadNASTRAN::load_file(const char                  *filename, 
                                   MBEntityHandle              &file_set, 
                                   const FileOptions           &options,
                                   const MBReaderIface::IDTag  *subset_list,
                                   int                         subset_list_length,
                                   const MBTag*                file_id_tag) {
  // at this time there is no support for reading a subset of the file
  if (subset_list && subset_list_length) {
    readMeshIface->report_error( "Reading subset of files not supported for NASTRAN." );
    return MB_UNSUPPORTED_OPERATION;
  }
  
  nodeId = elemId = 0;
  fileIDTag = file_id_tag;

  bool debug = false;
  if (debug) std::cout << "begin ReadNASTRAN::load_file" << std::endl;
  MBErrorCode result;
  FILE* file = fopen( filename, "r" );

  // create the file set
  result = MBI->create_meshset(MESHSET_SET, file_set);
  if(MB_SUCCESS != result) return result;

  // create tags
  MBTag id_tag, material_tag;
  result = MBI->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE, 
                           MB_TYPE_INTEGER, id_tag, 0);
  if(MB_SUCCESS!=result && MB_ALREADY_ALLOCATED!=result) return result;
  result = MBI->tag_create(MATERIAL_SET_TAG_NAME, sizeof(int), MB_TAG_DENSE, 
                           MB_TYPE_INTEGER, material_tag, 0);
  if(MB_SUCCESS!=result && MB_ALREADY_ALLOCATED!=result) return result;

  // IMPORTANT: these are the same as the MOAB enums. Currently only vertices,
  // tets, and prisms are supported. Implement more as needed.
  const char* const data_type_names[] = { "GRID",     // MBVERTEX
                                          "MBEDGE",
                                          "MBTRI",
                                          "MBQUAD",
                                          "MBPOLYGON",
                                          "CTETRA",   // MBTET
                                          "MBPYRAMID",
                                          "CPENTA",   // MBPRISM
                                          "MBKNIFE",
                                          "MBHEX",
                                          "MBPOLYHEDRON", 
                                          "MBENTITYSET",
                                          "MBMAXTYPE" };

  /* Count the entities of each type in the file. If the token is not
  matched assume that it is something else in the file (not an error).
  During the next read we will check for errors.
  This is used to allocate the node array. Let the tokenizer run out of 
  scope (and thus call fclose). There is no other way to rewind the 
  tokenizer. */
  int entity_count[MBMAXTYPE];
  for(int i=0; i<MBMAXTYPE; i++) entity_count[i] = 0;
  { FileTokenizer tokens( file, readMeshIface );
    while( !tokens.eof() ) {
      int data = tokens.match_token( data_type_names ) -1;
      if(-1 != data) entity_count[data]++;
    }
  }   
  if(debug) {
    for(int i=0; i<MBMAXTYPE; i++) {
      std::cout << "entity_count[" << i << "]=" << entity_count[i] << std::endl;
    }
  }

  /* If we assume that the node ids are continuous, then the node ids can be
  easily mapped to vertex handles. If they are not continuous then the node
  ids themselves will need to be used to create connectivity. It is very 
  slow to find each node by its id because this is a tag call. */
  bool node_ids_are_continuous = true;

  // Now that the number of vertices is known, create the vertices.
  MBEntityHandle start_vert = NULL;
  std::vector<double*> coord_arrays(3);
  result = readMeshIface->get_node_arrays( 3, entity_count[0], MB_START_ID,
                                           start_vert, coord_arrays );
  if(MB_SUCCESS != result) return result;
  if(0 == start_vert) return MB_FAILURE; // check for NULL
  int vert_index = 0;
  if(debug) std::cout << "allocated coord arrays" << std::endl;

  // Read the file again and build the entities. Each line in the file
  // represents one entity.
  file = fopen( filename, "r" );
  FileTokenizer tokens( file, readMeshIface );
  while( !tokens.eof() ) {
    // these have been enumerated to match MBEntityTypes
    int data = tokens.match_token( data_type_names ) -1;
    MBEntityType datatype = (MBEntityType)data;
    if(debug) {
      std::cout << "MBVERTEX=" << MBVERTEX << " datatype=" << datatype << std::endl;
    }

    // the eof was found
    if( tokens.eof() ) {
      break;

    // no matching entity was found
    } else if(-1 == datatype ) {
      return MB_NOT_IMPLEMENTED;

    // a vertex line was found
    } else if(MBVERTEX == datatype) {
      result = read_node(tokens, id_tag, file_set, debug, coord_arrays, 
                         vert_index, start_vert, node_ids_are_continuous ); 
      if(MB_SUCCESS != result) return result;

    // an element line was found
    } else {
      result = read_element(tokens, id_tag, material_tag, datatype, file_set, 
                            debug, start_vert, node_ids_are_continuous ); 
      if(MB_SUCCESS != result) return result;
    }
  } 
  return MB_SUCCESS;
}

/* It has been determined that this line is a vertex. Read the rest of
   the line and create the vertex. */
MBErrorCode ReadNASTRAN::read_node( FileTokenizer &tokens, 
                                    MBTag id_tag,
                                    const MBEntityHandle file_set,
                                    const bool debug, 
                                    std::vector<double*> coord_arrays, 
                                    int &vert_index, 
                                    const MBEntityHandle start_vert,
                                    bool &node_ids_are_continuous ) {
  // read the node's id (unique)
  MBErrorCode result;
  int id;
  if( !tokens.get_integers(1, &id) ) return MB_FAILURE;;

  // read the coordinates
  const char* coords_string = tokens.get_string();
  if (!coords_string) return MB_FAILURE;
  if(debug) std::cout << "coords_string=" << coords_string << std::endl; 
  double coords[3];
  result = parse_coords( coords_string, coords, debug );
  if(MB_SUCCESS != result) return result;

  // create the node
  MBEntityHandle vert = start_vert+vert_index;
  for(int i=0; i<3; i++) coord_arrays[i][vert_index] = coords[i];
  vert_index++;

  // check to see if the node ids are still continuous
  if(id != vert_index) {
    node_ids_are_continuous = false;
    std::cout << "node ids are not continuous" << std::endl;
  }

  // tag with the global id from the NASTRAN file
  result = MBI->tag_set_data( id_tag, &vert, 1, &id );
  if(MB_SUCCESS != result) return result;

  // add to the file set
  result = MBI->add_entities( file_set, &vert, 1 );
  if(MB_SUCCESS != result) return result;

  return MB_SUCCESS;
}

/* The coordinates are difficult to parse. Below are some samples. */
MBErrorCode ReadNASTRAN::parse_coords( const char *coords_string, 
                                       double     coords[],
                                       const bool debug ) {
  /* GRID           1       03.9804546.9052-15.6008-1    
     has the coordinates: ( 3.980454, 6.9052e-1, 5.6008e-1 )
     GRID      200005       04.004752-3.985-15.4955-1  
     has the coordinates: ( 4.004752, -3.985e-1, 5.4955e-1 ) */
  char x[9], y[8], z[8];
  strncpy(x, coords_string,                     9); //will copy first 4 characters
  strncpy(y, coords_string+sizeof(x),           8);
  strncpy(z, coords_string+sizeof(x)+sizeof(y), 8);
  if(debug) std::cout << "x=" << x << " y=" << y << " z=" << z << std::endl;

  // returns memory location of last match if found. returns null if not found.
  char *x_loc, *y_loc, *z_loc;
  x_loc = strrchr(x,'-'); // - x + 1;
  y_loc = strrchr(y,'-'); // - y + 1;
  z_loc = strrchr(z,'-'); // - z + 1;
  if(debug) std::cout << "loc=" << x_loc << " " << y_loc << " " << z_loc << std::endl;

  // copies( destination memory address, source memory address, size to copy)
  // copy only the exponent, if it exists.
  char exp_x[9], exp_y[8], exp_z[8];
  if(0 != x_loc) {
    memcpy(exp_x, x_loc, 9);
    coords[0] = atof(x) * pow(10, atof(exp_x));
  } else {
    coords[0] = atof(x);
  }
  if(0 != y_loc) {
    memcpy(exp_y, y_loc, 8 );
    coords[1] = atof(y) * pow(10, atof(exp_y));
  } else {
    coords[1] = atof(y);
  }
  if(0 != z_loc) {
    memcpy(exp_z, z_loc, 8-(int)(z_loc-z) );
    coords[2] = atof(z)*pow(10, atof(exp_z));;
  } else {
    coords[2] = atof(z);
  }
  //std::cout << "exp=" << exp_x << " " << exp_y << " " << exp_z << std::endl;
  if(debug) std::cout << "coords[]=" << coords[0] << " " << coords[1] 
                      << " " << coords[2] << std::endl;

  return MB_SUCCESS;
}

/* The type of element has already been identified. Read the rest of the
   line and create the element. Assume that all of the nodes have already
   been created. */
MBErrorCode ReadNASTRAN::read_element( FileTokenizer &tokens, 
                                       MBTag id_tag, 
                                       MBTag material_tag, 
                                       const MBEntityType element_type,
                                       const MBEntityHandle file_set, 
                                       const bool debug, 
                                       const MBEntityHandle start_vert,
                                       const bool node_ids_are_continuous ) {

  // read the element's id (unique) and material set
  MBErrorCode result;
  int id;
  bool res = tokens.get_integers(1, &id);
  if(false == res) return MB_FAILURE;
  int material;
  res = tokens.get_integers(1, &material);
  if(false == res) return MB_FAILURE;

  // the size of the connectivity array depends on the element type
  int n_conn = MBCN::VerticesPerEntity(element_type);
  
  // read the connected node ids from the file
  int conn[n_conn];
  res = tokens.get_integers(n_conn, conn);
  if(false == res) return MB_FAILURE;

  /* Get the vertex handles with those node ids. If the nodes do not have
     continuous ids we need to do a very slow tag search for global ids.
     The cannonical numbering between NASTRAN and MOAB is the same for 
     tets and prisms. */
  MBEntityHandle conn_verts[n_conn];
  for(int i=0; i<n_conn; i++) {
    if(node_ids_are_continuous) {  
      conn_verts[i] = start_vert + conn[i] -1;
    } else {
      void *vert_id_val[] = {&conn[i]};
      MBRange vert;
      result = MBI->get_entities_by_type_and_tag( file_set, MBVERTEX, &id_tag, 
                                                  vert_id_val, 1, vert );
      if(MB_SUCCESS != result) return result;
      if(1 != vert.size()) return MB_FAILURE; // only one vertex should have this global id 
      conn_verts[i] = vert.front();
    }
    if(debug) std::cout << "conn[" << i << "]=" << conn[i] << std::endl;
  }

  // Create the element and set the global id from the NASTRAN file
  MBEntityHandle element;
  result = MBI->create_element( element_type, conn_verts, n_conn, element );
  if(MB_SUCCESS != result) return result;
  result = MBI->tag_set_data( id_tag, &element, 1, &id );
  if(MB_SUCCESS != result) return result;

  // Create the material set for the element if it does not already exist.
  void *mat_val[] = {&material};
  MBRange material_set;
  result = MBI->get_entities_by_type_and_tag( file_set, MBENTITYSET, &material_tag, 
                                                mat_val, 1, material_set );
  if(MB_SUCCESS != result) return result;   
  if(0 == material_set.size()) {
    MBEntityHandle new_material_set;
    result = MBI->create_meshset( 0, new_material_set );
    if(MB_SUCCESS != result) return result;
    result = MBI->add_entities( file_set, &new_material_set, 1 );
    result = MBI->tag_set_data( material_tag, &new_material_set, 1, &material );
    if(MB_SUCCESS != result) return result;
    material_set.insert( new_material_set );
  } else if(1 < material_set.size()) {
    return MB_MULTIPLE_ENTITIES_FOUND;
  }
  result = MBI->add_entities( material_set.front(), &element, 1 );
  if(MB_SUCCESS != result) return result;

  return MB_SUCCESS;
}
