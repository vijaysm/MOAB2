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
                                   MBEntityHandle              file_set, 
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

  // create tags
  MBTag id_tag, material_tag;
  result = MBI->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE, 
                           MB_TYPE_INTEGER, id_tag, 0);
  if(MB_SUCCESS!=result && MB_ALREADY_ALLOCATED!=result) return result;
  result = MBI->tag_create(MATERIAL_SET_TAG_NAME, sizeof(int), MB_TAG_DENSE, 
                           MB_TYPE_INTEGER, material_tag, 0);
  if(MB_SUCCESS!=result && MB_ALREADY_ALLOCATED!=result) return result;

  // Count the entities of each type in the file. This is used to allocate the node array. 
  int entity_count[MBMAXTYPE];
  for(int i=0; i<MBMAXTYPE; i++) entity_count[i] = 0;
 
  /* Determine the line_format of the first line. Assume that the entire file 
     has the same format. */  
  std::string line;
  std::ifstream file (filename);
  getline(file,line);
  line_format format;
  result = determine_line_format( line, format );
  if(MB_SUCCESS != result) return result;

  /* Count the number of each entity in the file. This allows one to allocate
     a sequential array of vertex handles. */
  while(!file.eof()) {
    // Cut the line into fields as determined by the line format.
    // Use a vector to allow for an unknown number of tokens (continue lines).
    // Continue lines are not implemented.
    std::vector<std::string> tokens;
    tokens.reserve(10); // assume 10 fields to avoid extra vector resizing
    result = tokenize_line( line, format, tokens );
    if(MB_SUCCESS != result) return result;

    // Process the tokens of the line. The first token describes the entity type.
    MBEntityType type;
    result = determine_entity_type( tokens.front(), type );
    if(MB_SUCCESS != result) return result;
    entity_count[type]++;
    getline(file,line);
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

  // Read the file again to create the entities.
  file.clear();  // clear eof state from object
  file.seekg(0); // rewind file
  while (!file.eof()) {
    getline(file,line);

    // Cut the line into fields as determined by the line format.
    // Use a vector to allow for an unknown number of tokens (continue lines).
    // Continue lines are not implemented.
    std::vector<std::string> tokens;
    tokens.reserve(10); // assume 10 fields to avoid extra vector resizing
    result = tokenize_line( line, format, tokens );
    if(MB_SUCCESS != result) return result;

    // Process the tokens of the line. The first token describes the entity type.
    MBEntityType type;
    result = determine_entity_type( tokens.front(), type );
    if(MB_SUCCESS != result) return result;

    // Create the entity.
    if(MBVERTEX == type) {
      result = read_node(tokens, id_tag, file_set, debug, coord_arrays, 
			 vert_index, start_vert, node_ids_are_continuous ); 
      if(MB_SUCCESS != result) return result;
    } else {
      result = read_element(tokens, id_tag, material_tag, type, file_set, 
			    debug, start_vert, node_ids_are_continuous ); 
      if(MB_SUCCESS != result) return result;
    }
  }
  
  file.close();
  return MB_SUCCESS;
}

  /* Determine the type of NASTRAN line: small field, large field, or free field.
     small field: each line has 10 fields of 8 characters
     large field: 1x8, 4x16, 1x8. Field 1 must have an asterisk following the character string
     free field: each line entry must be separated by a comma
     Implementation tries to avoid more searches than necessary. */
MBErrorCode ReadNASTRAN::determine_line_format( const std::string line, 
                                                  line_format &format ) {
  std::string::size_type found_asterisk = line.find("*");
  if(std::string::npos != found_asterisk) {
      format = LARGE_FIELD;
      return MB_SUCCESS;
    } else {
      std::string::size_type found_comma = line.find(",");
      if(std::string::npos != found_comma) {
	format = FREE_FIELD;
	return MB_SUCCESS;
      } else {
	format = SMALL_FIELD;
	return MB_SUCCESS;
      }
    }
  }

  /* Tokenize the line. Continue-lines have not been implemented. */
MBErrorCode ReadNASTRAN::tokenize_line(const std::string line, const line_format format, 
				       std::vector<std::string> &tokens ) { 
  size_t line_size = line.size();
  switch(format) {
  case SMALL_FIELD: {
    // Expect 10 fields of 8 characters.
    // The sample file does not have all 10 fields in each line
    const int field_length = 8;
    unsigned int n_tokens = line_size / field_length;
    for(unsigned int i=0; i<n_tokens; i++) {
      tokens.push_back( line.substr(i*field_length,field_length) );
    }
    break; 
  } case LARGE_FIELD:
    return MB_NOT_IMPLEMENTED;
    break;
  case FREE_FIELD:
    return MB_NOT_IMPLEMENTED;
    break;
  default:
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}

MBErrorCode ReadNASTRAN::determine_entity_type( const std::string first_token, 
                                                  MBEntityType &type ) {
  if     (0==first_token.compare("GRID    ")) type = MBVERTEX;
  else if(0==first_token.compare("CTETRA  ")) type = MBTET;
  else if(0==first_token.compare("CPENTA  ")) type = MBPRISM;
  else if(0==first_token.compare("CHEXA   ")) type = MBHEX;
  else return MB_NOT_IMPLEMENTED;

  return MB_SUCCESS;
}

/* Some help from Jason:
   Nastran floats must contain a decimal point, may contain
   a leading '-' and may contain an exponent.  The 'E' is optional
   when specifying an exponent.  A '-' or '+' at any location other
   than the first position indicates an exponent.  For a positive
   exponent, either a '+' or an 'E' must be specified.  For a
   negative exponent, the 'E' is option and the '-' is always specified.
   Examples for the real value 7.0 from mcs2006 quick reference guide:
   7.0  .7E1  0.7+1  .70+1  7.E+0  70.-1

   From the test file created in SC/Tetra:
   GRID           1       03.9804546.9052-15.6008-1    
   has the coordinates: ( 3.980454, 6.9052e-1, 5.6008e-1 )
   GRID      200005       04.004752-3.985-15.4955-1  
   has the coordinates: ( 4.004752, -3.985e-1, 5.4955e-1 ) */
MBErrorCode ReadNASTRAN::get_real( const std::string token, double &real ) {
  std::string significand = token;
  std::string exponent = "0";

  // Cut off the first digit because a "-" could be here indicating a negative
  // number. Instead we are looking for a negative exponent.
  std::string back_token = token.substr(1);
  
  // A minus that is not the first digit is always a negative exponent
  std::string::size_type found_minus = back_token.find("-");
  if(std::string::npos != found_minus) {
    // separate the significand from the exponent at the "-"
    exponent = token.substr( found_minus+1 );
    significand = token.substr( 0, found_minus+1 );

    // if the significand has an "E", remove it
    if(std::string::npos != significand.find("E")) 
      // Assume the "E" is at the end of the significand.
      significand = significand.substr( 1, significand.size()-2 );

    // If a minus does not exist past the 1st digit, but an "E" or "+" does, then 
    // it is a positive exponent. First look for an "E",
  } else { 
    std::string::size_type found_E = token.find("E");
    if(std::string::npos != found_E) {
      significand = token.substr(0, found_E-1);
      exponent = token.substr(found_E+1);
      // if there is a "+" on the exponent, cut it off
      std::size_t found_plus = exponent.find("+");
      if(std::string::npos != found_plus) {
	exponent = exponent.substr(found_plus+1);
      }
    } else {
      // if there is a "+" on the exponent, cut it off
      std::size_t found_plus = token.find("+");
      if(std::string::npos != found_plus) {
	significand = token.substr(0, found_plus-1);
	exponent = token.substr(found_plus+1);
      }
    }
  }
    
  // now assemble the real number
  double signi = atof(significand.c_str());
  double expon = atof(exponent.c_str());

  if(HUGE_VAL==signi || HUGE_VAL==expon) return MB_FAILURE;
  real = signi * pow( 10, expon );
  return MB_SUCCESS;
}

/* It has been determined that this line is a vertex. Read the rest of
   the line and create the vertex. */
MBErrorCode ReadNASTRAN::read_node( const std::vector<std::string> tokens, 
				      MBTag id_tag,
				      const MBEntityHandle file_set,
				      const bool debug, 
				      std::vector<double*> coord_arrays, 
				      int &vert_index, 
				      const MBEntityHandle start_vert,
				      bool &node_ids_are_continuous ) {
  // read the node's id (unique)
  MBErrorCode result;
  int id = atoi(tokens[1].c_str());

  // read the node's coordinate system number
  // "0" or blank refers to the basic coordinate system.
  int coord_system = atoi(tokens[2].c_str());
  if(0 != coord_system) {
    std::cerr << "ReadNASTRAN: alternative coordinate systems not implemented" 
              << std::endl;
    return MB_NOT_IMPLEMENTED;    
  }

  // read the coordinates
  double coords[3];
  for(unsigned int i=0; i<3; i++) {
    result = get_real( tokens[i+3], coords[i] );
    if(MB_SUCCESS != result) return result;
    if(debug) std::cout << "read_node: coords[" << i << "]=" << coords[i] << std::endl;
  }

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

/* The type of element has already been identified. Read the rest of the
   line and create the element. Assume that all of the nodes have already
   been created. */
MBErrorCode ReadNASTRAN::read_element( const std::vector<std::string> tokens, 
				       MBTag id_tag, 
				       MBTag material_tag, 
				       const MBEntityType element_type,
				       const MBEntityHandle file_set, 
				       const bool debug, 
				       const MBEntityHandle start_vert,
				       const bool node_ids_are_continuous ) {

  // read the element's id (unique) and material set
  MBErrorCode result;
  int id = atoi(tokens[1].c_str());    
  int material = atoi(tokens[2].c_str());

  // the size of the connectivity array depends on the element type
  int n_conn = MBCN::VerticesPerEntity(element_type);
  
  // read the connected node ids from the file
  int conn[n_conn];
  for(int i=0; i<n_conn; i++) {
    conn[i] = atoi(tokens[3+i].c_str());
  }

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
