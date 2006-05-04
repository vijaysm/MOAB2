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


#ifndef MB_READ_UTIL_HPP
#define MB_READ_UTIL_HPP

#ifndef IS_BUILDING_MB
#error "MBReadUtil.hpp isn't supposed to be included into an application"
#endif

#include "MBReadUtilIface.hpp"

class MBCore;
class MBError;

class MBReadUtil : public MBReadUtilIface
{
private:
  //! pointer to the MBCore
  MBCore* mMB;
  MBError* mError;
public:

  //! constructor takes MBCore pointer
  MBReadUtil(MBCore* mdb, MBError* error_handler);

  //! destructor
  ~MBReadUtil(){}

  //! gets arrays for coordinate data from the MB
  MBErrorCode get_node_arrays(
      const int num_arrays,
      const int num_nodes, 
      const int preferred_start_id,
      const int preferred_start_proc,
      MBEntityHandle& actual_start_handle, 
      std::vector<double*>& arrays
      );

  //! get array for connectivity data from the MB
  MBErrorCode get_element_array(
      const int num_elements, 
      const int verts_per_element,
      const MBEntityType mdb_type,
      int preferred_start_id, 
      int preferred_start_proc, 
      MBEntityHandle& actual_start_handle, 
      MBEntityHandle*& array
      );

  /** Allocate storage for poly (polygon or polyhedron elements) 
   * 
   * Allocate storage for poly (polygon or polyhedron elements) and
   * return connectivity arrays for direct read into memory.
   *
   *\param num_poly            The number of polygons to allocate handles for.
   *\param conn_list_length    The total length of the connectivity list.
   *\param mdb_type            <code>MBPOLYGON</code> or <code>MBPOLYHEDRON</code>
   *\param preferred_start_id  Preferred integer id for first element
   *\param actual_start_handle Actual integer id for first element (returned)
   *\param last_index_array    Array of indices into <code>connectivity_array</code<
   *\param connectivity_array  The connectivity array
   *\author Jason Kraftcheck
   */
  MBErrorCode get_poly_element_array(
      const int num_poly, 
      const int conn_list_length,
      const MBEntityType mdb_type,
      const int preferred_start_id, 
      const int preferred_start_proc, 
      MBEntityHandle& actual_start_handle, 
      int*& last_index_array,
      MBEntityHandle*& connectivity_array
      );
 
  //! tell MB which elements have been added to the database
  MBErrorCode update_adjacencies(
      const MBEntityHandle start_handle,
      const int number_elements,
      const int number_vertices_per_element,
      const MBEntityHandle* conn_array);

  //! tell MB there was an error when reading the mesh
  //! it makes sense to have this as long as MBInterface has a load_mesh function
  MBErrorCode report_error( const std::string& error );

  MBErrorCode report_error( const char* error, ... );

};

#endif


