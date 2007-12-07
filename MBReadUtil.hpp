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
  
  //! Get Processor ID
  unsigned parallel_rank() const;

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
      const int preferred_start_id, 
      const int preferred_start_proc, 
      MBEntityHandle& actual_start_handle, 
      MBEntityHandle*& array
      );

    /**
     *\brief Gather entities related to those in the partition
     * Gather entities related to those in the input partition.  Related
     * means down-adjacent to, contained in, etc.
     * \param partition Entities for which to gather related entities
     * \param related_ents Related entities
     * \param all_sets If non-NULL, all sets in mesh instance are returned
     * in the pointed-to range
     */
  MBErrorCode gather_related_ents(MBRange &partition,
                                  MBRange &related_ents,
                                  MBRange *all_sets);
  
  MBErrorCode create_entity_sets(
    MBEntityID num_sets,
    const unsigned* set_flags,
    MBEntityID preffered_start_id,
    int preffered_start_proc,
    MBEntityHandle& actual_start_handle
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

  MBErrorCode report_error( const char* error, ... )
#ifdef __GNUC__
__attribute__((format(printf,2,3)))
#endif
  ;
};

#endif


