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

/**
 * \class  ReadHDF5
 * \brief  Read mesh from TSTT HDF5 file.
 * \author Jason Kraftcheck
 * \date   18 April 2004
 */

#ifndef READ_HDF5_HPP
#define READ_HDF5_HPP

#include <stdlib.h>
#include <list>
#include "mhdf.h"
#include "MBInterface.hpp"
#include "MBReadUtilIface.hpp"
#include "MBRange.hpp"
#include "MBReaderIface.hpp"

class MB_DLL_EXPORT ReadHDF5 : public MBReaderIface
{
public:

  static MBReaderIface* factory( MBInterface* );

  ReadHDF5( MBInterface* iface );
  
  virtual ~ReadHDF5();
  
  /** Export specified meshsets to file
   * \param filename     The filename to export.  Must end in <em>.mhdf</em>
   * \param export_sets  Array of handles to sets to export, or NULL to export all.
   * \param export_set_count Length of <code>export_sets</code> array.
   */
  MBErrorCode load_file( const char* filename,
                         const int* material_set_list,
                         int material_set_count  );

private:
  MBErrorCode init();
  
  //! Range of entities, grouped by type 
  struct ElemSet 
  {
    //! The range of entities.
    MBRange range;
    //! The type of the entities in the range
    MBEntityType type;
    //! The type handle for the mhdf library.
    const char* type2;
    //! The first Id allocated by the mhdf library.  Entities in range have sequential IDs.
    MBEntityHandle first_id;
  };
  
  //! The size of the data buffer (<code>dataBuffer</code>).
  const int bufferSize;
  //! A memory buffer to use for all I/O operations.
  char* dataBuffer;

  //! MBInterface pointer passed to constructor
  MBInterface* iFace;
  
  //! The file handle from the mhdf library
  mhdf_FileHandle filePtr;
  
  //! Cache pointer to read util
  MBReadUtilIface* readUtil;
  
  //! The list elements to export.
  std::list<ElemSet> elemList;
  
  //! The list of nodes to export
  ElemSet nodeSet;
  
  //! The list of sets to export
  ElemSet setSet;
  
  //! The type of an MBEntityHandle
  hid_t handleType;

  //! Get entity given ID
  MBEntityHandle handle_from_id( id_t id );
  
  //! Get range of entities given range of IDs
  MBErrorCode handles_from_ids( id_t start_id, id_t count, MBRange& merge );
  
  //! Get list of handles from list of IDs
  MBErrorCode handles_from_ids( id_t* id_array, id_t count, 
                                std::vector<MBEntityHandle>& append );
 
  //! Read node coordinates.
  MBErrorCode read_nodes( );
  
  //! Read element connectivity.
  MBErrorCode read_elems( const char* elem_group );
  
  //! Read poly(gons|hedra)
  MBErrorCode read_poly( const char* elem_group );
  
  //! Read sets
  MBErrorCode read_sets( );
  
  //! Read element adjacencies
  MBErrorCode read_adjacencies( ElemSet& for_this_set );
  
  //! Create tag and read all data.
  MBErrorCode read_tag( const char* name );
  
  //! Read dense tag for all entities 
  MBErrorCode read_dense_tag( MBTag tag_handle,
                              hid_t hdf_read_type,
                              size_t read_size,
                              bool is_handle_type );

  //! Read dense tag for specified entity set
  MBErrorCode read_dense_tag( ElemSet& set,
                              MBTag tag_handle,
                              hid_t hdf_read_type,
                              size_t read_size,
                              bool is_handle_type );
  
  //! Read sparse tag for all entities.
  MBErrorCode read_sparse_tag( MBTag tag_handle,
                               hid_t hdf_read_type,
                               size_t read_size,
                               bool is_handle_type );
                               
  MBErrorCode read_qa( MBEntityHandle& import_set_out );
                               
  MBErrorCode convert_id_to_handle( const ElemSet& elems_in_this_set,
                                    MBEntityHandle* in_out_array,
                                    size_t array_length );
                               
  MBErrorCode convert_id_to_handle( MBEntityHandle* in_out_array,
                                    size_t array_length );
                                    
  MBErrorCode convert_range_to_handle( const MBEntityHandle* ranges,
                                       size_t num_ranges,
                                       MBRange& merge );
};

#endif
