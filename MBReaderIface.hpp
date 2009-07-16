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
 *\class MBReaderIface
 *\brief Interface for mesh reader implementations.
 *\version 1.00
 *\date 2004-4-23
 *\author Jason Kraftcheck
 */

#ifndef MB_READER_IFACE_HPP
#define MB_READER_IFACE_HPP

#include "MBTypes.h"

class FileOptions;

class MBReaderIface
{
  public:
  
    virtual ~MBReaderIface() {}
    
      /** Struct used to specify subset of file to read */
    struct IDTag {
      const char* tag_name;  //!< Name of tag containing integer IDs
      const int* tag_values; //!< Array of integer ID values
      int num_tag_values;    //!< Length of tag_values array
    };
    
    /**
     *\brief Load mesh from a file.
     *
     * Method all readers must provide to import a mesh.
     *
     *\param file_name      The file to read.
     *\param file_set       Output: a new entity set containing all data read from file.
     *\param subset_list    An array of tag name and value sets specifying
     *                      the subset of the file to read.  If multiple
     *                      tags are specified, the sets that match all
     *                      tags (intersection) should be read.
     *\param subset_list_length The length of the 'subset_list' array.
     *\param file_id_tag    If specified, reader should store for each entity
     *                      it reads, a unique integer ID for this tag.
     *\author Jason Kraftcheck
     */
    virtual MBErrorCode load_file( const char* file_name,
                                   MBEntityHandle& file_set,
                                   const FileOptions& opts,
                                   const IDTag* subset_list = 0,
                                   int subset_list_length = 0,
                                   const MBTag* file_id_tag = 0 ) = 0;


    /**
     *\brief Read tag values from a file.
     *
     * Read the list if all integer tag values from the file for
     * a tag that is a single integer value per entity.
     *
     *\param file_name      The file to read.
     *\param tag_name       The tag for which to read values
     *\param tag_values_out Output: The list of tag values.
     *\param subset_list    An array of tag name and value sets specifying
     *                      the subset of the file to read.  If multiple
     *                      tags are specified, the sets that match all
     *                      tags (intersection) should be read.
     *\param subset_list_length The length of the 'subset_list' array.
     */
    virtual MBErrorCode read_tag_values( const char* file_name,
                                         const char* tag_name,
                                         const FileOptions& opts,
                                         std::vector<int>& tag_values_out,
                                         const IDTag* subset_list = 0,
                                         int subset_list_length = 0 ) = 0;
};

#endif

    
