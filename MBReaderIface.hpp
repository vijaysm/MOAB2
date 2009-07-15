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
    
    /**
     *\brief Load mesh from a file.
     *
     * Method all readers must provide to import a mesh.
     *
     *\param file_name      The file to read.
     *\param file_set       Output: a new entity set containing all data read from file.
     *\param set_tag_name   If only reading part of the file, the entities
     *                      to be read will be identified by their values
     *                      for this integer tag.
     *\param set_tag_values For the integer tag with the name indicated by
     *                      set_tag_name, the list of tag values for entities/sets
     *                      to read.
     *\param num_set_tag_values The length of the 'set_tag_values' array.
     *\param file_id_tag    If specified, reader should store for each entity
     *                      it reads, a unique integer ID for this tag.
     *\author Jason Kraftcheck
     */
    virtual MBErrorCode load_file( const char* file_name,
                                   MBEntityHandle& file_set,
                                   const FileOptions& opts,
                                   const char* set_tag_name = 0,
                                   const int* set_tag_values = 0,
                                   int num_set_tag_values = 0,
                                   const MBTag* file_id_tag = 0 ) = 0;

};

#endif

    
