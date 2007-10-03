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
     *\param file_name           The file to read.
     *\param file_set            Output: a new entity set containing all data read from file.
     *\param material_set_list   A list of material sets to read, or NULL
     *                           if the entire file is to be read.
     *\param material_set_list_len The length of <code>material_set_list</code>
     *\author Jason Kraftcheck
     */
    virtual MBErrorCode load_file( const char* file_name,
                                   MBEntityHandle& file_set,
                                   const FileOptions& opts,
                                   const int* material_set_list,
                                   const int material_set_list_len ) = 0;

};

#endif

    
