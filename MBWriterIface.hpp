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
 *\class MBWriterIface
 *\brief Interface for mesh writer implementations.
 *\version 1.00
 *\date 2004-4-23
 *\author Jason Kraftcheck
 */

#ifndef MB_WRITER_IFACE_HPP
#define MB_WRITER_IFACE_HPP

#include <vector>
#include <string>
#include "MBInterface.hpp"

class MBWriterIface
{
  public:
  
    virtual ~MBWriterIface() {}
    
    /**
     *\brief Export mesh to a file.
     *
     * Method all writers must provide to export a mesh.
     *
     *\param file_name      The name of the file to create.
     *\param overwrite      If false, reader should fail if the file already
     *                      exists.  
     *\param meshset_list   A list of meshsets to export, or NULL if the
     *                      whole mesh is to be exported.
     *\param num_sets       The length of <code>meshset_list</code> or zero
     *                      if the whole mesh is to be exported.
     *\param qa_records     File history metadata
     *\param requseted_output_dimension  The geometric dimension of the
     *                      output mesh (coord values per vertex.)  If
     *                      zero, the dimension of the mesh as returned
     *                      from MBInterface should be used.
     *\author Jason Kraftcheck
     */
    virtual MBErrorCode write_file( const char* file_name,
                                    const bool overwrite,
                                    const MBEntityHandle* meshset_list,
                                    const int num_sets,
                                    std::vector<std::string>& qa_records,
                                    int requested_output_dimension = 0 ) = 0;
};

#endif

    
    
