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

//-------------------------------------------------------------------------
// Filename      : ReadVtk.hpp
//
// Purpose       : Vtk reader
//
// Special Notes : Lots of code taken from verde implementation
//
//-------------------------------------------------------------------------

#ifndef READVTK_HPP
#define READVTK_HPP

#ifndef IS_BUILDING_MB
#error "ReadVtk.hpp isn't supposed to be included into an application"
#endif

#include "MBInterface.hpp"
#include "MBReaderIface.hpp"

class MBReadUtilIface;

class ReadVtk : public MBReaderIface
{
   
public:

  static MBReaderIface* factory( MBInterface* );

    //! load a file
  MBErrorCode load_file(const char *file_name,
                        const int* material_set_list,
                        const int num_material_sets );
  
    //! Constructor
  ReadVtk(MBInterface* impl = NULL);

   //! Destructor
  virtual ~ReadVtk() {}

private:

  MBReadUtilIface* readMeshIface;

  //------------member variables ------------//

    //! interface instance
  MBInterface* mdbImpl;
  
    //! file name
  std::string fileName;

    //! Meshset Handle for the mesh that is currently being read
  MBEntityHandle mCurrentMeshHandle;
};

#endif




