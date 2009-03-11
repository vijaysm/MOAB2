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
// Filename      : ReadCGM.hpp
//
// Purpose       : .sat, .step and .brep file reader
//
// Special Notes : Lots of code taken from cgm2moab implementation
//
// Creator       : Jane Hu
//
// Date          : 3/09
//
//-------------------------------------------------------------------------

#ifndef READCGM_HPP
#define READCGM_HPP

#ifndef CGM
#error "ReadCGM.hpp isn't supposed to be included without building CGM"
#endif

#include <string>
class MBReadUtilIface;
class GeomTopoTool;

#include "MBReaderIface.hpp"
class ReadCGM : public MBReaderIface
{

public:

  static MBReaderIface* factory( MBInterface* );

    //! load an ExoII file
  MBErrorCode load_file(const char *cgm_file_name,
                         MBEntityHandle& file_set,
                         const FileOptions& opts,
                         const int* blocks_to_load,
                         const int num_blocks);

   //! Constructor
   ReadCGM(MBInterface* impl = NULL);

   //! Destructor
  virtual ~ReadCGM();

private:

  MBReadUtilIface* readUtilIface;

  GeomTopoTool* myGeomTool;

  const char* get_geom_file_type( const char* filename );
  const char* get_geom_fptr_type( FILE* file );

  int is_cubit_file( FILE* file );
  int is_step_file( FILE* file );
  int is_iges_file( FILE* file );
  int is_acis_txt_file( FILE* file );
  int is_acis_bin_file( FILE* file );
  int is_occ_brep_file( FILE* file );


  //------------member variables ------------//

    //! interface instance
  MBInterface* mdbImpl;

    //! file name
  std::string cgmFile;

    //! Meshset Handle for the mesh that is currently being read
  MBEntityHandle mCurrentMeshHandle;

  MBTag geom_tag, id_tag, name_tag, category_tag;
};

#endif
