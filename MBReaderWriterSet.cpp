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

#include "MBCore.hpp"
#include "MBError.hpp"

#include "MBReaderWriterSet.hpp"
#include "MBReaderIface.hpp"
#include "MBWriterIface.hpp"

#include "ReadVtk.hpp"
#include "ReadSTL.hpp"
#include "ReadGmsh.hpp"
#include "Tqdcfr.hpp"

#include "WriteAns.hpp"
#include "WriteVtk.hpp"
#include "WriteGMV.hpp"
#include "WriteSTL.hpp"
#include "WriteGmsh.hpp"

#ifdef NETCDF_FILE
#  include "ReadNCDF.hpp"
#  include "WriteNCDF.hpp"
#  include "WriteSLAC.hpp"
#endif

#ifdef HDF5_FILE
#  include "ReadHDF5.hpp"
#  include "WriteHDF5.hpp"
#endif

MBReaderWriterSet::MBReaderWriterSet( MBCore* mdb, MBError* handler )
  : mbCore( mdb ), mbError( handler ) 
{
#ifdef HDF5_FILE
  const char* hdf5_list[] = { "h5m", "mhdf", NULL };
  register_factory(  ReadHDF5::factory, WriteHDF5::factory, "MOAB HDF5", hdf5_list );
#endif

#ifdef NETCDF_FILE
  const char* exo_list[] = { "exo", "exoII", "exo2", "g", "gen", NULL };
  register_factory( ReadNCDF::factory, WriteNCDF::factory, "Exodus II", exo_list );
#endif
  
  const char* vtk_list[] = { "vtk", NULL };
  register_factory( ReadVtk::factory, WriteVtk::factory, "Kitware VTK", vtk_list );
  
  const char* cub_list[] = { "cub", NULL };
  register_factory( Tqdcfr::factory, NULL, "Cubit", cub_list );

#ifdef NETCDF_FILE  
  const char* slac_list[] = { "slac", NULL };
  register_factory( NULL, WriteSLAC::factory, "SLAC", slac_list );
#endif

  const char* gmv_list[] = { "gmv", NULL };
  register_factory( NULL, WriteGMV::factory, "GMV", gmv_list );
  
  const char* ans_list[] = { "ans", NULL };
  register_factory( NULL, WriteAns::factory, "Ansys", ans_list );
  
  const char* gmsh_list[] = { "msh", "gmsh", NULL };
  register_factory( ReadGmsh::factory, WriteGmsh::factory, "Gmsh mesh file", gmsh_list );
  
  const char* stl_list[] = { "stl", NULL };
  register_factory( ReadSTL::ascii_instance, WriteSTL::ascii_instance, "Stereo Lithography File (STL)", stl_list );
  
  const char* stlb_list[] = { "stlb", NULL };
  register_factory( ReadSTL::binary_instance, WriteSTL::binary_instance, "Binary Stereo Lithography (STL)", stlb_list );
}


MBReaderWriterSet::~MBReaderWriterSet()
{
}

MBErrorCode MBReaderWriterSet::register_factory( reader_factory_t reader,
                                                 writer_factory_t writer,
                                                 const char* description,
                                                 const char** extensions )
{
  if (!reader && !writer)
    return MB_FAILURE;
  
    // count extensions and check for duplicates
  const char** iter;
  for (iter = extensions; *iter; ++iter)
  {
    iter_type h = handler_from_extension( *iter );
    if (h != end())
    {
      if (NULL != reader && h->have_reader())
        mbError->set_last_error( "Conflicting readers for file extension \"%s\":"
                                 " \"%s\" and \"%s\".",
                                 *iter, h->description().c_str(), description );
      else if(NULL != writer && h->have_writer())
        mbError->set_last_error( "Conflicting writers for file extension \"%s\":"
                                 " \"%s\" and \"%s\".",
                                 *iter, h->description().c_str(), description );
    }
  }
  handlerList.push_back( Handler(reader, writer, description, extensions, iter - extensions) );
  return MB_SUCCESS;
}    
  
MBReaderIface* MBReaderWriterSet::get_file_extension_reader( 
                                  const std::string& filename ) const
{
  std::string ext = extension_from_filename( filename );
  iter_type handler = handler_from_extension( ext, true, false );
  return handler == end() ? NULL : handler->make_reader(mbCore);
}

MBWriterIface* MBReaderWriterSet::get_file_extension_writer( 
                                  const std::string& filename ) const
{
  std::string ext = extension_from_filename( filename );
  iter_type handler = handler_from_extension( ext, false, true );
  return handler == end() ? NULL : handler->make_writer(mbCore);
}

std::string MBReaderWriterSet::extension_from_filename( 
                                 const std::string& filename )
{
  std::string::size_type idx = filename.find_last_of( "." );
  if (idx == std::string::npos)
    return std::string("");
  else
    return filename.substr( idx + 1 );
}

MBReaderWriterSet::Handler::Handler( reader_factory_t read_f, 
                                     writer_factory_t write_f,
                                     const char* desc, 
                                     const char** ext, 
                                     int num_ext )
 : mReader(read_f), mWriter(write_f), mDescription(desc), mExtensions(num_ext)
{
  for (int i = 0; i < num_ext; ++i)
    mExtensions[i] = ext[i];
}

#ifdef WIN32
#define strcasecmp(A,B) _stricmp( A, B )
#endif

MBReaderWriterSet::iter_type 
MBReaderWriterSet::handler_from_extension( const std::string& ext,
                                           bool with_reader,
                                           bool with_writer ) const
{
  iter_type iter;
  std::vector<std::string>::const_iterator siter;
  
    // try case-sensitive compare
  for (iter = begin(); iter != end(); ++iter)
  {
    if ((with_reader && !iter->have_reader()) ||
        (with_writer && !iter->have_writer()))
      continue;
      
    for (siter = iter->mExtensions.begin(); siter != iter->mExtensions.end(); ++siter)
      if (*siter == ext)
        return iter;
  }
  
    // try case-insensitive compare
  for (iter = begin(); iter != end(); ++iter)
  {
    if ((with_reader && iter->have_reader()) ||
        (with_writer && iter->have_writer()))
      continue;
 
    for (siter = iter->mExtensions.begin(); siter != iter->mExtensions.end(); ++siter)
      if (0 == strcasecmp( siter->c_str(), ext.c_str() ))
        return iter;
  }
  
  return end();
}
                                       
