#include "MBCore.hpp"
#include "MBError.hpp"

#include "MBReaderWriterSet.hpp"
#include "MBReaderIface.hpp"
#include "MBWriterIface.hpp"

#include "ReadNCDF.hpp"
#include "ReadVtk.hpp"
#include "Tqdcfr.hpp"

#include "WriteVtk.hpp"
#include "WriteGMV.hpp"
#include "WriteNCDF.hpp"
#include "WriteSLAC.hpp"

#ifdef HDF5_FILE
#  include "ReadHDF5.hpp"
#  include "WriteHDF5.hpp"
#endif

MBReaderWriterSet::MBReaderWriterSet( MBCore* mdb, MBError* handler )
  : mbCore( mdb ), mbError( handler ) 
{
  const char* exo_list[] = { "exo", "exoII", "exo2", "g", "gen", NULL };
  register_factory( ReadNCDF::factory, WriteNCDF::factory, "Exodus II", exo_list );
  
  const char* vtk_list[] = { "vtk", NULL };
  register_factory( ReadVtk::factory, WriteVtk::factory, "Kitware VTK", vtk_list );
  
  const char* cub_list[] = { "cub", NULL };
  register_factory( Tqdcfr::factory, NULL, "Cubit", cub_list );
  
  const char* slac_list[] = { "slac", NULL };
  register_factory( NULL, WriteSLAC::factory, "SLAC", slac_list );
  
  const char* gmv_list[] = { "gmv", NULL };
  register_factory( NULL, WriteGMV::factory, "GMV", gmv_list );
  
#ifdef HDF5_FILE
  const char* hdf5_list[] = { "h5m", "mhdf", NULL };
  register_factory(  ReadHDF5::factory, WriteHDF5::factory, "TSTT HDF5", hdf5_list );
#endif
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
                                       
