#include "MBCore.hpp"
#include "MBError.hpp"

#include "MBReaderWriterSet.hpp"
#include "MBReaderIface.hpp"
#include "MBWriterIface.hpp"

#include "ReadNCDF.hpp"
#include "ReadVtk.hpp"
#include "Tqdcfr.hpp"

#include "WriteGMV.hpp"
#include "WriteNCDF.hpp"
#include "WriteSLAC.hpp"

MBReaderWriterSet::MBReaderWriterSet( MBCore* mdb, MBError* handler )
  : mbCore( mdb ), mbError( handler ) 
{
  const char* exo_list[] = { "exo", "exoII", "exo2", "g", "gen", NULL };
  register_reader( ReadNCDF::factory, "Exodus II", exo_list );
  register_writer( WriteNCDF::factory, "Exodus II", exo_list );
  
  const char* vtk_list[] = { "vtk", NULL };
  register_reader( ReadVtk::factory, "Kitware VTK", vtk_list );
  
  const char* cub_list[] = { "cub", NULL };
  register_reader( Tqdcfr::factory, "Cubit", cub_list );
  
  const char* slac_list[] = { "slac", NULL };
  register_writer( WriteSLAC::factory, "SLAC", slac_list );
  
  const char* gmv_list[] = { "gmv", NULL };
  register_writer( WriteGMV::factory, "GMV", gmv_list );
}


MBReaderWriterSet::~MBReaderWriterSet()
{
}

MBErrorCode MBReaderWriterSet::register_reader( reader_factory_t function,
                                                const char* description,
                                                const char** extensions )
{
    // count extensions and check for duplicates
  const char** iter;
  for (iter = extensions; *iter; ++iter)
  {
    const Reader* r = find_reader( *iter );
    if (r)
    {
       mbError->set_last_error( "Conflicting readers for file extension \"%s\":"
                               " \"%s\" and \"%s\".",
                               *iter, r->description.c_str(), description );
    }
  }
  readerList.push_back( Reader(function, description, extensions, iter - extensions) );
  return MB_SUCCESS;
}    
  
MBErrorCode MBReaderWriterSet::register_writer( writer_factory_t function,
                                                const char* description,
                                                const char** extensions )
{
    // count extensions and check for duplicates
  const char** iter;
  for (iter = extensions; *iter; ++iter)
  {
    const Writer* r = find_writer( *iter );
    if (r)
    {
       mbError->set_last_error( "Conflicting writers for file extension \"%s\":"
                               " \"%s\" and \"%s\".",
                               *iter, r->description.c_str(), description );
    }
  }
  writerList.push_back( Writer(function, description, extensions, iter - extensions) );
  return MB_SUCCESS;
}    

MBReaderIface* MBReaderWriterSet::get_file_extension_reader( 
                                  const std::string& filename ) const
{
  std::string ext = extension_from_filename( filename );
  const Reader* reader = find_reader( ext );
  return reader ? reader->factory(mbCore) : NULL;
}

MBWriterIface* MBReaderWriterSet::get_file_extension_writer( 
                                  const std::string& filename ) const
{
  std::string ext = extension_from_filename( filename );
  const Writer* writer = find_writer( ext );
  return writer ? writer->factory(mbCore) : NULL;
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

MBErrorCode MBReaderWriterSet::get_all_readers( std::vector<MBReaderIface*>& list ) const
{
  list.resize( readerList.size() );
  for (unsigned int i = 0; i < readerList.size(); ++i)
    list[i] = readerList[i].factory( mbCore );
  return MB_SUCCESS;
}

MBReaderWriterSet::Reader::Reader( reader_factory_t fact, 
                                   const char* desc, 
                                   const char** ext, int num_ext )
 : factory(fact), description(desc), extensions(num_ext)
{
  for (int i = 0; i < num_ext; ++i)
    extensions[i] = ext[i];
}

MBReaderWriterSet::Writer::Writer( writer_factory_t fact, 
                                   const char* desc, 
                                   const char** ext, int num_ext )
 : factory(fact), description(desc), extensions(num_ext)
{
  for (int i = 0; i < num_ext; ++i)
    extensions[i] = ext[i];
}

#ifdef WIN32
#define strcasecmp(A,B) _stricmp( A, B )
#endif

const MBReaderWriterSet::Reader* 
MBReaderWriterSet::find_reader( const std::string& ext ) const
{
  std::vector<Reader>::const_iterator iter;
  std::vector<std::string>::const_iterator siter;
  
    // try case-sensitive compare
  for (iter = readerList.begin(); iter != readerList.end(); ++iter)
    for (siter = iter->extensions.begin(); siter != iter->extensions.end(); ++siter)
      if (*siter == ext)
        return &*iter;
      
    // try case-insensitive compare
  for (iter = readerList.begin(); iter != readerList.end(); ++iter)
    for (siter = iter->extensions.begin(); siter != iter->extensions.end(); ++siter)
      if (0 == strcasecmp( siter->c_str(), ext.c_str() ))
        return &*iter;
  
  return NULL;
}

const MBReaderWriterSet::Writer* 
MBReaderWriterSet::find_writer( const std::string& ext ) const
{
  std::vector<Writer>::const_iterator iter;
  std::vector<std::string>::const_iterator siter;
  
    // try case-sensitive compare
  for (iter = writerList.begin(); iter != writerList.end(); ++iter)
    for (siter = iter->extensions.begin(); siter != iter->extensions.end(); ++siter)
      if (*siter == ext)
        return &*iter;
      
    // try case-insensitive compare
  for (iter = writerList.begin(); iter != writerList.end(); ++iter)
    for (siter = iter->extensions.begin(); siter != iter->extensions.end(); ++siter)
      if (0 == strcasecmp( siter->c_str(), ext.c_str() ))
        return &*iter;
  
  return NULL;
}
    
