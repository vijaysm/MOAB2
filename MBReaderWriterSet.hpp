/**
 *\class MBReaderWriterSet
 *\brief Maintain list of readers and writers.
 *\version 1.00
 *\date 2004-4-23
 *\author Jason Kraftcheck
 */

#ifndef MB_READER_WRITER_SET_HPP
#define MB_READER_WRITER_SET_HPP

#include <vector>
#include <string>
#include <map>
#include "MBInterface.hpp"

class MBReaderIface;
class MBWriterIface;
class MBCore;
class MBError;

class MBReaderWriterSet
{

  public:
    
    typedef MBReaderIface* (*reader_factory_t)( MBInterface* );
    typedef MBWriterIface* (*writer_factory_t)( MBInterface* );
  
    MBReaderWriterSet( MBCore* mdb, MBError* handler );
  
    ~MBReaderWriterSet();
  
    MBErrorCode register_reader( reader_factory_t function,
                                 const char* description,
                                 const char** extensions );
    MBErrorCode register_writer( writer_factory_t function,
                                 const char* description,
                                 const char** extensions );
    
    /** 
     * Create a reader object for the passed file name 
     * according to the dot-extension of the file name.
     * Caller must delete the object when finished.
     * Returns null if no matching file extension.
     */
    MBReaderIface* get_file_extension_reader( const std::string& filename ) const;

    /** 
     * Create a writer object for the passed file name 
     * according to the dot-extension of the file name.
     * Caller must delete the object when finished.
     * Returns null if no matching file extension.
     */
    MBWriterIface* get_file_extension_writer( const std::string& filename ) const;
    
    /**
     * Return an instance of every reader.
     * Caller must delete each instance when done.
     */
    MBErrorCode get_all_readers( std::vector<MBReaderIface*>& list ) const;
    
    static std::string extension_from_filename( const std::string& filename );
    
  private:
  
    struct Reader {
      Reader( reader_factory_t fact, const char* desc, const char** ext, int num_ext );
      reader_factory_t factory;
      std::string description;
      std::vector<std::string> extensions;
    };
   
    struct Writer {
      Writer( writer_factory_t fact, const char* desc, const char** ext, int num_ext );
      writer_factory_t factory;
      std::string description;
      std::vector<std::string> extensions;
    };
    
    const Reader* find_reader( const std::string& extension ) const;
    const Writer* find_writer( const std::string& extension ) const;
  
    MBCore* mbCore;
    MBError* mbError;
  
    std::vector<Reader> readerList;
    std::vector<Writer> writerList;
};


#endif
