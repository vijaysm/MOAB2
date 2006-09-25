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
 *\class MBReaderWriterSet
 *\brief Maintain list of readers and writers.
 *\version 1.00
 *\date 2004-4-23
 *\author Jason Kraftcheck
 */

#ifndef MB_READER_WRITER_SET_HPP
#define MB_READER_WRITER_SET_HPP

#include <list>
#include <string>
#include "MBTypes.h"

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
    
    /**
     * Regiseter a reader and/or writer
     * Either factory function may be NULL, but not both.
     *
     *\param reader_fact  A factory method to create an instance of the reader
     *\param writer_fact  A factory method to create an instance of the reader
     *\param description  A short description of the file format.
     *\param extensions   A null-terminated list of file extensions
     */
    MBErrorCode register_factory( reader_factory_t reader_fact,
                                  writer_factory_t writer_fact,
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
     * Get the file extension from a file name
     */
    static std::string extension_from_filename( const std::string& filename );
  
    class Handler {
      
      friend class MBReaderWriterSet;
      
      public:
      
      Handler( reader_factory_t read_f,
               writer_factory_t write_f,
               const char* desc, 
               const char** ext, 
               int num_ext );
      
      inline const std::string& description() const { return mDescription; }
      inline void get_extensions( std::vector<std::string>& list_out ) const
        { list_out = mExtensions; }
      
      inline bool have_reader() const { return NULL != mReader; }
      inline bool have_writer() const { return NULL != mWriter; }
      
      inline MBReaderIface* make_reader( MBInterface* iface ) const
        { return have_reader() ? mReader(iface) : NULL; }
      
      inline MBWriterIface* make_writer( MBInterface* iface ) const
        { return have_writer() ? mWriter(iface) : NULL; }
      
      private:
      
      reader_factory_t mReader;
      writer_factory_t mWriter;
      
      std::string mDescription;
      std::vector<std::string> mExtensions;
    };
    
    typedef std::list<Handler>::const_iterator iter_type;
    
    inline iter_type begin() const { return handlerList.begin(); }
    
    inline iter_type end()   const { return handlerList.end();   }
    
    
    iter_type handler_from_extension( const std::string& extension,
                                      bool with_reader = false, 
                                      bool with_writer = false) const;
    
  private:
  
    MBCore* mbCore;
    MBError* mbError;
  
    std::list<Handler> handlerList;
};


#endif
