/**
 *\class MBReaderIface
 *\brief Interface for mesh reader implementations.
 *\version 1.00
 *\date 2004-4-23
 *\author Jason Kraftcheck
 */

#ifndef MB_READER_IFACE_HPP
#define MB_READER_IFACE_HPP

#include <vector>
#include <string>
#include "MBInterface.hpp"

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
     *\param material_set_list   A list of material sets to read, or NULL
     *                           if the entire file is to be read.
     *\param material_set_list_len The length of <code>material_set_list</code>
     *\author Jason Kraftcheck
     */
    virtual MBErrorCode load_file( const char* file_name,
                                   const int* material_set_list,
                                   const int material_set_list_len ) = 0;
};

#endif

    
