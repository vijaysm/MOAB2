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
 * \class WriteSTL
 * \brief ASCII and Binary Stereo Lithography File writers.
 * \author Jason Kraftcheck
 *
 * This writer will write only the MBTRI elements in the mesh.  It
 * will not decompose other 2-D elements into triangles, nor will
 * it skin the mesh or do any other high-level operation to generate
 * triangles from 3-D elements.  
 *
 * Binary files will be written with a little-endian byte order by
 * default.  The byte order can be controlled by creating an integer
 * tag named "__STL_BYTE_ORDER" and setting the global/mesh value to
 * 0 for little endian or 1 for big endian.
 */

#ifndef WRITE_STL_HPP
#define WRITE_STL_HPP

#include "MBForward.hpp"
#include "MBWriterIface.hpp"

#include <stdio.h>

class MBWriteUtilIface;

class WriteSTL : public MBWriterIface
{
 
public:
  
    //! factory method for binary STL writer
  static MBWriterIface* binary_instance( MBInterface* );
    //! factory method for ASCII STL writer.
  static MBWriterIface* ascii_instance( MBInterface* );

   //! Constructor
  WriteSTL(MBInterface *impl);

   //! Destructor
  virtual ~WriteSTL();
  
    //! writes out a file
  MBErrorCode write_file(const char *file_name,
                         const bool overwrite,
                         const MBEntityHandle *output_list,
                         const int num_sets,
                         std::vector<std::string>& qa_list,
                         int export_dimension);  

protected:
  
    //! Write list of triangles to an STL file.  
    //! Subclasses provide format-specific implementations 
    //! of this function.
  virtual MBErrorCode write_triangles( FILE* file,
                                       const char header[82],
                                       const MBRange& triangles ) = 0;

    //! Allow subclasses to request that file be opened in "binary"
    //! mode (a windows thing).
  virtual bool need_binary_io() const = 0;

    //! Given an array of vertex coordinates for a triangle,
    //! pass back individual point coordinates as floats and 
    //! calculate triangle normal.
  MBErrorCode get_triangle_data( const double vtx_coords[9],
                                 float v1[3],
                                 float v2[3],
                                 float v3[3],
                                 float n[3] );
                                       
    //! interface instance
  MBInterface *mbImpl;
  MBWriteUtilIface* mWriteIface;
  
private:

    //! Construct 80-byte, null-terminated description string from
    //! qa_list.  Unused space in header will be null-char padded.
  MBErrorCode make_header( char header[82], std::vector<std::string>& qa_list );
  
    //! Get triangles to write from input array of entity sets.  If
    //! no sets, gets all triangles.
  MBErrorCode get_triangles( const MBEntityHandle* set_array,
                             int set_array_length,
                             MBRange& triangles );  
  
    //! Open a file, respecting passed overwrite value and
    //! subclass-specified value for need_binary_io().
  FILE* open_file( const char* name, bool overwrite );
};


//! Specialize WriteSTL for writing ASCII STL Files.
class WriteASCIISTL : public WriteSTL
{
 
public:

  WriteASCIISTL(MBInterface *impl) : WriteSTL(impl) {}

  virtual ~WriteASCIISTL() {}

protected:

  virtual bool need_binary_io() const { return false; }
  
  virtual MBErrorCode write_triangles( FILE* file,
                                       const char header[81],
                                       const MBRange& triangles );
};


//! Specialize WriteSTL for writing binary STL Files.
class WriteBinarySTL : public WriteSTL
{
 
public:

  WriteBinarySTL(MBInterface *impl) : WriteSTL(impl) {}

  virtual ~WriteBinarySTL() {}

protected:

  virtual bool need_binary_io() const { return true; }
  
  virtual MBErrorCode write_triangles( FILE* file,
                                       const char header[81],
                                       const MBRange& triangles );
};

#endif
