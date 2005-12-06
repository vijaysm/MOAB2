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
 * \class ReadSTL
 * \brief ASCII and Binary Stereo Lithography File readers.
 * \author Jason Kraftcheck
 *
 * STL files contain no connectivity infomration.  Each triangle
 * is specified as by three sets of single-precision coordinate
 * triples.  This reader does not use ANY tolerance when comparing
 * vertex locations to recover connectivity.  The points must be
 * EXACTLY equal (including the sign on zero values.)  If the file
 * was written by an application which represented connectivity
 * explicitly, there is no reason for the vertex coordinates to 
 * be anything other than exactly equal.   
 *
 * For binary STL files, the defacto standard is that they be written
 * with a little-endian byte order.  The reader will attempt to 
 * determine the byte order automatically, and if it is ambiguous, 
 * will default to little-endian.  The byte ordering may be forced by
 * by creating an integer tag named "__STL_BYTE_ORDER" and setting a 
 * global/mesh value for the tag as 1 for big-endian or 0 for 
 * little-endian.
 *
 * For binary files, this reader relies on the file size to determine
 * the validity of the file and may use it in guessing the byte order.
 * This should not be an issue, as the file size can be determined
 * exactly from the number of triangles for a valid file.  However, if
 * for some reason the file is readable even though it is invalid (e.g.
 * it is some hybrid file with STL data in the beginning and some app-
 * specific data appended to the end of the file) the check on the file
 * size can be disabled by giving the reader a something other than a
 * regular file to read from.  For example, on Unix-like systems, have
 * the reader read from a FIFO instead of a file:
 *   mkfifo /tmp/fifo.stlb
 *   cat my_binary_file.stlb > /tmp/fifo.stlb 
 * and instruct the MOAB-based application to read from /tmp/fifo.stlb
 */


#ifndef READ_STL_HPP
#define READ_STL_HPP

#include <stdio.h>
#include "MBInterface.hpp"
#include "MBReaderIface.hpp"

class MBReadUtilIface;

// Base class for binary and ASCII readers
class ReadSTL : public MBReaderIface
{
   
public:

    //! factory method for binary STL reader
  static MBReaderIface* binary_instance( MBInterface* );
    //! factory method for ASCII STL reader.
  static MBReaderIface* ascii_instance( MBInterface* );

    //! Generic file loading code for both binary and ASCII readers.
    //! Calls reader-specific read_triangles function to do actual I/O.
  MBErrorCode load_file(const char *file_name,
                        const int* material_set_list,
                        const int num_material_sets );
  
    //! Constructor
  ReadSTL(MBInterface* impl = NULL);

   //! Destructor
  virtual ~ReadSTL();
                                 
  static long get_file_size( FILE* file );

    // An object to hold vertex coordinates, and an operator
    // for storing them in a STL tree-based container.
  struct Point {
    float coords[3];
    bool operator<( const Point& ) const;
  };
    // Data returned by read_triangles.
  struct Triangle {
    Point points[3];
  };

protected:
  
    // I/O-specific part of reader.  Read list of triangles
    // from file.
  virtual MBErrorCode read_triangles( const char* name, 
                                      std::vector<Triangle>& tris_out 
                                    ) = 0;

  MBReadUtilIface* readMeshIface;

    //! interface instance
  MBInterface* mdbImpl;

private:

    //! Generic file loading code for both binary and ASCII readers.
    //! Calls reader-specific read_triangles function to do actual I/O.
  MBErrorCode load_file_impl(const char *file_name,
                             const int* material_set_list,
                             const int num_material_sets );

    //! Meshset Handle for the mesh that is currently being read
  MBEntityHandle mCurrentMeshHandle;
};


//! Specialize ReadSTL for reading ASCII STL files
class ReadASCIISTL : public ReadSTL
{
public:
  
  ReadASCIISTL(MBInterface* impl = NULL) : ReadSTL(impl) {}

  virtual ~ReadASCIISTL() {}
  
protected:

  MBErrorCode read_triangles( const char* name, std::vector<Triangle>& tris_out );
};


//! Specialize ReadSTL for reading binary STL files
class ReadBinarySTL : public ReadSTL
{
public:

  ReadBinarySTL(MBInterface* impl = NULL) : ReadSTL(impl) {}

  virtual ~ReadBinarySTL() { }
  
protected:
  
  MBErrorCode read_triangles( const char* name, std::vector<Triangle>& tris_out );
};

#endif
