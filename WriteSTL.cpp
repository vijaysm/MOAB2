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
 */


#include "WriteSTL.hpp"
#include "MBCN.hpp"

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <inttypes.h>
#include <unistd.h>
#include <math.h>

#ifdef _MSC_VER /* windows */
#  include <io.h>
#  define O_BINARY _O_BINARY
#  define O_CREAT  _O_CREAT
#  define O_EXCL   _O_EXCL
#  define O_TRUNC  _O_TRUNC
#  define O_WRONLY _O_WRONLY
#  define fdopen(A,B)  _fdopen( (A), (B) )
#  define open(A,B,C)  _open( (A), (B), (C) )
#else  /* posix */
#  include <fcntl.h>
#  define _S_IREAD  (S_IRUSR|S_IRGRP|S_IROTH)
#  define _S_IWRITE (S_IWUSR|S_IWGRP|S_IWOTH)
#endif

inline static uint32_t byte_swap( uint32_t value )
{
  return ((value & 0xFF000000) >> 24) |
         ((value & 0x00FF0000) >>  8) |
         ((value & 0x0000FF00) <<  8) |
         ((value & 0X000000FF) << 24);
}


inline static float byte_swap( float value )
{
  uint32_t bytes = byte_swap( *(uint32_t*)&value );
  return *(float*)&bytes;
}


inline static bool is_platform_little_endian()
{
  static const unsigned int one = 1;
  static const bool little = !*((char*)&one);
  return little;
}


MBWriterIface *WriteSTL::ascii_instance( MBInterface* iface )
  { return new WriteASCIISTL( iface ); }

MBWriterIface *WriteSTL::binary_instance( MBInterface* iface )
  { return new WriteBinarySTL( iface ); }

WriteSTL::WriteSTL(MBInterface *impl) 
    : mbImpl(impl)
{
  impl->query_interface("MBWriteUtilIface", reinterpret_cast<void**>(&mWriteIface));
}

WriteSTL::~WriteSTL() 
{
  mbImpl->release_interface("MBWriteUtilIface", mWriteIface);
}


MBErrorCode WriteSTL::write_file(const char *file_name, 
                                 const bool overwrite,
                                 const MBEntityHandle *ent_handles,
                                 const int num_sets,
                                 std::vector<std::string>& qa_list, 
                                 int  )
{
  char header[81];
  MBRange triangles;
  MBErrorCode rval;
  
  rval = make_header( header, qa_list );
  if (MB_SUCCESS != rval)
    return rval;
  
  rval = get_triangles( ent_handles, num_sets, triangles );
  if (MB_SUCCESS != rval)
    return rval;
    
  FILE* file = open_file( file_name, overwrite );
  if (!file)
    return MB_FILE_DOES_NOT_EXIST; 
  
  rval = write_triangles( file, header, triangles );
  fclose( file );
  return rval;
}


FILE* WriteSTL::open_file( const char* name, bool overwrite )
{
    // Open file with write access, and create it if it doesn't exist.
  int flags = O_WRONLY|O_CREAT;
    // Select behavior if the named file already exists.  If 
    // overwrite is true, truncate the file.  If it is false,
    // make the call to open() fail.
  if (overwrite)
    flags |= O_TRUNC;
  else
    flags |= O_EXCL;
    // If platform defines a "binary" bit in the file access
    // flags (i.e. we're building on windows), then set it
    // if we're writing a binary file.
#ifdef O_BINARY
  if (need_binary_io())
    flags |= O_BINARY;
#endif

    // Give everyone read and write, but not execute, permision.
    // These are not the final permisions for the file.  Permissions
    // are removed according to the user's umask.  All we want to
    // say here is that the executable bits should not be set because
    // this isn't an executable file.  Everything else is a user
    // preference and should be left up to the umask.
  int creat_mode = _S_IREAD|_S_IWRITE;

    // Open the file.
  int fd = open( name, flags, creat_mode );
  if (fd < 0)
  {
    mWriteIface->report_error( "%s: %s\n", name, strerror(errno) );
    return 0;
  }
  FILE* result = fdopen( fd, need_binary_io() ? "wb": "w" );
  if (!result)
    close( fd );
  
  return result;
}

MBErrorCode WriteSTL::make_header( char header[81], std::vector<std::string>& qa_list )
{
  memset( header, 0, 81 );
  
  std::string result;
  for (std::vector<std::string>::iterator i = qa_list.begin(); i != qa_list.end(); ++i)
  {
    result += " ";
    result += *i;
  }
  
  size_t len = result.size();
  if (len > 80)
    len = 80;
  memcpy( header, result.c_str(), len );
  
  return MB_SUCCESS;
}

MBErrorCode WriteSTL::get_triangles( const MBEntityHandle* set_array,
                                     int set_array_length,
                                     MBRange& triangles )
{
  if (!set_array || set_array_length == 0)
  {
    return mbImpl->get_entities_by_type( 0, MBTRI, triangles );
  }
  
  const MBEntityHandle* iter = set_array;
  const MBEntityHandle* end = iter + set_array_length;
  for (; iter != end; ++iter)
  {
    MBRange r;
    MBErrorCode rval = mbImpl->get_entities_by_type( *iter, MBTRI, r, true );
    if (MB_SUCCESS != rval)
      return rval;
    triangles.merge( r );
  }
  
  return MB_SUCCESS;
}

MBErrorCode WriteSTL::get_triangle_data( const double coords[9],
                                         float v1[3],
                                         float v2[3],
                                         float v3[3],
                                         float n[3] )
{
  float e1[3], e2[3];
  v1[0] = coords[0];
  v1[1] = coords[1];
  v1[2] = coords[2];
  v2[0] = coords[3];
  v2[1] = coords[4];
  v2[2] = coords[5];
  v3[0] = coords[6];
  v3[1] = coords[7];
  v3[2] = coords[8];
  e1[0] = v2[0] - v1[0];
  e1[1] = v2[1] - v1[1];
  e1[2] = v2[2] - v1[2];
  e2[0] = v3[0] - v1[0];
  e2[1] = v3[1] - v1[1];
  e2[2] = v3[2] - v1[2];
  n[0] = e1[1]*e2[2] - e1[2]*e2[1];
  n[1] = e1[2]*e2[0] - e1[0]*e2[2];
  n[2] = e1[0]*e2[1] - e1[1]*e2[0];
  float inv_len = 1.0f / sqrtf( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );
  n[0] *= inv_len;
  n[1] *= inv_len;
  n[2] *= inv_len;
  return MB_SUCCESS;
}


MBErrorCode WriteASCIISTL::write_triangles( FILE* file,
                                            const char header[81],
                                            const MBRange& triangles )
{
  const char solid_name[] = "MOAB";
  
  char myheader[81] = "solid ";
  strcat( myheader, solid_name );
  strncat( myheader, header, 80 );
  
  if (EOF == fputs( myheader, file ) || EOF == fputs( "\n", file ))
    return MB_FILE_WRITE_ERROR;
  
  MBErrorCode rval;
  double coords[9];
  float v1[3], v2[3], v3[3];
  float n[3];
  for (MBRange::const_iterator iter = triangles.begin();
       iter != triangles.end(); ++iter)
  {
    const MBEntityHandle* conn;
    int num_vtx;
    
    rval = mbImpl->get_connectivity( *iter, conn, num_vtx );
    if (MB_SUCCESS != rval)
      return rval;
    if (num_vtx != 3)
      return MB_FAILURE;
    
    rval = mbImpl->get_coords( conn, 3, coords );
    if (MB_SUCCESS != rval)
      return rval;
    
    rval = get_triangle_data( coords, v1, v2, v3, n );
    if (MB_SUCCESS != rval)
      return rval;
   
    fprintf( file, "facet normal %e %e %e\n", n[0], n[1], n[2] );
    fprintf( file,"outer loop\n" );
    fprintf( file,"vertex %e %e %e\n", v1[0], v1[1], v1[2] );
    fprintf( file,"vertex %e %e %e\n", v2[0], v2[1], v2[2] );
    fprintf( file,"vertex %e %e %e\n", v3[0], v3[1], v3[2] );
    fprintf( file,"endloop\n" );
    fprintf( file,"endfacet\n" );
  }
  
  fprintf( file,"endsolid %s\n", solid_name );
  return MB_SUCCESS;
}

struct BinTri
{
  float normal[3];
  float vertex1[3];
  float vertex2[3];
  float vertex3[3];
  char pad[2];
};

static inline void byte_swap( float vect[3] )
{
  vect[0] = byte_swap( vect[0] );
  vect[1] = byte_swap( vect[1] );
  vect[2] = byte_swap( vect[2] );
}

MBErrorCode WriteBinarySTL::write_triangles( FILE* file,
                                             const char header[81],
                                             const MBRange& triangles )
{
  MBErrorCode rval;
  if (fwrite( header, 80, 1, file ) != 1)
    return MB_FILE_WRITE_ERROR;
  
  bool swap_bytes = !is_platform_little_endian();  // default to little endian

    // Check for tag specifying file byte order
  MBTag bo_tag = 0;
  rval = mbImpl->tag_get_handle( "__STL_BYTE_ORDER", bo_tag );
  if (MB_SUCCESS == rval)
  {
    int value;
    rval = mbImpl->tag_get_data( bo_tag, 0, 1, &value );
    if (MB_SUCCESS != rval) 
      return rval;
    bool is_file_little_endian = (0 == value);
    swap_bytes = (is_platform_little_endian() != is_file_little_endian);
  } 
  else if (MB_TAG_NOT_FOUND != rval)
    return rval;
  
  int32_t count = swap_bytes ? byte_swap(triangles.size()) : triangles.size();
  if (fwrite( &count, 4, 1, file ) != 1)
    return MB_FILE_WRITE_ERROR;

  double coords[9];
  BinTri tri;
  for (MBRange::const_iterator iter = triangles.begin();
       iter != triangles.end(); ++iter)
  {
    const MBEntityHandle* conn;
    int num_vtx;
    
    rval = mbImpl->get_connectivity( *iter, conn, num_vtx );
    if (MB_SUCCESS != rval)
      return rval;
    if (num_vtx != 3)
      return MB_FAILURE;
    
    rval = mbImpl->get_coords( conn, 3, coords );
    if (MB_SUCCESS != rval)
      return rval;
    
    rval = get_triangle_data( coords, tri.vertex1, tri.vertex2, tri.vertex3, tri.normal );
    if (MB_SUCCESS != rval)
      return rval;
    
    if (swap_bytes)
    {
      byte_swap( tri.normal );
      byte_swap( tri.vertex1 );
      byte_swap( tri.vertex2 );
      byte_swap( tri.vertex3 );
    }
   
    if (1 != fwrite( &tri, 50, 1, file ))
      return MB_FILE_WRITE_ERROR;
  }
  
  return MB_SUCCESS;
}

