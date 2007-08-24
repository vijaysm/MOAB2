/*
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

/**\file FileOptions.cpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2007-08-21
 */

#include "FileOptions.hpp"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

const char DEFAULT_SEPARATOR = ';';

static inline bool strempty( const char* s ) { return !*s; }

FileOptions::FileOptions( const char* str )
  : mData(0)
{
    // if option string is null, just return
  if (!str)
    return;
  
    // check if alternate separator is specified
  char separator[2] = { DEFAULT_SEPARATOR, '\0' };
  if (*str == DEFAULT_SEPARATOR) {
    ++str;
    if (strempty(str))
      return;
    separator[0] = *str;
    ++str;
  }
  
    // don't bother allocating copy of input string if
    // input string is empty.
  if (!strempty(str))
  {
       // tokenize at separator character
    mData = strdup( str );
    for (char* i = strtok( mData, separator ); i; i = strtok( 0, separator )) 
      if (!strempty(i)) // skip empty strings
        mOptions.push_back( i );
  }
}

FileOptions::FileOptions( const FileOptions& copy ) :
  mData(0), mOptions( copy.mOptions.size() )
{
  if (!copy.mOptions.empty()) {
    const char* last = copy.mOptions.back();
    const char* endptr = last + strlen(last) + 1;
    size_t len = endptr - copy.mData;
    mData = (char*)malloc( len );
    memcpy( mData, copy.mData, len );
    for (size_t i = 0; i < mOptions.size(); ++i)
      mOptions[i] = mData + (copy.mOptions[i] - copy.mData);
  }
}

FileOptions& FileOptions::operator=( const FileOptions& copy )
{
  free( mData );
  mData = 0;
  mOptions.resize( copy.mOptions.size() );

  if (!copy.mOptions.empty()) {
    const char* last = copy.mOptions.back();
    const char* endptr = last + strlen(last) + 1;
    size_t len = endptr - copy.mData;
    mData = (char*)malloc( len );
    memcpy( mData, copy.mData, len );
    for (size_t i = 0; i < mOptions.size(); ++i)
      mOptions[i] = mData + (copy.mOptions[i] - copy.mData);
  }
    
  return *this;
}

FileOptions::~FileOptions()
{
  free( mData );
}

MBErrorCode FileOptions::get_null_option( const char* name ) const
{
  const char* s;
  MBErrorCode rval = get_option( name, s );
  if (MB_SUCCESS != rval)
    return rval;
  return strempty(s) ? MB_SUCCESS : MB_TYPE_OUT_OF_RANGE;
}

MBErrorCode FileOptions::get_int_option( const char* name, int& value ) const
{
  const char* s;
  MBErrorCode rval = get_option( name, s );
  if (MB_SUCCESS != rval)
    return rval;
  
    // empty string
  if (strempty(s))
    return MB_TYPE_OUT_OF_RANGE;
  
    // parse value
  char* endptr;
  long int pval = strtol( s, &endptr, 0 );
  if (!strempty(endptr)) // syntax error
    return MB_TYPE_OUT_OF_RANGE;
  
    // check for overflow (parsing long int, returning int)
  value = pval;
  if (pval != (long int)value)
    return MB_TYPE_OUT_OF_RANGE;
  
  return MB_SUCCESS;
}

MBErrorCode FileOptions::get_real_option ( const char* name, double& value ) const
{
  const char* s;
  MBErrorCode rval = get_option( name, s );
  if (MB_SUCCESS != rval)
    return rval;
  
    // empty string
  if (strempty(s))
    return MB_TYPE_OUT_OF_RANGE;
  
    // parse value
  char* endptr;
  value = strtod( s, &endptr );
  if (!strempty(endptr)) // syntax error
    return MB_TYPE_OUT_OF_RANGE;
  
  return MB_SUCCESS;
}

MBErrorCode FileOptions::get_str_option( const char* name, std::string& value ) const
{
  const char* s;
  MBErrorCode rval = get_option( name, s );
  if (MB_SUCCESS != rval)
    return rval;
  if (strempty(s))
    return MB_TYPE_OUT_OF_RANGE;
  value = s;
  return MB_SUCCESS;
}

MBErrorCode FileOptions::get_option( const char* name, std::string& value ) const
{
  const char* s;
  MBErrorCode rval = get_option( name, s );
  if (MB_SUCCESS != rval)
    return rval;
  
  value = s;
  return MB_SUCCESS;
}  

MBErrorCode FileOptions::get_option( const char* name, const char*& value ) const
{
  std::vector<const char*>::const_iterator i;
  for (i = mOptions.begin(); i != mOptions.end(); ++i) {
    const char* opt = *i;
    if (compare( name, opt )) {
      value = opt + strlen(name);
        // if compare returned true, next char after option
        // name must be either the null char or an equals symbol.
      if (*value == '=') 
        ++value;
        
      return MB_SUCCESS;
    }
  }
  
  return MB_ENTITY_NOT_FOUND;
}

bool FileOptions::compare( const char* name, const char* option )
{
  while (!strempty(name) && toupper(*name) == toupper(*option)) {
    ++name;
    ++option;
  }
   // match if name matched option for length of name,
   // and option either matched entirely or matches up to
   // and equals sign.
  return strempty(name) && (strempty(option) || *option == '=');
}

void FileOptions::get_options( std::vector<std::string>& list ) const
{
  list.clear();
  list.resize( mOptions.size() );
  std::copy( mOptions.begin(), mOptions.end(), list.begin() );
}

#ifdef TEST

#include <iostream>

#define CHECK(A) \
  if (MB_SUCCESS != (A)) { \
    std::cerr << "Failure at line " << __LINE__ << ": error code " << (A) << std::endl; \
    return 1; \
  }

#define EQUAL(A,B) \
  if (A != B) { \
    std::cerr << "Failure at line " << __LINE__ << ": expected " << (B) << " but got " << (A) << std::endl; \
    return 2; \
  }

int main()
{
  FileOptions tool( "INT1=1;NUL1;STR1=ABC;DBL1=1.0;dbl2=2.0;DBL3=3.0;INT2=2;nul2;NUL3;INT3=3;str2=once upon a time;str3==fubar=;;" );

  std::string s;
  int i;
  double d;
  MBErrorCode rval;
  
    // test basic get_option method without deleting entry
  rval = tool.get_option( "STR1", s );
  CHECK(rval);
  EQUAL( s, "ABC" );
  
    // test basic get_option method again, this time deleting the entry
  rval = tool.get_option( "STR1", s );
  CHECK(rval);
  EQUAL( s, "ABC" );
  
    // test basig get_option method with a null option
  rval = tool.get_option( "NUL2", s );
  CHECK( rval );
  EQUAL( s.empty(), true );

  
    // test null option
  rval = tool.get_null_option( "nul1" );
  CHECK( rval );
  
    // try null option method on non-null value
  rval = tool.get_null_option( "INT1" ) ;
  EQUAL( rval, MB_TYPE_OUT_OF_RANGE) ;
  

    // test integer option
  rval = tool.get_int_option( "int1", i );
  CHECK( rval );
  EQUAL( i, 1 );
  
  rval = tool.get_int_option( "int2", i );
  CHECK( rval );
  EQUAL( i, 2 );
  
    // test integer option on non-integer value
  rval = tool.get_int_option( "dbl2", i );
  EQUAL( rval, MB_TYPE_OUT_OF_RANGE );
  
    // test integer option on null value
  rval = tool.get_int_option( "NUL3", i);
  EQUAL( rval, MB_TYPE_OUT_OF_RANGE );
  
    // test double option
  rval = tool.get_real_option( "dbl1", d );
  CHECK( rval );
  EQUAL( d, 1.0 );
  
  rval = tool.get_real_option( "dbl2", d );
  CHECK( rval );
  EQUAL( d, 2.0 );
  
  rval = tool.get_real_option( "int3", d );
  CHECK( rval );
  EQUAL( d, 3.0 );
  
    // test real option on non-real value
  rval = tool.get_real_option( "str2", d );
  EQUAL( rval, MB_TYPE_OUT_OF_RANGE );
  
  
    // test real option on null value
  rval = tool.get_real_option( "NUL3", d );
  EQUAL( rval, MB_TYPE_OUT_OF_RANGE );
  
    // test get a simple string option
  rval = tool.get_str_option( "DBL3", s );
  CHECK( rval );
  EQUAL( s, "3.0" );
  
    // test get a string with spaces
  rval = tool.get_str_option("STR2", s );
  CHECK( rval );
  EQUAL( s, "once upon a time" );
  
    // try to get a string value for a null option
  rval = tool.get_str_option( "nul3", s );
  EQUAL( rval, MB_TYPE_OUT_OF_RANGE );
  
    // test options using generic get_option method
    
  rval = tool.get_option( "NUL3", s );
  CHECK( rval );
  EQUAL( s.empty(), true );
  
  rval = tool.get_option( "STR3", s );
  CHECK( rval );
  EQUAL( s, "=fubar=" );
  
    // test size of options string
  unsigned l = tool.size();
  EQUAL( l, 12u );
  
  
    // test alternate separator
  
  FileOptions tool2( ";+OPT1=ABC+OPT2=" );
  l = tool2.size();
  EQUAL( l, 2 );
  
  rval = tool2.get_option( "opt1", s );
  CHECK( rval );
  EQUAL( s, "ABC" );
  
  rval = tool2.get_option( "opt2", s );
  CHECK( rval );
  bool e = s.empty();
  EQUAL( e, true );
  
  l = tool2.size();
  EQUAL( l, 2 );
  
    
    // test empty options string
    
  FileOptions tool3( ";;;;" );
  e = tool3.empty();
  EQUAL( e, true );
  l = tool3.size();
  EQUAL( l, 0 );
  
  FileOptions tool4(NULL);
  e = tool4.empty();
  EQUAL( e, true );
  l = tool4.size();
  EQUAL( l, 0 );
  
  FileOptions tool5(";+");
  e = tool5.empty();
  EQUAL( e, true );
  l = tool5.size();
  EQUAL( l, 0 );
  
    // test copy constructor
  
  FileOptions tool6( tool2 );
  
  rval = tool6.get_option( "opt1", s );
  CHECK( rval );
  EQUAL( s, "ABC" );
  
  rval = tool6.get_option( "opt2", s );
  CHECK( rval );
  e = s.empty();
  EQUAL( e, true );
  
  l = tool6.size();
  EQUAL( l, 2 );
  
  FileOptions tool7( tool5 );
  e = tool7.empty();
  EQUAL( e, true );
  l = tool7.size();
  EQUAL( l, 0 );
  
    // test assignment operator
  
  FileOptions tool8( tool2 );
  tool8 = tool;
  EQUAL( tool8.size(), tool.size() );
    
  return 0;
}

#endif
