/*
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
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

#include "moab/FileOptions.hpp"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <iostream>

using namespace moab;

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
  FileOptions tool( "INT1=1;NUL1;STR1=ABC;DBL1=1.0;dbl2=2.0;DBL3=3.0;INT2=2;nul2;NUL3;INT3=3;str2=once upon a time;str3==fubar=;;INTS=1-3,5,6;DBLS=1.0,2.0, 3.0;STRS=var1, var2_var2;STRS2=" );

  std::string s;
  int i;
  double d;
  ErrorCode rval;
  
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
  
    // We haven't looked at all of the options yet
  EQUAL( false, tool.all_seen() );
  rval = tool.get_unseen_option( s );
  CHECK( rval );
  EQUAL( s, "str3" );
  
    // test options using generic get_option method
    
  rval = tool.get_option( "NUL3", s );
  CHECK( rval );
  EQUAL( s.empty(), true );
  
  rval = tool.get_option( "STR3", s );
  CHECK( rval );
  EQUAL( s, "=fubar=" );
  
    // test size of options string
  unsigned l = tool.size();
  EQUAL( l, 16u );
  
    // test ints option
  std::vector<int> ivals;
  rval = tool.get_ints_option("INTS", ivals);
  CHECK( rval );
  EQUAL(5, ivals.size());
  EQUAL(1, ivals[0]);
  EQUAL(2, ivals[1]);
  EQUAL(3, ivals[2]);
  EQUAL(5, ivals[3]);
  EQUAL(6, ivals[4]);

    // test dbls option
  std::vector<double> vals;
  rval = tool.get_reals_option("DBLS", vals);
  CHECK( rval );
  EQUAL(3, vals.size());
  EQUAL(1.0, vals[0]);
  EQUAL(2.0, vals[1]);
  EQUAL(3.0, vals[2]);
  
    // test strs option
  std::vector<std::string> svals;
  rval = tool.get_strs_option("STRS", svals);
  CHECK( rval );
  EQUAL(2, svals.size());
  EQUAL("var1", svals[0]);
  EQUAL("var2_var2", svals[1]);
  
  svals.clear();
  rval = tool.get_strs_option("STRS2", svals);
  EQUAL( MB_TYPE_OUT_OF_RANGE, rval );
  
    // We requested every option
  EQUAL( true, tool.all_seen() );
  rval = tool.get_unseen_option( s );
  EQUAL( MB_ENTITY_NOT_FOUND, rval );
  
    // test alternate separator
  
  FileOptions tool2( ";+OPT1=ABC+OPT2=" );
  l = tool2.size();
  EQUAL( l, 2 );

    // We haven't looked at all of the options yet
  EQUAL( false, tool2.all_seen() );
  rval = tool2.get_unseen_option( s );
  CHECK( rval );
  EQUAL( s, "OPT1" );
   
  rval = tool2.get_option( "opt1", s );
  CHECK( rval );
  EQUAL( s, "ABC" );
  
  rval = tool2.get_option( "opt2", s );
  CHECK( rval );
  bool e = s.empty();
  EQUAL( e, true );
  
  l = tool2.size();
  EQUAL( l, 2 );
  
    // We requested every option
  EQUAL( true, tool2.all_seen() );
  rval = tool2.get_unseen_option( s );
  EQUAL( MB_ENTITY_NOT_FOUND, rval );
  
    
    // test empty options string
    
  FileOptions tool3( ";;;;" );
  e = tool3.empty();
  EQUAL( e, true );
  l = tool3.size();
  EQUAL( l, 0 );
  EQUAL( true, tool3.all_seen() );
  
  FileOptions tool4(NULL);
  e = tool4.empty();
  EQUAL( e, true );
  l = tool4.size();
  EQUAL( l, 0 );
  EQUAL( true, tool4.all_seen() );
  
  FileOptions tool5(";+");
  e = tool5.empty();
  EQUAL( e, true );
  l = tool5.size();
  EQUAL( l, 0 );
  EQUAL( true, tool5.all_seen() );
  
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

