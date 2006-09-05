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
 * \class MBAffineXform
 * \brief Define an affine transformatino
 * \author Jason Kraftcheck (kraftche@cae.wisc.edu)
 * \date August, 2006
 */

#include "MBAffineXform.hpp"
#include <assert.h>

#ifndef TEST

const char* const AFFINE_XFORM_TAG_NAME = "AFFINE_TRANSFORM";

MBErrorCode MBAffineXform::get_tag( MBTag& tag_out,
                                    MBInterface* interface,
                                    const char* tagname )
{
  assert( sizeof(MBAffineXform) == 12*sizeof(double) );
  
  if (!tagname)
    tagname = AFFINE_XFORM_TAG_NAME;
 
  MBErrorCode rval;
  
  rval = interface->tag_get_handle( tagname, tag_out );
  
  if (MB_TAG_NOT_FOUND == rval) 
    return interface->tag_create( tagname, 
                                  sizeof(MBAffineXform),
                                  MB_TAG_DENSE,
                                  MB_TYPE_DOUBLE,
                                  tag_out,
                                  0 );
  if (MB_SUCCESS != rval)
    return rval;
  
  int size;
  rval = interface->tag_get_size( tag_out, size );
  if (MB_SUCCESS != rval || size != sizeof(MBAffineXform))
    return MB_FAILURE;
  
  MBDataType type;
  rval = interface->tag_get_data_type( tag_out, type );
  if (MB_SUCCESS != rval || type != MB_TYPE_DOUBLE)
    return MB_FAILURE;
  
  return MB_SUCCESS;
}

#else // ******************* Unit Test code ***********************

#include <iostream>
#define ASSERT_VECTORS_EQUAL(A, B) assert_vectors_equal( (A), (B), #A, #B, __LINE__ )
#define ASSERT_DOUBLES_EQUAL(A, B) assert_doubles_equal( (A), (B), #A, #B, __LINE__ )
#define ASSERT(B) assert_bool( (B), #B, __LINE__ )

const double TOL = 1e-6;

int error_count = 0;

void assert_vectors_equal( const MBCartVect& a, const MBCartVect& b, 
                           const char* sa, const char* sb,
                           int lineno )
{
  if (fabs(a[0] - b[0]) > TOL ||
      fabs(a[1] - b[1]) > TOL ||
      fabs(a[2] - b[2]) > TOL) {
    std::cerr << "Assertion failed at line " << lineno << std::endl
              << "\t" << sa << " == " << sb << std::endl
              << "\t[" << a[0] << ", " << a[1] << ", " << a[2] << "] == ["
              << b[0] << ", " << b[1] << ", " << b[2] << "]" << std::endl;
    ++error_count;
  }
}

void assert_doubles_equal( double a, double b, const char* sa, const char* sb, int lineno )
{
  if (fabs(a - b) > TOL) {
    std::cerr << "Assertion failed at line " << lineno << std::endl
              << "\t" << sa << " == " << sb << std::endl
              << "\t" << a << " == " << b << std::endl;
    ++error_count;
  }
}

void assert_bool( bool b, const char* sb, int lineno )
{
  if (!b) {
    std::cerr << "Assertion failed at line " << lineno << std::endl
              << "\t" << sb << std::endl;
    ++error_count;
  }
}


const MBCartVect point1( 0.0, 0.0, 0.0 ), point2( 3.5, 1000, -200 );
const MBCartVect vect1( 0.0, 0.0, -100.0 ), vect2( 1.0, 0.0, 1.0 );


void test_none()
{
    // default xform should do nothing.
  MBCartVect output;
  MBAffineXform none;
  none.xform_point( point1.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, point1 );
  none.xform_point( point2.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, point2 );
  none.xform_vector( vect1.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, vect1 );
  none.xform_vector( vect2.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, vect2 );
}

void test_translation()
{
  MBCartVect offset( 1.0, 2.0, 3.0 );
  MBCartVect output;
  
  MBAffineXform move = MBAffineXform::translation( offset.array() );
  
  // test that points are moved by offset
  move.xform_point( point1.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, point1 + offset );
  move.xform_point( point2.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, point2 + offset );
  
  // vectors should not be changed by a translation
  move.xform_vector( vect1.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, vect1 );
  move.xform_vector( vect2.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, vect2 );
}

void test_rotation()
{
  MBCartVect output;
  
  // rotate 90 degress about Z axis
  
  MBAffineXform rot = MBAffineXform::rotation( M_PI/2.0, MBCartVect(0,0,1).array() );
  ASSERT_DOUBLES_EQUAL( rot.matrix().determinant(), 1.0 );
  
  rot.xform_point( point1.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, point1 ); // origin not affected by transform
  
  MBCartVect expectedz( -point2[1], point2[0], point2[2] ); // in first quadrant
  rot.xform_point( point2.array(), output.array() );
  ASSERT_DOUBLES_EQUAL( output.length(), point2.length() );
  ASSERT_VECTORS_EQUAL( output, expectedz );
  
  rot.xform_vector( vect1.array(), output.array() );
  ASSERT_DOUBLES_EQUAL( output.length(), vect1.length() );
  ASSERT_VECTORS_EQUAL( output, vect1 );
  
  rot.xform_vector( vect2.array(), output.array() );
  ASSERT_DOUBLES_EQUAL( output.length(), vect2.length() );
  ASSERT_VECTORS_EQUAL( output, MBCartVect( 0, 1, 1 ) );
  
  // rotate 90 degress about Y axis
  
  rot = MBAffineXform::rotation( M_PI/2.0, MBCartVect(0,1,0).array() );
  ASSERT_DOUBLES_EQUAL( rot.matrix().determinant(), 1.0 );
  
  rot.xform_point( point1.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, point1 ); // origin not affected by transform
  
  MBCartVect expectedy( point2[2], point2[1], -point2[0] ); // in second quadrant
  rot.xform_point( point2.array(), output.array() );
  ASSERT_DOUBLES_EQUAL( output.length(), point2.length() );
  ASSERT_VECTORS_EQUAL( output, expectedy );
  
  rot.xform_vector( vect1.array(), output.array() );
  ASSERT_DOUBLES_EQUAL( output.length(), vect1.length() );
  ASSERT_VECTORS_EQUAL( output, MBCartVect(-100,0,0) );
  
  rot.xform_vector( vect2.array(), output.array() );
  ASSERT_DOUBLES_EQUAL( output.length(), vect2.length() );
  ASSERT_VECTORS_EQUAL( output, MBCartVect( 1, 0, -1 ) );
  
  // rotate 90 degress about X axis
  
  rot = MBAffineXform::rotation( M_PI/2.0, MBCartVect(1,0,0).array() );
  ASSERT_DOUBLES_EQUAL( rot.matrix().determinant(), 1.0 );
  
  rot.xform_point( point1.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, point1 ); // origin not affected by transform
  
  MBCartVect expectedx( point2[0], -point2[2], point2[1] ); // in third quadrant
  rot.xform_point( point2.array(), output.array() );
  ASSERT_DOUBLES_EQUAL( output.length(), point2.length() );
  ASSERT_VECTORS_EQUAL( output, expectedx );
  
  rot.xform_vector( vect1.array(), output.array() );
  ASSERT_DOUBLES_EQUAL( output.length(), vect1.length() );
  ASSERT_VECTORS_EQUAL( output, MBCartVect(0,100,0) );
  
  rot.xform_vector( vect2.array(), output.array() );
  ASSERT_DOUBLES_EQUAL( output.length(), vect2.length() );
  ASSERT_VECTORS_EQUAL( output, MBCartVect( 1, -1, 0 ) );
  
  // rotate 180 degrees about vector in XY plane
  
  rot = MBAffineXform::rotation( M_PI, MBCartVect( 1, 1, 0 ).array() );
  ASSERT_DOUBLES_EQUAL( rot.matrix().determinant(), 1.0 );
  
  rot.xform_point( point1.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, point1 ); // origin not affected by transform
  
  rot.xform_point( point2.array(), output.array() );
  ASSERT_DOUBLES_EQUAL( output.length(), point2.length() );
  ASSERT_VECTORS_EQUAL( output, MBCartVect( point2[1], point2[0], -point2[2] ) );
  
  rot.xform_vector( vect1.array(), output.array() );
  ASSERT_DOUBLES_EQUAL( output.length(), vect1.length() );
  ASSERT_VECTORS_EQUAL( output, -vect1 ); // vector is in xy plane
  
  rot.xform_vector( vect2.array(), output.array() );
  ASSERT_DOUBLES_EQUAL( output.length(), vect2.length() );
  ASSERT_VECTORS_EQUAL( output, MBCartVect( 0, 1, -1 ) );
}

MBCartVect refl( const MBCartVect& vect, const MBCartVect& norm )
{
  MBCartVect n(norm);
  n.normalize();
  double d = vect % n;
  return vect - 2 * d * n;
}

void test_reflection()
{
  MBCartVect output;
  
  // reflect about XY plane
  MBAffineXform ref = MBAffineXform::reflection( MBCartVect( 0, 0, 1 ).array() );
  ASSERT_DOUBLES_EQUAL( ref.matrix().determinant(), -1.0 );
  ref.xform_point( point1.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, point1 );
  ref.xform_point( point2.array(), output.array() );
  ASSERT_DOUBLES_EQUAL( output.length(), point2.length() );
  ASSERT_VECTORS_EQUAL( output, MBCartVect(point2[0],point2[1],-point2[2]) );
  ref.xform_vector( vect1.array(), output.array() );
  ASSERT_DOUBLES_EQUAL( output.length(), vect1.length() );
  ASSERT_VECTORS_EQUAL( output, -vect1 );
  ref.xform_vector( vect2.array(), output.array() );
  ASSERT_DOUBLES_EQUAL( output.length(), vect2.length() );
  ASSERT_VECTORS_EQUAL( output, MBCartVect(1,0,-1) );
  
  // reflect about arbitrary palne
  MBCartVect norm( 3, 2, 1 );
  ref = MBAffineXform::reflection( norm.array() );
  ASSERT_DOUBLES_EQUAL( ref.matrix().determinant(), -1.0 );
  ref.xform_point( point1.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, point1 );
  ref.xform_point( point2.array(), output.array() );
  ASSERT_DOUBLES_EQUAL( output.length(), point2.length() );
  ASSERT_VECTORS_EQUAL( output, refl(point2,norm) );
  ref.xform_vector( vect1.array(), output.array() );
  ASSERT_DOUBLES_EQUAL( output.length(), vect1.length() );
  ASSERT_VECTORS_EQUAL( output, refl(vect1,norm) );
  ref.xform_vector( vect2.array(), output.array() );
  ASSERT_DOUBLES_EQUAL( output.length(), vect2.length() );
  ASSERT_VECTORS_EQUAL( output, refl(vect2,norm) );
}

void test_scale()
{
  MBCartVect output;
  
  //scale in X only
  MBAffineXform scale = MBAffineXform::scale( MBCartVect(2,1,1).array() );
  scale.xform_point( point1.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, MBCartVect(2*point1[0],point1[1],point1[2]) );
  scale.xform_point( point2.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, MBCartVect(2*point2[0],point2[1],point2[2]) );
  scale.xform_vector( vect1.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, MBCartVect(2*vect1[0],vect1[1],vect1[2]) );
  scale.xform_vector( vect2.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, MBCartVect(2*vect2[0],vect2[1],vect2[2]) );
  
  // scale in all
  scale = MBAffineXform::scale( MBCartVect(0.5,0.5,0.5).array() );
  scale.xform_point( point1.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, 0.5*point1 );
  scale.xform_point( point2.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, 0.5*point2 );
  scale.xform_vector( vect1.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, 0.5*vect1 );
  scale.xform_vector( vect2.array(), output.array() );
  ASSERT_VECTORS_EQUAL( output, 0.5*vect2 );
}

void test_accumulate()
{
  MBCartVect indiv, accum;
  
  // build an group of transforms.  make sure translation is somewhere in the middle
  MBAffineXform move, scal, rot1, rot2, refl;
  move = MBAffineXform::translation( MBCartVect( 5, -5, 1 ).array() );
  scal = MBAffineXform::scale( MBCartVect( 1, 0.5, 2 ).array() );
  rot1 = MBAffineXform::rotation( M_PI/3, MBCartVect( 0.5, 0.5, 1 ).array() );
  rot2 = MBAffineXform::rotation( M_PI/4, MBCartVect( 1.0, 0.0, 0.0 ).array() );
  refl = MBAffineXform::reflection( MBCartVect( -1, -1, 0 ).array() );
  MBAffineXform accu;
  accu.accumulate( scal );
  accu.accumulate( rot1 );
  accu.accumulate( move );
  accu.accumulate( refl );
  accu.accumulate( rot2 );
  
  accu.xform_point( point1.array(), accum.array() );
  scal.xform_point( point1.array(), indiv.array() );
  rot1.xform_point( indiv.array() );
  move.xform_point( indiv.array() );
  refl.xform_point( indiv.array() );
  rot2.xform_point( indiv.array() );
  ASSERT_VECTORS_EQUAL( accum, indiv );
  
  accu.xform_point( point2.array(), accum.array() );
  scal.xform_point( point2.array(), indiv.array() );
  rot1.xform_point( indiv.array() );
  move.xform_point( indiv.array() );
  refl.xform_point( indiv.array() );
  rot2.xform_point( indiv.array() );
  ASSERT_VECTORS_EQUAL( accum, indiv );
  
  accu.xform_vector( vect1.array(), accum.array() );
  scal.xform_vector( vect1.array(), indiv.array() );
  rot1.xform_vector( indiv.array() );
  move.xform_vector( indiv.array() );
  refl.xform_vector( indiv.array() );
  rot2.xform_vector( indiv.array() );
  ASSERT_VECTORS_EQUAL( accum, indiv );
  
  accu.xform_vector( vect2.array(), accum.array() );
  scal.xform_vector( vect2.array(), indiv.array() );
  rot1.xform_vector( indiv.array() );
  move.xform_vector( indiv.array() );
  refl.xform_vector( indiv.array() );
  rot2.xform_vector( indiv.array() );
  ASSERT_VECTORS_EQUAL( accum, indiv );
}

void test_inversion() {
  MBCartVect result;
  
  // build an group of transforms.  make sure translation is somewhere in the middle
  MBAffineXform move, scal, rot1, rot2, refl;
  move = MBAffineXform::translation( MBCartVect( 5, -5, 1 ).array() );
  scal = MBAffineXform::scale( MBCartVect( 1, 0.5, 2 ).array() );
  rot1 = MBAffineXform::rotation( M_PI/3, MBCartVect( 0.5, 0.5, 1 ).array() );
  rot2 = MBAffineXform::rotation( M_PI/4, MBCartVect( 1.0, 0.0, 0.0 ).array() );
  refl = MBAffineXform::reflection( MBCartVect( -1, -1, 0 ).array() );
  MBAffineXform acc;
  acc.accumulate( scal );
  acc.accumulate( rot1 );
  acc.accumulate( move );
  acc.accumulate( refl );
  acc.accumulate( rot2 );
  
  MBAffineXform inv = acc.inverse();
  
  acc.xform_point( point1.array(), result.array() );
  inv.xform_point( result.array() );
  ASSERT_VECTORS_EQUAL( point1, result );
  
  acc.xform_point( point2.array(), result.array() );
  inv.xform_point( result.array() );
  ASSERT_VECTORS_EQUAL( point2, result );
  
  acc.xform_vector( vect1.array(), result.array() );
  inv.xform_vector( result.array() );
  ASSERT_VECTORS_EQUAL( vect1, result );
  
  acc.xform_vector( vect2.array(), result.array() );
  inv.xform_vector( result.array() );
  ASSERT_VECTORS_EQUAL( vect2, result );
}
  
void test_is_reflection()
{
  MBAffineXform refl1, refl2, scale;
  refl1 = MBAffineXform::reflection( MBCartVect( -1, -1, 0).array() );
  refl2 = MBAffineXform::reflection( MBCartVect(  1,  0, 0).array() );
  scale = MBAffineXform::scale( MBCartVect( -1, 1, 1 ).array() );
  
  ASSERT( refl1.reflection() );
  ASSERT( refl2.reflection() );
  ASSERT( scale.reflection() );
  
  MBAffineXform inv1, inv2, inv3;
  inv1 = refl1.inverse();
  inv2 = refl2.inverse();
  inv3 = scale.inverse();
  
  ASSERT( inv1.reflection() );
  ASSERT( inv2.reflection() );
  ASSERT( inv3.reflection() );
  
  refl1.accumulate( refl2 );
  refl2.accumulate( scale );
  ASSERT( ! refl1.reflection() );
  ASSERT( ! refl2.reflection() );
  
  MBAffineXform rot, mov;
  rot = MBAffineXform::rotation( M_PI/4, MBCartVect(1,1,1).array() );
  mov = MBAffineXform::translation( MBCartVect(-5,6,7).array() );
  ASSERT( !rot.reflection() );
  ASSERT( !mov.reflection() );
}

int main()
{
  test_none();
  test_translation();
  test_rotation();
  test_reflection();
  test_scale();
  test_accumulate();
  test_inversion();
  test_is_reflection();
  return error_count;
}

#endif  // #ifdef TEST
