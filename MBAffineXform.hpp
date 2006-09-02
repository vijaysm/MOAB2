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

#ifndef MB_AFFINE_XFORM_HPP
#define MB_AFFINE_XFORM_HPP

#include "MBCartVect.hpp"
#include "MBMatrix3.hpp"

class MBAffineXform
{
  public:
  
    inline MBAffineXform();
    
    inline MBAffineXform( const double* three_by_three, 
                          const double* translation );
                          
    inline MBAffineXform( const MBMatrix3& mat, const MBCartVect& off );
    
    /** move */
    static inline MBAffineXform translation( const double* vector );
    /** rotate about axis through origin */
    static inline MBAffineXform rotation( double radians, const double* axis );
    /** reflect about plane through origin */
    static inline MBAffineXform reflection( const double* plane_normal );
    /** scale about origin */
    static inline MBAffineXform scale( const double* fractions );
    
    inline void accumulate( const MBAffineXform& other );
    
    inline void xform_point( const double* input, double* output ) const;
    inline void xform_point( double* in_out ) const;
    
    inline void xform_vector( const double* input, double* output ) const;
    inline void xform_vector( double* in_out ) const;
    
    MBAffineXform inverse() const;
    
    static MBErrorCode get_tag( MBTag& tag_handle_out,
                                MBInterface* moab,
                                const char* tagname = 0 );
    
    const MBMatrix3& matrix() const { return mMatrix; }
    const MBCartVect& offset() const { return mOffset; }

  private:
  
    MBMatrix3 mMatrix;
    MBCartVect mOffset;
};

inline MBAffineXform::MBAffineXform()
  : mMatrix(1.0), mOffset(0.0) 
  {}

inline MBAffineXform::MBAffineXform( const double* three_by_three, 
                                     const double* translation )
 : mMatrix(three_by_three), mOffset(translation)
 {}

inline MBAffineXform::MBAffineXform( const MBMatrix3& mat, const MBCartVect& off )
  : mMatrix(mat), mOffset(off)
  {}

inline MBAffineXform MBAffineXform::translation( const double* vector )
{
  return MBAffineXform( MBMatrix3(1.0), MBCartVect(vector) );
}

inline MBAffineXform MBAffineXform::rotation( double angle, const double* axis )
{
  MBCartVect a(axis);
  a.normalize();

  const double c = cos(angle);
  const double s = sin(angle);
  const MBMatrix3 m1(    c,   -a[2]*s, a[1]*s,
                       a[2]*s,   c,   -a[0]*s,
                      -a[1]*s, a[0]*s,   c    );
  return MBAffineXform( m1 + (1.0-c)*outer_product( a, a ), MBCartVect(0.0) );
}

inline MBAffineXform MBAffineXform::reflection( const double* plane_normal )
{
  double i = plane_normal[0];
  double j = plane_normal[1];
  double k = plane_normal[2];
  MBMatrix3 m( j*j+k*k-i*i,    -2.0*i*j,    -2.0*i*k,
                  -2.0*i*j, i*i+k*k-j*j,    -2.0*j*k,
                  -2.0*i*k,    -2.0*j*k, i*i+j*j-k*k );
  m *= 1.0 / (i*i + j*j + k*k); //normalize
  return MBAffineXform( m, MBCartVect(0.0) );
}

inline MBAffineXform MBAffineXform::scale( const double* f )
{
  return MBAffineXform( MBMatrix3( MBCartVect(f) ), MBCartVect( 0.0 ) );
}

inline void MBAffineXform::accumulate( const MBAffineXform& other )
{
  mMatrix = other.mMatrix * mMatrix;
  other.xform_point( mOffset.array() );
}

inline void MBAffineXform::xform_point( const double* input, double* output ) const
{
  xform_vector( input, output );
  output[0] += mOffset[0];
  output[1] += mOffset[1];
  output[2] += mOffset[2];
}

inline void MBAffineXform::xform_point( double* in_out ) const
{
  xform_vector( in_out );
  in_out[0] += mOffset[0];
  in_out[1] += mOffset[1];
  in_out[2] += mOffset[2];
} 

inline void MBAffineXform::xform_vector( const double* input, double* output ) const
{
  output[0] = input[0]*mMatrix[0][0] + input[1]*mMatrix[0][1] + input[2]*mMatrix[0][2];
  output[1] = input[0]*mMatrix[1][0] + input[1]*mMatrix[1][1] + input[2]*mMatrix[1][2];
  output[2] = input[0]*mMatrix[2][0] + input[1]*mMatrix[2][1] + input[2]*mMatrix[2][2];
}

inline void MBAffineXform::xform_vector( double* in_out ) const
{
  double input[] = { in_out[0], in_out[1], in_out[2] };
  xform_vector( input, in_out );
}

inline MBAffineXform MBAffineXform::inverse() const
{
  MBMatrix3 m = mMatrix.inverse();
  return MBAffineXform( m, m * -mOffset );
}


#endif
