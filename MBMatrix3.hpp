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

/**\file MBMatrix3.hpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2006-07-18
 */

#ifndef MB_MATRIX_3_HPP
#define MB_MATRIX_3_HPP

#include "MBCartVect.hpp"
#include "MBInterface.hpp"
#include <iosfwd>

class MBMatrix3 
{
  double d[9];
public:

  inline MBMatrix3() 
    {}
  
  inline MBMatrix3( double diagonal )
    { 
      d[0] = d[4] = d[8] = diagonal;
      d[1] = d[2] = d[3] = 0;
      d[5] = d[6] = d[7] = 0;
    }
  
  inline MBMatrix3( const MBCartVect& row0,
                    const MBCartVect& row1,
                    const MBCartVect& row2 )
    {
      row0.get( d );
      row1.get( d+3 );
      row2.get( d+6 );
    }
  
  inline MBMatrix3( const double* v )
    { 
      d[0] = v[0]; d[1] = v[1]; d[2] = v[2];
      d[3] = v[3]; d[4] = v[4]; d[5] = v[5]; 
      d[6] = v[6]; d[7] = v[7]; d[8] = v[8];
    }
  
  inline MBMatrix3( double v00, double v01, double v02,
                    double v10, double v11, double v12,
                    double v20, double v21, double v22 )
    {
      d[0] = v00; d[1] = v01; d[2] = v02;
      d[3] = v10; d[4] = v11; d[5] = v12;
      d[6] = v20; d[7] = v21; d[8] = v22;
    }
  
  inline MBMatrix3( const MBMatrix3& m )
    {
      d[0] = m.d[0]; d[1] = m.d[1]; d[2] = m.d[2];
      d[3] = m.d[3]; d[4] = m.d[4]; d[5] = m.d[5];
      d[6] = m.d[6]; d[7] = m.d[7]; d[8] = m.d[8];
    }
   
  inline MBMatrix3& operator=( const MBMatrix3& m )
    {
      d[0] = m.d[0]; d[1] = m.d[1]; d[2] = m.d[2];
      d[3] = m.d[3]; d[4] = m.d[4]; d[5] = m.d[5];
      d[6] = m.d[6]; d[7] = m.d[7]; d[8] = m.d[8];
      return *this;
    }
  
  inline MBMatrix3& operator=( const double* v )
    { 
      d[0] = v[0]; d[1] = v[1]; d[2] = v[2];
      d[3] = v[3]; d[4] = v[4]; d[5] = v[5]; 
      d[6] = v[6]; d[7] = v[7]; d[8] = v[8];
      return *this;
    }

  inline double* operator[]( unsigned i )
    { return d + 3*i; }
  inline const double* operator[]( unsigned i ) const
    { return d + 3*i; }
  inline double& operator()(unsigned r, unsigned c)
    { return d[3*r+c]; }
  inline double operator()(unsigned r, unsigned c) const
    { return d[3*r+c]; }
  
  inline MBMatrix3& operator+=( const MBMatrix3& m )
    {
      d[0] += m.d[0]; d[1] += m.d[1]; d[2] += m.d[2];
      d[3] += m.d[3]; d[4] += m.d[4]; d[5] += m.d[5];
      d[6] += m.d[6]; d[7] += m.d[7]; d[8] += m.d[8];
      return *this;
    }
  
  inline MBMatrix3& operator-=( const MBMatrix3& m )
    {
      d[0] -= m.d[0]; d[1] -= m.d[1]; d[2] -= m.d[2];
      d[3] -= m.d[3]; d[4] -= m.d[4]; d[5] -= m.d[5];
      d[6] -= m.d[6]; d[7] -= m.d[7]; d[8] -= m.d[8];
      return *this;
    }
  
  inline MBMatrix3& operator*=( double s )
    {
      d[0] *= s; d[1] *= s; d[2] *= s;
      d[3] *= s; d[4] *= s; d[5] *= s;
      d[6] *= s; d[7] *= s; d[8] *= s;
      return *this;
    }
  
  inline MBMatrix3& operator/=( double s )
    {
      d[0] /= s; d[1] /= s; d[2] /= s;
      d[3] /= s; d[4] /= s; d[5] /= s;
      d[6] /= s; d[7] /= s; d[8] /= s;
      return *this;
    }
  
  inline MBMatrix3& operator*=( const MBMatrix3& m );
  
  inline double determinant() const;
  
  inline MBMatrix3 inverse() const;
  
  inline MBMatrix3 transpose() const;
  
};

inline MBMatrix3 operator+( const MBMatrix3& a, const MBMatrix3& b )
  { return MBMatrix3(a) += b; }

inline MBMatrix3 operator-( const MBMatrix3& a, const MBMatrix3& b )
  { return MBMatrix3(a) -= b; }

inline MBMatrix3 operator*( const MBMatrix3& a, const MBMatrix3& b )
{
  return MBMatrix3( a(0,0) * b(0,0) + a(0,1) * b(1,0) + a(0,2) * b(2,0),
                    a(0,0) * b(0,1) + a(0,1) * b(1,1) + a(0,2) * b(2,1),
                    a(0,0) * b(0,2) + a(0,1) * b(1,2) + a(0,2) * b(2,2),
                    a(1,0) * b(0,0) + a(1,1) * b(1,0) + a(1,2) * b(2,0),
                    a(1,0) * b(0,1) + a(1,1) * b(1,1) + a(1,2) * b(2,1),
                    a(1,0) * b(0,2) + a(1,1) * b(1,2) + a(1,2) * b(2,2),
                    a(2,0) * b(0,0) + a(2,1) * b(1,0) + a(2,2) * b(2,0),
                    a(2,0) * b(0,1) + a(2,1) * b(1,1) + a(2,2) * b(2,1),
                    a(2,0) * b(0,2) + a(2,1) * b(1,2) + a(2,2) * b(2,2) );
}

inline MBMatrix3& MBMatrix3::operator*=( const MBMatrix3& m )
  { return *this = MBMatrix3(*this) * m; }

inline MBMatrix3 outer_product( const MBCartVect& u,
                                const MBCartVect& v )
{
  return MBMatrix3( u[0] * v[0], u[0] * v[1], u[0] * v[2],
                    u[1] * v[0], u[1] * v[1], u[1] * v[2],
                    u[2] * v[0], u[2] * v[1], u[2] * v[2] );
}

inline double MBMatrix3::determinant() const
{
  return d[0] * d[4] * d[8] 
       + d[1] * d[5] * d[6]
       + d[2] * d[3] * d[7]
       - d[0] * d[5] * d[7]
       - d[1] * d[3] * d[8]
       - d[2] * d[4] * d[6];
}

inline MBMatrix3 MBMatrix3::inverse() const
{
  double i = 1.0 / determinant();
  return MBMatrix3( i * (d[4] * d[8] - d[5] * d[7]),
                    i * (d[2] * d[7] - d[8] * d[1]),
                    i * (d[1] * d[5] - d[4] * d[2]),
                    i * (d[5] * d[6] - d[8] * d[3]),
                    i * (d[0] * d[8] - d[6] * d[2]),
                    i * (d[2] * d[3] - d[5] * d[0]),
                    i * (d[3] * d[7] - d[6] * d[4]),
                    i * (d[1] * d[6] - d[7] * d[0]),
                    i * (d[0] * d[4] - d[3] * d[1]) );
}

inline MBMatrix3 MBMatrix3::transpose() const
{
  return MBMatrix3( d[0], d[3], d[6],
                    d[1], d[4], d[7],
                    d[2], d[5], d[8] );
}

                         
MBErrorCode EigenDecomp( const MBMatrix3& a, 
                         double Eigenvalues[3],
                         MBCartVect Eigenvectors[3] );

std::ostream& operator<<( std::ostream&, const MBMatrix3& );

#endif
