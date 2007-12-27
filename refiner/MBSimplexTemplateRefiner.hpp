/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2007 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

/** \class MBSimplexTemplateRefiner
  *
  * This is a concrete subclass of MBEntityRefiner that implements
  * refinement using templates applied to simplices.
  * Entities that are not simplices are divided into tetrahedra,
  * triangles, or lines before being processed.
  * Points are passed through unchanged.
  *
  * \author David Thompson
  * \author Philippe Pebay
  *
  * \date 24 December 2007
  */
#ifndef MB_SIMPLEXTEMPLATEREFINER_H
#define MB_SIMPLEXTEMPLATEREFINER_H

#include "MBEntityRefiner.hpp"

class MB_DLL_EXPORT MBSimplexTemplateRefiner : public MBEntityRefiner
{
public:
  MBSimplexTemplateRefiner( MBInterface* mesh );
  virtual ~MBSimplexTemplateRefiner();

  virtual bool refine_entity( MBEntityHandle entity );
  virtual unsigned long get_heap_size_bound( int max_recursions ) const { return 48 * 4 * ( 1 << max_recursions ); }

protected:
  static int* template_index;
  static int* templates;

  void refine_0_simplex( double* v0, const void* t0 );
  bool refine_1_simplex( int max_depth, double* v0, void* t0, double* v1, void* t1 );
  bool refine_2_simplex();
  bool refine_3_simplex();
};
#endif // MB_SIMPLEXTEMPLATEREFINER_H

