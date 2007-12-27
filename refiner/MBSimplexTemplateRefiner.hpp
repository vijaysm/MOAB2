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
  static int template_index[64][2];
  static int permutations_from_index[24][14];
  static int templates[];

  void refine_0_simplex( const double* v0, const void* t0 );
  bool refine_1_simplex( int max_depth,
    const double* v0, const void* t0, const double* v1, const void* t1 );
  bool refine_2_simplex( int max_depth, int move,
    const double* v0, const void* t0, const double* v1, const void* t1, const double* v2, const void* t2 );
  bool refine_3_simplex( int max_depth,
                         double* v0, void* t0, 
                         double* v1, void* t1, 
                         double* v2, void* t2,
                         double* v3, void* t3 );
  static bool compare_Hopf_cross_string_dist( const double* v00, const double* v01, const double* v10, const double* v11 );
  void evaluate_tags_at_facepoint( const double* c0, const void* t0,
				   const double* c1, const void* t1,
				   const double* c2, const void* t2,
				   const double* cm, void* tm ) const;

  int best_tets( int* alternates, double*[14], int, int ) { return alternates[0]; }
};
#endif // MB_SIMPLEXTEMPLATEREFINER_H

