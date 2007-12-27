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

/** \class MBEdgeSizeEvaluator
  *
  * This is an abstract class that embodies the rule used during edge-based mesh
  * refinement to decide whether an edge should be subdivided or not.
  * Subclasses must implement the pure virtual evaluate_edge() function.
  *
  * \author David Thompson
  *
  * \date 19 November 2007
  */
#ifndef MB_EDGESIZEEVALUATOR_H
#define MB_EDGESIZEEVALUATOR_H

#include "MBTypes.h" // for MB_DLL_EXPORT

#include <vector>

class MBInterface;

class MB_DLL_EXPORT MBEdgeSizeEvaluator
{
public:
  MBEdgeSizeEvaluator( MBInterface* );
  virtual ~MBEdgeSizeEvaluator();

  virtual bool evaluate_edge(
    const double* p0, const void* t0,
    double* p1, void* t1,
    const double* p2, const void* t2 ) = 0;

  void reset_vertex_tags();
  int add_vertex_tag( MBTag tag_handle );
  int get_vertex_tag_size() { return this->vertex_size; }
  void evaluate_tags_at_midpoint(
    const double* c0, const void* t0,
    const double* cm, void* tm,
    const double* c1, const void* t1 ) const;

protected:
  std::vector< std::pair< MBTag, int > > vertex_tags;
  int vertex_size;
  MBInterface* mesh;
};

#endif // MB_EDGESIZEEVALUATOR_H
