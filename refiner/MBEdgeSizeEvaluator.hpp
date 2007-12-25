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
  /// Construct an evaluator.
  MBEdgeSizeEvaluator( MBInterface* );
  /// Destruction is virtual so subclasses may clean up after refinement.
  virtual ~MBEdgeSizeEvaluator();

  /** \brief Returns true if the edge \a p0 - \a p2 should be subdivided, false otherwise.
    *
    * The arguments \a p0, \a p1, and \a p2 are all pointers to arrays of 6 doubles each
    * while the arguments \a t0, \a t1, and \a t2 are all pointers to arrays of tag data
    * defined at the corresponding point. While the endpoints \a p0 and \a p2 are
    * immutable, the mid-edge point coordinates \a p1 and tag data \a t1 may be altered by
    * evaluate_edge(). Altered values will be ignored if evaluate_edge() returns false.
    * Be careful to ensure that all calls to evaluate_edge() perform identical modifications
    * given identical input values!
    *
    * A list of tags passed in \a t0, \a t1, and \a t2 is stored in the vertexTags member.
    * The vertexSize member stores the total length of data associated with each pointer (in bytes).
    * Subclasses may access vertexTags and vertexSize directly; the refiner uses public methods to
    * populate vertexTags before evaluate_edge() is called.
    */
  virtual bool evaluate_edge(
    const double* p0, const void* t0,
    double* p1, void* t1,
    const double* p2, const void* t2 ) = 0;

  /// Clear the list of tag values that will appear past the vertex coordinates in \a p0, \a p1, and \a p2.
  void reset_vertex_tags();
  /** Add a tag to the list of tag values that will appear past the vertex coordinates.
    * The return value is the offset into each vertex coordinate pointer (\a p0, \a p1, \a p2) where the
    * tag value(s) will be stored.
    */
  int add_vertex_tag( MBTag tag_handle );
  /// Return the number of bytes to allocate for tag data per point.
  int get_vertex_tag_size() { return this->vertexSize; }

protected:
  std::vector< std::pair< MBTag, int > > vertexTags;
  int vertexSize;
  MBInterface* mesh;
};

#endif // MB_EDGESIZEEVALUATOR_H
