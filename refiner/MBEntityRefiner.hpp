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

/** \class MBEntityRefiner
  *
  * This is an abstract class that contains the method used for per-entity refinement.
  * Subclasses must implement the pure virtual refine_entity() function and
  * may implement the vertices_per_split() function.
  * This class constructor requires a non-NULL pointer to a mesh so that, given an
  * entity handle, it can look up vertex coordinates and tags to prepare arguments for
  * the refine_entity() method.
  *
  * Although the MBMeshRefiner class may not initially support it, entity refiners
  * are required to support some level of recursion.
  * The maximum number of recursive calls allowed may be set with
  * MBEntityRefiner::set_maximum_number_of_subdivisions().
  * As a convenience, some of the framework for recursion is provided by the
  * MBEntityRefiner class.
  *
  * Specifically, MBEntityRefiner stores a pair of heap arrays
  * to hold edge midpoint vertex coordinates and tag values pre-allocated to the
  * maximum recursion depth so that no repeated allocation and deallocation
  * needs to take place during refinement.
  * To use these heaps, subclasses should call reset_heap_pointers() upon entry to
  * MBEntityRefiner::refine_entity().
  * Then, when the edge size evaluator requires an edge to be split, subclasses
  * should call heap_coord_storage() and heap_tag_storage() to obtain pointers as
  * required.
  *
  * \author David Thompson
  * \author Philippe Pebay
  *
  * \date 24 December 2007
  */
/** \class MBEntityRefinerOutputFunctor
  *
  * This is an abstract class used by MBEntityRefiner to output entities that are the product of refinement.
  * The parenthesis operator is overloaded with two forms:
  * one used for appending a vertex to an entity,
  * the other used to finalize the creation of the entity by specifying its type.
  *
  * \author David Thompson
  * \author Philippe Pebay
  *
  * \date 26 December 2007
  */
#ifndef MB_ENTITYREFINER_H
#define MB_ENTITYREFINER_H

#include "MBTypes.h" // for MB_DLL_EXPORT

#include <vector>

class MBInterface;
class MBEdgeSizeEvaluator;

class MB_DLL_EXPORT MBEntityRefinerOutputFunctor
{
public:
  virtual void operator () ( const double* vcoords, const void* vtags ) = 0;
  virtual void operator () ( MBEntityType etyp ) = 0;
};

class MB_DLL_EXPORT MBEntityRefiner
{
public:
  MBEntityRefiner( MBInterface* );
  virtual ~MBEntityRefiner();

  virtual bool refine_entity( MBEntityHandle ) = 0;
  virtual unsigned long get_heap_size_bound( int max_recursions ) const = 0;

  virtual bool set_edge_size_evaluator( MBEdgeSizeEvaluator* );
  MBEdgeSizeEvaluator* get_edge_size_evaluator() { return this->edge_size_evaluator; }

  virtual bool set_output_functor( MBEntityRefinerOutputFunctor* func_obj );
  MBEntityRefinerOutputFunctor* get_output_functor() { return this->output_functor; }

  virtual bool set_minimum_number_of_subdivisions( int mn );
  int get_minimum_number_of_subdivisions() const { return this->minimum_number_of_subdivisions; }

  virtual bool set_maximum_number_of_subdivisions( int mx );
  int get_maximum_number_of_subdivisions() const { return this->maximum_number_of_subdivisions; }

protected:
  MBInterface* mesh;
  MBEdgeSizeEvaluator* edge_size_evaluator;
  MBEntityRefinerOutputFunctor* output_functor;
  int minimum_number_of_subdivisions;
  int maximum_number_of_subdivisions;
  std::vector<double> coord_heap;
  std::vector<double>::iterator current_coord;
  std::vector<char> tag_heap;
  std::vector<char>::iterator current_tag;

  void update_heap_size();
  void reset_heap_pointers();
  double* heap_coord_storage();
  void* heap_tag_storage();
  void evaluate_tags_at_midpoint( const double* c0, const void* t0, double* cm, void* tm, const double* c1, const void* t1 ) const;
};

#endif // MB_ENTITYREFINER_H
