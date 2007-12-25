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
  * This is an abstract class that contains the method used for per-entity
  * refinement
  * Subclasses must implement the pure virtual refine_entity() function.
  *
  * \author David Thompson
  * \author Philippe Pebay
  *
  * \date 19 November 2007
  */
#ifndef MB_ENTITYREFINER_H
#define MB_ENTITYREFINER_H

#include "MBTypes.h" // for MB_DLL_EXPORT

#include <vector>

class MBInterface;
class MBEdgeSizeEvaluator;

class MB_DLL_EXPORT MBEntityRefiner
{
public:
  /// Construct an entity refiner.
  MBEntityRefiner( MBInterface* );
  /// Destruction is virtual so subclasses may clean up after refinement.
  virtual ~MBEntityRefiner();

  virtual bool refine_entity( MBEntityHandle ) = 0;

  bool set_edge_size_evaluator( MBEdgeSizeEvaluator* );
  MBEdgeSizeEvaluator* get_edge_size_evaluator() { return this->edge_size_evaluator; };

protected:
  MBInterface* mesh;
  MBEdgeSizeEvaluator* edge_size_evaluator;
};

#endif // MB_ENTITYREFINER_H
