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

/** \class MBMeshRefiner
  *
  * This is an class that contains the method used for mesh refinement.
  *
  * \author Philippe Pebay
  * \author David Thompson
  *
  * \date 19 November 2007
  */
#ifndef MB_MESHREFINER_H
#define MB_MESHREFINER_H

#include "MBTypes.h" // for MB_DLL_EXPORT

#include <vector>

class MBInterface;
class MBEntityRefiner;

class MB_DLL_EXPORT MBMeshRefiner
{
public:
  /// Construct a mesh refiner.
  MBMeshRefiner( MBInterface* );
  /// Destruction is virtual so subclasses may clean up after refinement.
  virtual ~MBMeshRefiner();

  virtual bool refine_mesh();

  bool set_entity_refiner( MBEntityRefiner* );
  MBEntityRefiner* get_entity_refiner() { return this->entity_refiner; };

protected:
  MBInterface* mesh;
  MBEntityRefiner* entity_refiner;
};

#endif // MB_MESHREFINER_H
