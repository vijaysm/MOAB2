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

/** \class MBSimplexTemplateTagAssigner
  *
  * This is an class that embodies the process of assigning tag
  * values to new vertices based on some pre-existing neighbors in a 
  * simplicial mesh.
  *
  * \author David Thompson
  * \author Philippe Pebay
  *
  * \date 28 December 2007
  */
#ifndef MB_SIMPLEXTEMPLATETAGASSIGNER_H
#define MB_SIMPLEXTEMPLATETAGASSIGNER_H

#include "MBMeshRefiner.hpp"

class MB_DLL_EXPORT MBSimplexTemplateTagAssigner
{
public:
  MBSimplexTemplateTagAssigner( MBMeshRefiner* );
  virtual ~MBSimplexTemplateTagAssigner();

  virtual void operator()( const void* ta, const void* tb, void* tp );
  virtual void operator()( const void* ta, const void* tb, const void* tc, void* tp );

protected:
  MBMeshRefiner* mesh_refiner;
};

#endif // MB_SIMPLEXTEMPLATETAGASSIGNER_H
