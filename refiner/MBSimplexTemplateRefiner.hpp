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
  *
  * \author David Thompson
  * \author Philippe Pebay
  *
  * \date 24 December 2007
  */
#ifndef MB_SIMPLEXTEMPLATEREFINER_H
#define MB_SIMPLEXTEMPLATEREFINER_H

#include "MBEntityRefiner.h"

class MB_DLL_EXPORT MBSimplexTemplateRefiner : public MBEntityRefiner
{
public:
  /// Construct a template refiner.
  MBSimplexTemplateRefiner();
  /// Destruction is virtual so subclasses may clean up after refinement.
  virtual ~MBSimplexTemplateRefiner();

  virtual bool refine_entity( MBEntityHandle );

protected:
  static int* template_index;
  static int* templates;
};
#endif // MB_SIMPLEXTEMPLATEREFINER_H

