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

#include "MBTypes.h" // for MB_DLL_EXPORT

class MBRefinerTagManager;
class MBSimplexTemplateRefiner;

class MB_DLL_EXPORT MBSimplexTemplateTagAssigner
{
public:
  MBSimplexTemplateTagAssigner( MBSimplexTemplateRefiner* );
  virtual ~MBSimplexTemplateTagAssigner();
  
  virtual void operator () ( const double* c0, const void* t0, MBEntityHandle h0,
                             const double* cm, void* tm,
                             const double* c1, const void* t1, MBEntityHandle h1 );
  virtual void operator () ( const void* t0,
                             const void* t1,
                             const void* t2,
                             void* tp );
  virtual void set_tag_manager( MBRefinerTagManager* tmgr );

protected:
  MBSimplexTemplateRefiner* mesh_refiner;
  MBRefinerTagManager* tag_manager;
};

#endif // MB_SIMPLEXTEMPLATETAGASSIGNER_H
