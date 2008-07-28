/*
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

/**\class MBMeshRefiner
  *\brief Refine a mesh using a streaming operation.
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
#include "MBRange.hpp"

#include <vector>

class MBInterface;
class MBEntityRefiner;
class MBParallelComm;
class MBRefinerTagManager;
class MBMeshOutputFunctor;

class MB_DLL_EXPORT MBMeshRefiner
{
public:
  MBMeshRefiner( MBInterface* imesh, MBInterface* omesh );
  virtual ~MBMeshRefiner();

  bool set_entity_refiner( MBEntityRefiner* );
  MBEntityRefiner* get_entity_refiner() { return this->entity_refiner; }

  bool set_comm( MBParallelComm* c ) { if ( ! c || this->comm == c ) return false; this->comm = c; return true; }
  MBParallelComm* get_comm() { return this->comm; }

  MBRefinerTagManager* get_tag_manager() { return this->tag_manager; }
  const MBRefinerTagManager* get_tag_manager() const { return this->tag_manager; }
  void reset_vertex_tags();
  int add_vertex_tag( MBTag tag_handle );

  virtual bool refine( MBRange& );

protected:
  MBInterface* mesh_in;
  MBInterface* mesh_out;
  MBEntityRefiner* entity_refiner;
  MBRefinerTagManager* tag_manager;
  MBMeshOutputFunctor* output_functor;
  MBParallelComm* comm;
};

#endif // MB_MESHREFINER_H
