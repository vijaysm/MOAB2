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

/**\class MBMeshOutputFunctor
  *\brief Implements the abstract MBEntityRefinerOutputFunctor class.
  *
  * This class is a concrete implementation of the MBEntityRefinerOutputFunctor.
  * It creates new vertices and regions in a new or existing mesh as
  * the input entities are streamed through the refiner.
  *
  * \author David Thompson
  * \author Philippe Pebay
  *
  * \date 28 July 2008
  */
#ifndef MB_MESHOUTPUTFUNCTOR_HPP
#define MB_MESHOUTPUTFUNCTOR_HPP

#include "MBTypes.h"
#include "MBEntityRefiner.hpp"
#include "MBProcessSet.hpp"

#include <vector>
#include <map>

#include <string.h>

class MBSplitVerticesBase;
class MBParallelComm;

class MBMeshOutputFunctor : public MBEntityRefinerOutputFunctor
{
public:
  MBMeshOutputFunctor( MBRefinerTagManager* tag_mgr );
  ~MBMeshOutputFunctor();

  void print_vert_crud( MBEntityHandle vout, int nvhash, MBEntityHandle* vhash, const double* vcoords, const void* vtags );
  void assign_global_ids( MBParallelComm* comm );

  void assign_tags( MBEntityHandle vhandle, const void* vtags );

  virtual MBEntityHandle operator () ( MBEntityHandle vhash, const double* vcoords, const void* vtags );
  virtual MBEntityHandle operator () ( int nvhash, MBEntityHandle* vhash, const double* vcoords, const void* vtags );
  virtual void operator () ( MBEntityHandle h );
  virtual void operator () ( MBEntityType etyp );

  MBInterface* mesh_in;
  MBInterface* mesh_out;
  bool input_is_output;
  std::vector<MBSplitVerticesBase*> split_vertices;
  std::vector<MBSplitVerticesBase*> new_entities;
  std::vector<MBEntityHandle> elem_vert;
  MBRefinerTagManager* tag_manager;
  MBEntityHandle destination_set;
  std::map<MBProcessSet,int> proc_partition_counts;
};

#endif // MB_MESHOUTPUTFUNCTOR_HPP
