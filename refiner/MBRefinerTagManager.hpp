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

/** \class MBRefinerTagManager
  *
  * This a class that manages which tags an edge refiner should include
  * on output vertices created during mesh refinement.
  * The 
  *
  * \author David Thompson
  *
  * \date 12 June 2008
  */
#ifndef MB_REFINERTAGMANAGER_H
#define MB_REFINERTAGMANAGER_H

#include "MBTypes.h" // for MB_DLL_EXPORT

#include <vector>

class MBInterface;

class MB_DLL_EXPORT MBRefinerTagManager
{
public:
  MBRefinerTagManager( MBInterface* in_mesh, MBInterface* out_mesh );
  virtual ~MBRefinerTagManager();

  void reset_vertex_tags();
  int add_vertex_tag( MBTag tag_handle );
  int get_vertex_tag_size() const { return this->vertex_size; }
  int get_number_of_vertex_tags() const { return this->input_vertex_tags.size(); }

  void create_output_tags();

  void get_input_vertex_tag( int i, MBTag& tag, int& byte_offset );
  void get_output_vertex_tag( int i, MBTag& tag, int& byte_offset );

  MBInterface* get_input_mesh() { return this->input_mesh; }
  MBInterface* get_output_mesh() { return this->output_mesh; }

  MBTag shared_proc() { return this->tag_psproc; }
  MBTag shared_procs() { return this->tag_psprocs; }

protected:
  std::vector< std::pair< MBTag, int > > input_vertex_tags;
  std::vector< std::pair< MBTag, int > > output_vertex_tags;
  int vertex_size;
  MBInterface* input_mesh;
  MBInterface* output_mesh;
  MBTag tag_pstatus; // Handle for PARALLEL_STATUS on mesh_in
  MBTag tag_psprocs; // Handle for PARALLEL_SHARED_PROCS on mesh_in
  MBTag tag_psproc;  // Handle for PARALLEL_SHARED_PROC on mesh_in
  MBTag tag_pshands; // Handle for PARALLEL_SHARED_HANDLES on mesh_in
  MBTag tag_pshand;  // Handle for PARALLEL_SHARED_HANDLE on mesh_in
};

#endif // MB_REFINERTAGMANAGER_H
