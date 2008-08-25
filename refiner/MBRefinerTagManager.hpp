/*
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

/**\class MBRefinerTagManager
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

#include "MBProcessSet.hpp"

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

  MBTag input_parallel_status() { return this->tag_ipstatus; }
  MBTag input_shared_proc() { return this->tag_ipsproc; }
  MBTag input_shared_procs() { return this->tag_ipsprocs; }

  int get_input_gids( int n, const MBEntityHandle* ents, std::vector<int>& gids );
  int get_output_gids( int n, const MBEntityHandle* ents, std::vector<int>& gids );
  int set_gid( MBEntityHandle ent, int gid );

  void set_sharing( MBEntityHandle ent_handle, MBProcessSet& procs );
  void get_common_processes( int num, const MBEntityHandle* src, MBProcessSet& common_shared_procs, bool on_output_mesh = true );

protected:
  std::vector< std::pair< MBTag, int > > input_vertex_tags;
  std::vector< std::pair< MBTag, int > > output_vertex_tags;
  int vertex_size;
  MBInterface* input_mesh;
  MBInterface* output_mesh;
  MBTag tag_ipstatus; // Handle for PARALLEL_STATUS on mesh_in
  MBTag tag_ipsprocs; // Handle for PARALLEL_SHARED_PROCS on mesh_in
  MBTag tag_ipsproc;  // Handle for PARALLEL_SHARED_PROC on mesh_in
  MBTag tag_ipshands; // Handle for PARALLEL_SHARED_HANDLES on mesh_in
  MBTag tag_ipshand;  // Handle for PARALLEL_SHARED_HANDLE on mesh_in
  MBTag tag_igid;     // Handle for global IDs on mesh_in
  MBTag tag_opstatus; // Handle for PARALLEL_STATUS on mesh_out
  MBTag tag_opsprocs; // Handle for PARALLEL_SHARED_PROCS on mesh_out
  MBTag tag_opsproc;  // Handle for PARALLEL_SHARED_PROC on mesh_out
  MBTag tag_opshands; // Handle for PARALLEL_SHARED_HANDLES on mesh_out
  MBTag tag_opshand;  // Handle for PARALLEL_SHARED_HANDLE on mesh_out
  MBTag tag_ogid;     // Handle for global IDs on mesh_out
  int rank;
  std::vector<int> shared_procs_in; // Used to hold procs sharing an input vert.
  std::vector<int> shared_procs_out; // Used to hold procs sharing an output entity.
  MBProcessSet current_shared_procs; // Holds process list as it is being accumulated
};

#endif // MB_REFINERTAGMANAGER_H
