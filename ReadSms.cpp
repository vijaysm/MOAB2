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

/**
 * \class ReadSms
 * \brief Sms (http://www.geuz.org/sms) file reader
 * \author Jason Kraftcheck
 */

#include "ReadSms.hpp"
#include "FileTokenizer.hpp" // for file tokenizer
#include "MBInternals.hpp"
#include "MBInterface.hpp"
#include "MBReadUtilIface.hpp"
#include "MBRange.hpp"
#include "MBTagConventions.hpp"
#include "MBParallelConventions.h"
#include "MBCN.hpp"

#include <errno.h>
#include <string.h>
#include <map>
#include <set>
#include <iostream>

#define CHECK(a) if (MB_SUCCESS != result) {      \
      std::cerr << a << std::endl;                \
      return result;                              \
    }
    

MBReaderIface* ReadSms::factory( MBInterface* iface )
  { return new ReadSms(iface); }

ReadSms::ReadSms(MBInterface* impl)
    : mdbImpl(impl)
{
  void* ptr = 0;
  mdbImpl->query_interface("MBReadUtilIface", &ptr);
  readMeshIface = reinterpret_cast<MBReadUtilIface*>(ptr);
}

ReadSms::~ReadSms()
{
  if (readMeshIface) {
    mdbImpl->release_interface("MBReadUtilIface", readMeshIface);
    readMeshIface = 0;
  }
}

MBErrorCode ReadSms::load_file( const char* filename, 
                                MBEntityHandle& file_set,
                                const FileOptions& ,
                                const char* name,
                                const int*, const int )
{
  if (name) {
    readMeshIface->report_error( "Reading subset of files not supported for Sms." );
    return MB_UNSUPPORTED_OPERATION;
  }

  mCurrentMeshHandle = 0;
  const MBErrorCode result = load_file_impl( filename );
  
    // If file read has failed, destroy anything that was
    // created during the read.
  if (MB_SUCCESS != result && mCurrentMeshHandle)
  {
    MBRange entities;
    mdbImpl->get_entities_by_handle( mCurrentMeshHandle, entities );
    entities.insert( mCurrentMeshHandle );
    mdbImpl->delete_entities( entities );
    mCurrentMeshHandle = 0;
  }
  
  file_set = mCurrentMeshHandle;
  return result;
}

MBErrorCode ReadSms::load_file_impl( const char* filename )
{
  bool warned = false;
  
  MBErrorCode result = mdbImpl->tag_get_handle( GLOBAL_ID_TAG_NAME, globalId );
  if (MB_TAG_NOT_FOUND == result)
    result = mdbImpl->tag_create( GLOBAL_ID_TAG_NAME,
                                  sizeof(int), MB_TAG_SPARSE,
                                  MB_TYPE_INTEGER, globalId, 0 );
  CHECK("Failed to create gid tag.");
    
  result = mdbImpl->tag_get_handle("PARAMETER_COORDS", paramCoords );
  double dum_params[] = {0.0, 0.0, 0.0};
  if (MB_TAG_NOT_FOUND == result) {
    result = mdbImpl->tag_create("PARAMETER_COORDS",
                                 3*sizeof(double), MB_TAG_DENSE,
                                 MB_TYPE_DOUBLE, paramCoords, 0 );
  }
  CHECK("Failed to create param coords tag.");
  if (MB_SUCCESS != result)
    return result;
    
  result = mdbImpl->tag_get_handle(GEOM_DIMENSION_TAG_NAME, geomDimension);
  if (MB_TAG_NOT_FOUND == result) {
    int dum_dim = -1;
    result = mdbImpl->tag_create(GEOM_DIMENSION_TAG_NAME,
                                 sizeof(int), MB_TAG_SPARSE,
                                 MB_TYPE_INTEGER, geomDimension, &dum_dim);
  }
  CHECK("Failed to create geom dim tag.");
  if (MB_SUCCESS != result)
    return result;
    
    // Open file
  FILE* file_ptr = fopen( filename, "r" );
  if (!file_ptr)
  {
    readMeshIface->report_error( "%s: %s\n", filename, strerror(errno) );
    return MB_FILE_DOES_NOT_EXIST;
  }

  char line[256];
  int file_type;
  fscanf(file_ptr, "%s %d", line, &file_type);

  if (3 == file_type) {
    result = read_parallel_info(file_ptr);
    CHECK("Failed to read parallel info.");
  }

  int nregions, nfaces, nedges, nvertices, npoints;
  fscanf(file_ptr, "%d %d %d %d %d", &nregions, &nfaces, &nedges,
         &nvertices, &npoints);

    // create the vertices
  std::vector<double*> coord_arrays;
  MBEntityHandle vstart = 0;
  result = readMeshIface->get_node_arrays( 3, nvertices, MB_START_ID, 
                                           vstart, coord_arrays );
  CHECK("Failed to get node arrays.");
  if (MB_SUCCESS != result)
    return result;
  
  MBEntityHandle this_gent, new_handle;
  std::vector<MBEntityHandle> gentities[4];
  int gent_id, dum_int;
  int gent_type, num_connections;
  
  for(int i = 0; i < nvertices; i++)
  {
    fscanf(file_ptr, "%d", &gent_id); 
    if (!gent_id) continue;

    fscanf(file_ptr,"%d %d %lf %lf %lf", &gent_type, &num_connections,
           coord_arrays[0]+i, coord_arrays[1]+i, coord_arrays[2]+i);
    
    result = get_set(gentities, gent_type, gent_id, geomDimension, this_gent);

    new_handle = vstart + i;
    result = mdbImpl->add_entities(this_gent, &new_handle, 1);
    CHECK("Adding vertex to geom set failed.");
    if (MB_SUCCESS != result) return result;

    switch(gent_type)
    {
      case 1:
          fscanf(file_ptr, "%le", dum_params);
          result = mdbImpl->tag_set_data(paramCoords, &new_handle, 1, dum_params);
          CHECK("Failed to set param coords tag for vertex.");
          if (MB_SUCCESS != result) return result;
          break;
      case 2:
          fscanf(file_ptr, "%le %le %d", dum_params, dum_params+1, &dum_int);
          dum_params[2] = dum_int;
          result = mdbImpl->tag_set_data(paramCoords, &new_handle, 1, dum_params);
          CHECK("Failed to set param coords tag for vertex.");
          if (MB_SUCCESS != result) return result;
          break;
          
      default: break;
    }
  } // end of reading vertices

// *******************************
//	Read Edges
// *******************************

  int vert1, vert2, num_pts;
  std::vector<MBEntityHandle> everts(2);
  MBEntityHandle estart, *connect;
  result = readMeshIface->get_element_array(nedges, 2, MBEDGE, 1, estart, connect);
  CHECK("Failed to create array of edges.");
  if (MB_SUCCESS != result) return result;

  for(int i = 0; i < nedges; i++)
  {
    fscanf(file_ptr,"%d",&gent_id);
    if (!gent_id) continue;

    fscanf(file_ptr, "%d %d %d %d %d", &gent_type, &vert1, &vert2, 
           &num_connections, &num_pts);
    connect[0] = vstart + vert1 - 1;
    connect[1] = vstart + vert2 - 1;
    if (num_pts > 1 && !warned) {
      std::cout << "Warning: num_points > 1 not supported; choosing last one." << std::endl;
      warned = true;
    }

    result = get_set(gentities, gent_type, gent_id, geomDimension, this_gent);
    CHECK("Problem getting geom set for edge.");

    new_handle = estart + i;
    result = mdbImpl->add_entities(this_gent, &new_handle, 1);
    CHECK("Failed to add edge to geom set.");
    if (MB_SUCCESS != result) return result;

    connect += 2;

    for(int j = 0; j < num_pts; j++) {
      switch(gent_type) {
        case 1: 
            fscanf(file_ptr, "%le", dum_params);
            result = mdbImpl->tag_set_data(paramCoords, &new_handle, 1, dum_params);
            CHECK("Failed to set param coords tag for edge.");
            if (MB_SUCCESS != result) return result;
            break;
        case 2: 
            fscanf(file_ptr, "%le %le %d", dum_params, dum_params+1, &dum_int);
            dum_params[2] = dum_int;
            result = mdbImpl->tag_set_data(paramCoords, &new_handle, 1, dum_params);
            CHECK("Failed to set param coords tag for edge.");
            if (MB_SUCCESS != result) return result;
            break;
        default: 
            break;
      }
    }

  }  // end of reading edges

// *******************************
//	Read Faces
// *******************************
  std::vector<MBEntityHandle> bound_ents, bound_verts, new_faces;
  int bound_id;
  MBRange shverts;
  new_faces.resize(nfaces);
  int num_bounding;
    
  for(int i = 0; i < nfaces; i++)
  {
    fscanf(file_ptr, "%d", &gent_id);
    if(!gent_id) continue;

    fscanf(file_ptr,"%d %d", &gent_type, &num_bounding);

    result = get_set(gentities, gent_type, gent_id, geomDimension, this_gent);
    CHECK("Problem getting geom set for face.");

    bound_ents.resize(num_bounding+1);
    bound_verts.resize(num_bounding);
    for(int j = 0; j < num_bounding; j++) {
      fscanf(file_ptr, "%d ", &bound_id);
      if (0 > bound_id) bound_id = abs(bound_id);
      assert(0 < bound_id && bound_id <= nedges);
      bound_ents[j] = estart + abs(bound_id) - 1;
    }

      // convert edge-based model to vertex-based one
    for (int j = 0; j < num_bounding; j++) {
      if (j == num_bounding-1) bound_ents[j+1] = bound_ents[0];
      result = mdbImpl->get_adjacencies(&bound_ents[j], 2, 0, false, shverts);
      CHECK("Failed to get vertices bounding edge.");
      if (MB_SUCCESS != result) return result;
      assert(shverts.size() == 1);
      bound_verts[j] = *shverts.begin();
      shverts.clear();
    }

    result = mdbImpl->create_element((MBEntityType)(MBTRI+num_bounding-3),
                                     &bound_verts[0], bound_verts.size(), 
                                     new_faces[i]);
    CHECK("Failed to create edge.");
    if (MB_SUCCESS != result) return result;

    result = mdbImpl->add_entities(this_gent, &new_faces[i], 1);
    CHECK("Failed to add edge to geom set.");
    if (MB_SUCCESS != result) return result;

    fscanf(file_ptr, "%d", &num_pts);
    if(!num_pts) continue;

    for(int j = 0; j < num_pts; j++) {
      switch(gent_type) {
        case 1: 
            fscanf(file_ptr, "%le", dum_params);
            result = mdbImpl->tag_set_data(paramCoords, &new_faces[i], 1, dum_params);
            CHECK("Failed to set param coords tag for face.");
            if (MB_SUCCESS != result) return result;
            break;
        case 2: 
            fscanf(file_ptr, "%le %le %d", dum_params, dum_params+1, &dum_int);
            dum_params[2] = dum_int;
            result = mdbImpl->tag_set_data(paramCoords, &new_faces[i], 1, dum_params);
            CHECK("Failed to set param coords tag for face.");
            if (MB_SUCCESS != result) return result;
            break;
        default: 
            break;
      }
    }

  } // end of reading faces


// *******************************
//	Read Regions
// *******************************
  int sense[MB_MAX_SUB_ENTITIES];
  bound_verts.resize(MB_MAX_SUB_ENTITIES);

  for(int i = 0; i < nregions; i++)
  {
    fscanf(file_ptr, "%d", &gent_id);
    if (!gent_id) continue;
    result = get_set(gentities, 3, gent_id, geomDimension, this_gent);
    CHECK("Couldn't get geom set for region.");
    fscanf(file_ptr, "%d", &num_bounding);
    bound_ents.resize(num_bounding);
    for(int j = 0; j < num_bounding; j++) {
      fscanf(file_ptr, "%d ", &bound_id);
      assert(abs(bound_id) < (int)new_faces.size()+1 && bound_id);
      sense[j] = (bound_id < 0) ? -1 : 1;
      bound_ents[j] = new_faces[abs(bound_id)-1];
    }

    MBEntityType etype;
    result = readMeshIface->get_ordered_vertices(&bound_ents[0], sense, 
                                                 num_bounding,
                                                 3, &bound_verts[0], etype);
    CHECK("Failed in get_ordered_vertices.");
        
      // make the element
    result = mdbImpl->create_element(etype, &bound_verts[0], 
                                     MBCN::VerticesPerEntity(etype), new_handle);
    CHECK("Failed to create region.");
    if (MB_SUCCESS != result) return result;

    result = mdbImpl->add_entities(this_gent, &new_handle, 1);
    CHECK("Failed to add region to geom set.");
    if (MB_SUCCESS != result) return result;

    fscanf(file_ptr, "%d ", &dum_int);

  } // end of reading regions

  return MB_SUCCESS;
}

MBErrorCode ReadSms::get_set(std::vector<MBEntityHandle> *sets,
                             int set_dim, int set_id,
                             MBTag dim_tag,
                             MBEntityHandle &this_set) 
{
  MBErrorCode result = MB_SUCCESS;
  
  if ((int)sets[set_dim].size() <= set_id || 
      !sets[set_dim][set_id]) {
    if ((int)sets[set_dim].size() <= set_id) 
      sets[set_dim].resize(set_id+1, 0);
      
    if (!sets[set_dim][set_id]) {
      result = mdbImpl->create_meshset(MESHSET_SET, 
                                       sets[set_dim][set_id]);
      if (MB_SUCCESS != result) return result;
      result = mdbImpl->tag_set_data(globalId, 
                                     &sets[set_dim][set_id], 1,
                                     &set_id);
      if (MB_SUCCESS != result) return result;
      result = mdbImpl->tag_set_data(dim_tag,
                                     &sets[set_dim][set_id], 1,
                                     &set_dim);
      if (MB_SUCCESS != result) return result;
    }
  }

  this_set = sets[set_dim][set_id];

  return result;
}

MBErrorCode ReadSms::read_parallel_info(FILE *file_ptr) 
{
  return MB_FAILURE;
    /*
    // read partition info
  int nparts, part_id, num_ifaces, num_corner_ents;
  fscanf(file_ptr, "%d %d %d %d", &nparts, &part_id, &num_ifaces, &num_corner_ents);
  
    // read interfaces
  int iface_id, iface_dim, iface_own, num_iface_corners;
  MBEntityHandle this_iface;
  std::vector<int> *iface_corners;
  for (int i = 0; i < num_ifaces; i++) {
    fscanf(file_ptr, "%d %d %d %d", &iface_id, &iface_dim, &iface_own,
           &num_iface_corners);
    
    result = get_set(iface_id, iface_dim, iface_own, this_iface);
    CHECK("Failed to make iface set.");
    
      // read the corner ids and store them on the set for now
    iface_corners = new std::vector<int>(num_iface_corners);
    for (int j = 0; j < num_iface_corners; j++)
      fscanf(file_ptr, "%d", &(*iface_corners)[j]);

    result = tag_set_data(ifaceCornerTag, &this_iface, 1,
                          &iface_corners);
    CHECK("Failed to set iface corner tag.");

  }

    // interface data has been read
  return MB_SUCCESS;
    */
}

