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


#ifdef WIN32
#ifdef _DEBUG
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif
#endif


#include "WriteVtk.hpp"
#include "VtkUtil.hpp"
#include "ExoIIUtil.hpp"

#include <utility>
#include <algorithm>
#include <time.h>
#include <string>
#include <vector>
#include <stdio.h>
#include <iostream>

#include "MBInterface.hpp"
#include "MBRange.hpp"
#include "MBCN.hpp"
#include "assert.h"
#include "MBInternals.hpp"
#include "ExoIIUtil.hpp"
#include "MBTagConventions.hpp"

#define INS_ID(stringvar, prefix, id) \
sprintf(stringvar, prefix, id)

MBWriterIface *WriteVtk::factory( MBInterface* iface )
  { return new WriteVtk( iface ); }

WriteVtk::WriteVtk(MBInterface *impl) 
    : mbImpl(impl), mCurrentMeshHandle(0)
{
  assert(impl != NULL);

  std::string iface_name = "MBWriteUtilIface";
  impl->query_interface(iface_name, reinterpret_cast<void**>(&mWriteIface));

  // initialize in case tag_get_handle fails below
  //! get and cache predefined tag handles
  int dum_val = 0;
  MBErrorCode result = impl->tag_get_handle(MATERIAL_SET_TAG_NAME,  mMaterialSetTag);
  if (MB_TAG_NOT_FOUND == result)
    mMaterialSetTag = 0;
  
  result = impl->tag_get_handle(HAS_MID_NODES_TAG_NAME, mHasMidNodesTag);
  if (MB_TAG_NOT_FOUND == result) {
    int dum_val_array[] = {0, 0, 0, 0};
    result = impl->tag_create(HAS_MID_NODES_TAG_NAME, 4*sizeof(int), MB_TAG_SPARSE, mHasMidNodesTag,
                              dum_val_array);
  }
  
  result = impl->tag_get_handle(GLOBAL_ID_TAG_NAME, mGlobalIdTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_SPARSE, mGlobalIdTag,
                              &dum_val);
  
  dum_val = -1;
  result = impl->tag_get_handle("__matSetIdTag", mMatSetIdTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create("__matSetIdTag", sizeof(int), MB_TAG_DENSE, mMatSetIdTag,
                              &dum_val);
  

  impl->tag_create("WriteVtk element mark", 1, MB_TAG_BIT, mEntityMark, NULL);

}

WriteVtk::~WriteVtk() 
{
  if (oFile.is_open())
    oFile.close();
  
  std::string iface_name = "MBWriteUtilIface";
  mbImpl->release_interface(iface_name, mWriteIface);

  mbImpl->tag_delete(mEntityMark);
}

void WriteVtk::reset_matset(std::vector<WriteVtk::MaterialSetData> &matset_info)
{
  std::vector<WriteVtk::MaterialSetData>::iterator iter;
  
  for (iter = matset_info.begin(); iter != matset_info.end(); iter++)
  {
    delete (*iter).elements;
  }
}

MBErrorCode WriteVtk::write_file(const char *file_name, 
                                 const bool overwrite,
                                 const MBEntityHandle *ent_handles,
                                 const int num_sets,
                                 std::vector<std::string>&, int )
{
  std::vector<MBEntityHandle> matsets, entities;

    // separate into material sets, dirichlet sets, neumann sets

  if (num_sets == 0) {
      // default to all material sets
    MBRange this_range;
    if (mMaterialSetTag) {
      mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &mMaterialSetTag, NULL, 1, this_range);
      std::copy(this_range.begin(), this_range.end(), std::back_inserter(matsets));
      this_range.clear();
    }
      // If no material sets, will do entire mesh
  }
  else {
    std::copy(ent_handles, ent_handles+num_sets, std::back_inserter(matsets));
  }

  std::vector<WriteVtk::MaterialSetData> matset_info;

  MeshInfo mesh_info;
  
  matset_info.clear();
  if(gather_mesh_information(mesh_info, matset_info, matsets) != MB_SUCCESS)
  {
    reset_matset(matset_info);
    return MB_FAILURE;
  }

  // try to open the file after gather mesh info succeeds
  MBErrorCode result = open_file(file_name, overwrite);
  if (MB_FAILURE == result) {
    reset_matset(matset_info);
    return result;
  }

  result = initialize_file(mesh_info);
  if(result != MB_SUCCESS)
  {
    reset_matset(matset_info);
    return MB_FAILURE;
  }

  if( write_nodes(mesh_info.nodes) != MB_SUCCESS )
  {
    reset_matset(matset_info);
    return MB_FAILURE;
  }

  if( write_matsets(mesh_info, matset_info) )
  {
    reset_matset(matset_info);
    return MB_FAILURE;
  }

  oFile.sync();
  oFile.close();
  
  return MB_SUCCESS;
}

MBErrorCode WriteVtk::gather_mesh_information(MeshInfo &mesh_info,
                                              std::vector<WriteVtk::MaterialSetData> &matset_info,
                                              std::vector<MBEntityHandle> &matsets )
{
  MBErrorCode rval;
  
  // clean out the bits for the element mark
  mbImpl->tag_delete(mEntityMark);
  mbImpl->tag_create("WriteVtk element mark", 1, MB_TAG_BIT, mEntityMark, NULL);

  if (matsets.empty())
  {
      // if no input sets, do entire mesh (root set)
    rval = gather_mesh_information(matset_info, 0);
    if (MB_SUCCESS != rval)
      return MB_FAILURE;
  }
  else
  {
      // get information for each set
    for (std::vector<MBEntityHandle>::const_iterator iter = matsets.begin();
         iter != matsets.end(); ++iter)
    {
      rval = gather_mesh_information(matset_info, *iter);
      if (MB_SUCCESS != rval)
        return MB_FAILURE;
    }
  }
  
    // Now get nodes from elements, count total elements, etc.
  mesh_info.num_elements = 0;
  mesh_info.num_matsets = matset_info.size();
  for (std::vector<WriteVtk::MaterialSetData>::iterator iter = matset_info.begin();
       iter != matset_info.end(); ++iter)
  {
    mesh_info.num_elements += iter->number_elements;
    rval = mWriteIface->gather_nodes_from_elements( *(iter->elements), mEntityMark, mesh_info.nodes);
    if (MB_SUCCESS != rval)
      return rval;
  }
  mesh_info.num_nodes = mesh_info.nodes.size();

  return MB_SUCCESS;
}

MBErrorCode WriteVtk::gather_mesh_information(std::vector<WriteVtk::MaterialSetData> &matset_info,
                                              MBEntityHandle matset )
{
    // Get elements of highest dimension in entity set
  MBRange elements;
  for (int dim = 3; dim > 0 && elements.empty(); --dim)
    mbImpl->get_entities_by_dimension( matset, dim, elements, true );
  
    // Now subdivide the range by element type
  MBRange::iterator iter, last;
  for (iter = elements.begin(); iter != elements.end(); iter = last)
  {
      // Find beginning of next element type in range
    int err;
    MBEntityType type = TYPE_FROM_HANDLE( *iter );
    MBEntityHandle handle = CREATE_HANDLE( type+1, 0, err );
    last = elements.lower_bound( iter, elements.end(), handle );
    
      // Create MaterialSetData containing the elements
    WriteVtk::MaterialSetData matset_data;
    matset_data.elements = new MBRange;
    matset_data.elements->merge( iter, last );
    matset_data.number_elements = matset_data.elements->size();
    matset_data.moab_type = type;
    matset_data.number_nodes_per_element = MBCN::VerticesPerEntity( type );
    
      // Add to list
    matset_info.push_back( matset_data );
  }
  
  return MB_SUCCESS;
}


MBErrorCode WriteVtk::write_nodes(const MBRange& nodes)
{
  const int num_nodes = nodes.size();
  
  //see if should transform coordinates
  MBErrorCode result;
  MBTag trans_tag;
  result = mbImpl->tag_get_handle( MESH_TRANSFORM_TAG_NAME, trans_tag);
  bool transform_needed = true;
  if( result == MB_TAG_NOT_FOUND )
    transform_needed = false;

  int dimension;
  mbImpl->get_dimension(dimension);

  std::vector<double*> coord_arrays(3);
  coord_arrays[0] = new double[num_nodes];
  coord_arrays[1] = new double[num_nodes];
  coord_arrays[2] = NULL;

  if( transform_needed || dimension == 3 ) 
    coord_arrays[2] = new double[num_nodes];
 
  result = mWriteIface->get_node_arrays(dimension, num_nodes, nodes, 
                                        mGlobalIdTag, 0, coord_arrays);
  if(result != MB_SUCCESS)
  {
    delete [] coord_arrays[0];
    delete [] coord_arrays[1];
    if(coord_arrays[2]) delete [] coord_arrays[2];
    return result;
  }

  if( transform_needed )
  {
    double trans_matrix[16]; 
    result = mbImpl->tag_get_data( trans_tag, NULL, 0, trans_matrix ); 
    if (MB_SUCCESS != result) {
      mWriteIface->report_error("Couldn't get transform data.");
      return result;
    }
      
    for( int i=0; i<num_nodes; i++)
    {

      double vec1[3];
      double vec2[3];

      vec2[0] =  coord_arrays[0][i];
      vec2[1] =  coord_arrays[1][i];
      vec2[2] =  coord_arrays[2][i];

      for( int row=0; row<3; row++ )
      {
        vec1[row] = 0.0;
        for( int col = 0; col<3; col++ )
        {
          vec1[row] += ( trans_matrix[ (row*4)+col ] * vec2[col] );
        }
      }

      coord_arrays[0][i] = vec1[0];
      coord_arrays[1][i] = vec1[1];
      coord_arrays[2][i] = vec1[2];

    }
  }


  // write the nodes 
  oFile << "POINTS " << num_nodes << " float" << std::endl;

  if (NULL != coord_arrays[2]) {
    for( int i=0; i<num_nodes; i++)
    {
      oFile << coord_arrays[0][i] << " " 
            << coord_arrays[1][i] << " "
            << coord_arrays[2][i] << std::endl;
    }
  }
  else {
    for( int i=0; i<num_nodes; i++)
    {
      oFile << coord_arrays[0][i] << " " 
            << coord_arrays[1][i] << " "
            << "0.0" << std::endl;
    }
  }

  // clean up
  delete [] coord_arrays[0];
  delete [] coord_arrays[1];
  if(coord_arrays[2]) 
    delete [] coord_arrays[2];

  return MB_SUCCESS;

}

MBErrorCode WriteVtk::write_matsets(MeshInfo &,
                   std::vector<WriteVtk::MaterialSetData> &matset_data )
{
  unsigned int i;
  std::vector<int> connect;
  const MBEntityHandle *connecth;
  int num_connecth;
  MBErrorCode result;
  int total_elems = 0, total_ints = 0;

    // don't usually have anywhere near 31 nodes per element
  connect.reserve(31);
  MBRange::iterator rit;

    // three loops, first to count, second to write connectivity, 
    // third to write cell types
  WriteVtk::MaterialSetData matset;
  for (i = 0; i < matset_data.size(); i++) {
      // make sure this element type is supported
    if (-1 == VtkUtil::vtkElemType[matset_data[i].moab_type]) {
      mWriteIface->report_error("Warning: won't be able to write out MOAB type %s, no"
                                " corresponding vtk type.", 
                                MBCN::EntityTypeName(matset_data[i].moab_type));
    }
    
    else {
      int tot_e = matset_data[i].elements->size();
      total_elems += tot_e;
      total_ints += (matset_data[i].number_nodes_per_element + 1) * tot_e;
    }
  }

  oFile << "CELLS " << total_elems << " " << total_ints << std::endl;

  for (i = 0; i < matset_data.size(); i++) {
    matset = matset_data[i];
    if (-1 == VtkUtil::vtkElemType[matset.moab_type]) continue;
    
      // 1st element is always the # nodes
    connect[0] = matset.number_nodes_per_element;
      
    for (rit = matset.elements->begin(); rit != matset.elements->end(); rit++) {
      
        // get the connectivity of this element
      result = mbImpl->get_connectivity(*rit, connecth, num_connecth);
      if (MB_SUCCESS != result) return result;
      
        // get the vertex ids
      result = mbImpl->tag_get_data(mGlobalIdTag, connecth, num_connecth, &connect[1]);
      if (MB_SUCCESS != result) return result;
      
        // write the data; <= because we have an extra member, # nodes
      for (int j = 0; j <= matset.number_nodes_per_element; j++) 
        oFile << connect[j] << " ";
      
      oFile << std::endl;
    }
  }

  oFile << "CELL_TYPES " << total_elems << std::endl;
  for (i = 0; i < matset_data.size(); i++) {
    int elem_type = VtkUtil::vtkElemType[matset_data[i].moab_type];
    if (-1 == elem_type) continue;
    for (int j = 0; j < matset_data[i].number_elements; j++)
      oFile << elem_type << std::endl;
  }
  
  return MB_SUCCESS;
}

MBErrorCode WriteVtk::initialize_file(MeshInfo &)
{
    // perform the initializations
  oFile << "# vtk DataFile Version 2.0" << std::endl;
  oFile << "MOAB vtk data file" << std::endl;
  oFile << "ASCII" << std::endl;
  oFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // wait until actual write operations to count how many and actually write 

  return MB_SUCCESS;
}


MBErrorCode WriteVtk::open_file(const char* filename, const bool overwrite)
{
    // check the file name
  if (NULL == strstr(filename, ".vtk")) {
    mWriteIface->report_error("Vtk files should end with .vtk ");
    return MB_FAILURE;
  }

    // not a valid filname
  if(strlen(filename) == 0)
  {
    mWriteIface->report_error("Output filename not specified");
    return MB_FAILURE;
  }

  if (overwrite) 
    oFile.open(filename, std::ios::out | std::ios::trunc);
  else {
    oFile.open(filename, std::ios::in);
    if (!oFile.fail()) {
        // didn't fail, which means it's already there, which is an error
      mWriteIface->report_error("Output filename already exists and overwrite not allowed.");
      oFile.close();
      return MB_FAILURE;
    }
    oFile.clear();
    oFile.open(filename, std::ios::out);
  }

    // file couldn't be opened
  if(!oFile.is_open())
  {
    mWriteIface->report_error("Cannot open %s", filename);
    return MB_FAILURE;
  }
   
  return MB_SUCCESS;
}


  
