
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
#include <fstream>

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
    : mbImpl(impl), oFile(NULL), mCurrentMeshHandle(0)
{
  assert(impl != NULL);

  std::string iface_name = "MBWriteUtilIface";
  impl->query_interface(iface_name, reinterpret_cast<void**>(&mWriteIface));

  // initialize in case tag_get_handle fails below
  //! get and cache predefined tag handles
  int dum_val = 0;
  MBErrorCode result = impl->tag_get_handle(MATERIAL_SET_TAG_NAME,  mMaterialSetTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(MATERIAL_SET_TAG_NAME, sizeof(int), MB_TAG_SPARSE, mMaterialSetTag,
                              &dum_val);
  
  result = impl->tag_get_handle(DIRICHLET_SET_TAG_NAME, mDirichletSetTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(DIRICHLET_SET_TAG_NAME, sizeof(int), MB_TAG_SPARSE, mDirichletSetTag,
                              &dum_val);
  
  result = impl->tag_get_handle(NEUMANN_SET_TAG_NAME,   mNeumannSetTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(NEUMANN_SET_TAG_NAME, sizeof(int), MB_TAG_SPARSE, mNeumannSetTag,
                              &dum_val);
  
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
                                 const bool /*overwrite*/,
                                 const MBEntityHandle *ent_handles,
                                 const int num_sets,
                                 std::vector<std::string>&, int )
{
  assert(0 != mMaterialSetTag &&
         0 != mNeumannSetTag &&
         0 != mDirichletSetTag);

  std::vector<MBEntityHandle> matsets, dirsets, neusets, entities;

    // separate into material sets, dirichlet sets, neumann sets

  if (num_sets == 0) {
      // default to all defined sets
    MBRange this_range;
    mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &mMaterialSetTag, NULL, 1, this_range);
    std::copy(this_range.begin(), this_range.end(), std::back_inserter(matsets));
    this_range.clear();
    mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &mDirichletSetTag, NULL, 1, this_range);
    std::copy(this_range.begin(), this_range.end(), std::back_inserter(dirsets));
    this_range.clear();
    mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &mNeumannSetTag, NULL, 1, this_range);
    std::copy(this_range.begin(), this_range.end(), std::back_inserter(neusets));
  }
  else {
    int dummy;
    for (const MBEntityHandle *iter = ent_handles; iter < ent_handles+num_sets; iter++) 
    {
      if (MB_SUCCESS == mbImpl->tag_get_data(mMaterialSetTag, &(*iter), 1, &dummy))
        matsets.push_back(*iter);
      else if (MB_SUCCESS == mbImpl->tag_get_data(mDirichletSetTag, &(*iter), 1, &dummy))
        dirsets.push_back(*iter);
      else if (MB_SUCCESS == mbImpl->tag_get_data(mNeumannSetTag, &(*iter), 1, &dummy))
        neusets.push_back(*iter);
    }
  }
  
    // if there is nothing to write just return.
  if (matsets.empty()) {
      mWriteIface->report_error("Nobody to write.  I come back tomoddow.");
      return MB_FILE_WRITE_ERROR;
  }

  if (!dirsets.empty() || !neusets.empty())
    mWriteIface->report_error("Warning, can't write dirichlet nor neumann sets to this format.");

  std::vector<WriteVtk::MaterialSetData> matset_info;
  std::vector<WriteVtk::DirichletSetData> dirset_info;
  std::vector<WriteVtk::NeumannSetData> neuset_info;

  MeshInfo mesh_info;
  
  matset_info.clear();
  if(gather_mesh_information(mesh_info, matset_info, neuset_info, dirset_info,
                             matsets, neusets, dirsets) != MB_SUCCESS)
  {
    reset_matset(matset_info);
    return MB_FAILURE;
  }

  // try to open the file after gather mesh info succeeds
  MBErrorCode result = open_file(file_name);
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

  if( write_nodes(mesh_info.num_nodes, mesh_info.nodes, mesh_info.num_dim) != MB_SUCCESS )
  {
    reset_matset(matset_info);
    return MB_FAILURE;
  }

  if( write_matsets(mesh_info, matset_info, neuset_info) )
  {
    reset_matset(matset_info);
    return MB_FAILURE;
  }

  oFile.flush();
  oFile.close();
  
  return MB_SUCCESS;
}

MBErrorCode WriteVtk::gather_mesh_information(MeshInfo &mesh_info,
                                                   std::vector<WriteVtk::MaterialSetData> &matset_info,
                                                   std::vector<WriteVtk::NeumannSetData> &,
                                                   std::vector<WriteVtk::DirichletSetData> &,
                                                   std::vector<MBEntityHandle> &matsets,
                                                   std::vector<MBEntityHandle> &,
                                                   std::vector<MBEntityHandle> &)
{

  std::vector<MBEntityHandle>::iterator vector_iter, end_vector_iter;

  mesh_info.num_nodes = 0;
  mesh_info.num_elements = 0;
  mesh_info.num_matsets = 0;
  
  int id = 0;

  vector_iter= matsets.begin();
  end_vector_iter = matsets.end();

  mesh_info.num_matsets = matsets.size();

  std::vector<MBEntityHandle> parent_meshsets;

  // clean out the bits for the element mark
  mbImpl->tag_delete(mEntityMark);
  mbImpl->tag_create("WriteVtk element mark", 1, MB_TAG_BIT, mEntityMark, NULL);

  int highest_dimension_of_element_matsets = 0;

  for(vector_iter = matsets.begin(); vector_iter != matsets.end(); vector_iter++)
  {
       
    WriteVtk::MaterialSetData matset_data;
    matset_data.elements = new MBRange;

    //for the purpose of qa records, get the parents of these matsets 
    if( mbImpl->get_parent_meshsets( *vector_iter, parent_meshsets ) != MB_SUCCESS )
      return MB_FAILURE;

    // get all Entity Handles in the mesh set
    MBRange dummy_range;
    mbImpl->get_entities_by_handle(*vector_iter, dummy_range, true );

      // find the dimension of the last entity in this range
    MBRange::iterator entity_iter = dummy_range.end();
    entity_iter = dummy_range.end();
    entity_iter--;
    int this_dim = MBCN::Dimension(TYPE_FROM_HANDLE(*entity_iter));
    entity_iter = dummy_range.begin();
    while (entity_iter != dummy_range.end() &&
           MBCN::Dimension(TYPE_FROM_HANDLE(*entity_iter)) != this_dim)
      entity_iter++;
    
    if (entity_iter != dummy_range.end())
      std::copy(entity_iter, dummy_range.end(), mb_range_inserter(*(matset_data.elements)));

    assert(matset_data.elements->begin() == matset_data.elements->end() ||
           MBCN::Dimension(TYPE_FROM_HANDLE(*(matset_data.elements->begin()))) == this_dim);
    
    // get the matset's id
    if(mbImpl->tag_get_data(mMaterialSetTag, &(*vector_iter), 1, &id) != MB_SUCCESS ) {
      mWriteIface->report_error("Couldn't get matset id from a tag for an element matset.");
      return MB_FAILURE;
    }
    
    matset_data.id = id; 
    matset_data.number_attributes = 0;
 
     // iterate through all the elements in the meshset
    MBRange::iterator elem_range_iter, end_elem_range_iter;
    elem_range_iter = matset_data.elements->begin();
    end_elem_range_iter = matset_data.elements->end();

      // get the entity type for this matset, verifying that it's the same for all elements
      // THIS ASSUMES HANDLES SORT BY TYPE!!!
    MBEntityType entity_type = TYPE_FROM_HANDLE(*elem_range_iter);
    end_elem_range_iter--;
    if (entity_type != TYPE_FROM_HANDLE(*(end_elem_range_iter++))) {
      mWriteIface->report_error("Entities in matset %i not of common type", id);
      return MB_FAILURE;
    }

    int dimension = MBCN::Dimension(entity_type);

    if( dimension > highest_dimension_of_element_matsets )
      highest_dimension_of_element_matsets = dimension;

    matset_data.moab_type = mbImpl->type_from_handle(*(matset_data.elements->begin()));
    if (MBMAXTYPE == matset_data.moab_type) return MB_FAILURE;
    
    std::vector<MBEntityHandle> tmp_conn;
    mbImpl->get_connectivity(&(*(matset_data.elements->begin())), 1, tmp_conn);
    matset_data.element_type = 
      ExoIIUtil::get_element_type_from_num_verts(tmp_conn.size(), entity_type, dimension);
    
    if (matset_data.element_type == EXOII_MAX_ELEM_TYPE) {
      mWriteIface->report_error("Element type in matset %i didn't get set correctly", id);
      return MB_FAILURE;
    }
    
    matset_data.number_nodes_per_element = ExoIIUtil::VerticesPerElement[matset_data.element_type];

    // number of nodes for this matset
    matset_data.number_elements = matset_data.elements->size();

    // total number of elements
    mesh_info.num_elements += matset_data.number_elements;

    // get the nodes for the elements
    mWriteIface->gather_nodes_from_elements(*matset_data.elements, mEntityMark, mesh_info.nodes);

    matset_info.push_back( matset_data );
  
  }
 

  //if user hasn't entered dimension, we figure it out
  if( mesh_info.num_dim == 0 )
  {
    //never want 1 or zero dimensions
    if( highest_dimension_of_element_matsets < 2 )
      mesh_info.num_dim = 3;
    else
      mesh_info.num_dim = highest_dimension_of_element_matsets;
  }

  if (mesh_info.num_dim < 3) {
      // check against mesh interface value
    int tmp_dim;
    MBErrorCode result = mbImpl->get_dimension(tmp_dim);
    if (MB_SUCCESS == result && mesh_info.num_dim < (unsigned int) tmp_dim)
      mesh_info.num_dim = tmp_dim;
  }

  MBRange::iterator range_iter, end_range_iter;
  range_iter = mesh_info.nodes.begin();
  end_range_iter = mesh_info.nodes.end();

  mesh_info.num_nodes = mesh_info.nodes.size(); 

    // (no neusets or dirsets in this format)
  return MB_SUCCESS;
}

MBErrorCode WriteVtk::write_nodes(const int num_nodes, const MBRange& nodes, const int dimension)
{
  //see if should transform coordinates
  MBErrorCode result;
  MBTag trans_tag;
  result = mbImpl->tag_get_handle( MESH_TRANSFORM_TAG_NAME, trans_tag);
  bool transform_needed = true;
  if( result == MB_TAG_NOT_FOUND )
    transform_needed = false;

  int num_coords_to_fill = transform_needed ? 3 : dimension;

  std::vector<double*> coord_arrays(3);
  coord_arrays[0] = new double[num_nodes];
  coord_arrays[1] = new double[num_nodes];
  coord_arrays[2] = NULL;

  if( num_coords_to_fill == 3 ) 
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
  std::cout << "POINTS " << num_nodes << " float" << std::endl;

  if (NULL != coord_arrays[2]) {
    for( int i=0; i<num_nodes; i++)
    {
      std::cout << coord_arrays[0][i] << " " 
            << coord_arrays[1][i] << " "
            << coord_arrays[2][i] << std::endl;
    }
  }
  else {
    for( int i=0; i<num_nodes; i++)
    {
      std::cout << coord_arrays[0][i] << " " 
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
                                    std::vector<WriteVtk::MaterialSetData> &matset_data,
                                    std::vector<WriteVtk::NeumannSetData> &)
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

  std::cout << "CELLS " << total_elems << " " << total_ints << std::endl;

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
        std::cout << connect[j] << " ";
      
      std::cout << std::endl;
    }
  }

  std::cout << std::endl << "CELL_TYPES " << total_elems << std::endl;
  for (i = 0; i < matset_data.size(); i++) {
    int elem_type = VtkUtil::vtkElemType[matset_data[i].moab_type];
    if (-1 == elem_type) continue;
    for (int j = 0; j < matset_data[i].number_elements; j++)
      std::cout << elem_type << std::endl;
  }
  
  return MB_SUCCESS;
}

MBErrorCode WriteVtk::initialize_file(MeshInfo &)
{
    // perform the initializations
  std::cout << "# vtk DataFile Version 2.0" << std::endl;
  std::cout << "MOAB vtk data file" << std::endl;
  std::cout << "ASCII" << std::endl;
  std::cout << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // wait until actual write operations to count how many and actually write 

  return MB_SUCCESS;
}


MBErrorCode WriteVtk::open_file(const char* filename)
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

  oFile.open(filename, std::ofstream::out | std::ofstream::trunc);

    // file couldn't be opened
  if(!oFile.is_open())
  {
    mWriteIface->report_error("Cannot open %s", filename);
    return MB_FAILURE;
  }
   
  return MB_SUCCESS;
}


  
