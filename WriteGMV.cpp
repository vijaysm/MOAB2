
#ifdef WIN32
#ifdef _DEBUG
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif
#endif

#include "WriteGMV.hpp"

#include "MBInterface.hpp"
#include "MBInternals.hpp"
#include "MBRange.hpp"
#include "MBCN.hpp"
#include <fstream>

const char *WriteGMV::gmvTypeNames[] = {
  "",
  "line",
  "tri",
  "quad",
  "",
  "tet",
  "pyramid",
  "prism",
  "",
  "hex",
  "",
  ""
};

WriteGMV::WriteGMV(MBInterface *impl) 
    : mbImpl(impl), mCurrentMeshHandle(0)
{
  assert(impl != NULL);

  std::string iface_name = "MBWriteUtilIface";
  impl->query_interface(iface_name, reinterpret_cast<void**>(&mWriteIface));

  // initialize in case tag_get_handle fails below
  mMaterialSetTag  = 0;
  mDirichletSetTag = 0;
  mNeumannSetTag   = 0;
  mHasMidNodesTag  = 0;
  mGeomDimensionTag= 0;
  mGlobalIdTag= 0;

  //! get and cache predefined tag handles
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
  
  result = impl->tag_get_handle(GLOBAL_ID_TAG_NAME,             mGlobalIdTag);
  if (MB_TAG_NOT_FOUND == result)
    result = impl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_SPARSE, mGlobalIdTag,
                              &dum_val);
  

  impl->tag_create("__WriteGMV element mark", 1, MB_TAG_BIT, mEntityMark, NULL);

}

WriteGMV::~WriteGMV() 
{
  std::string iface_name = "MBWriteUtilIface";
  mbImpl->release_interface(iface_name, mWriteIface);

  mbImpl->tag_delete(mEntityMark);

}

MBErrorCode WriteGMV::write_file(const char *file_name,
                                 const MBEntityHandle output_set,
                                 const int user_dimension) 
{
    // general function for writing a mesh
  
  MBRange dum_range, elements, all_verts;
  MBErrorCode result;
  MBEntityType otype;

    // get elements to be output; iterate over type because of special handling
    // for polyhedra
  for (otype = MBCN::TypeDimensionMap[user_dimension].first;
       otype <= MBCN::TypeDimensionMap[user_dimension].second; otype++) {
    dum_range.clear();
    if (otype != MBPOLYHEDRON)
      result = mbImpl->get_entities_by_type(output_set, otype, dum_range, true);
    else
        // if type requested is polyhedra, we'll really be outputting polygons, so
        // put polygons in range
      result = mbImpl->get_entities_by_type(output_set, MBPOLYGON, dum_range, true);
    if (MB_SUCCESS != result) return result;

    std::copy(dum_range.begin(), dum_range.end(), mb_range_inserter(elements));
  }
  
    // gather the vertices in these elements
  result = mbImpl->get_adjacencies(elements, 0, false, all_verts, MBInterface::UNION);
  if (MB_SUCCESS != result) return result;
  
  int num_verts = all_verts.size();
  
    // allocate coordinate arrays and put pointers to them in a list
  double *xcoord = new double[num_verts];
  double *ycoord = new double[num_verts];
  double *zcoord = new double[num_verts];
  std::vector<double*> coord_arrays;
  coord_arrays.push_back(xcoord);
  coord_arrays.push_back(ycoord);
  coord_arrays.push_back(zcoord);
  
    // fill them in, writing id tags at the same time
  result = mWriteIface->get_node_arrays(3, num_verts, all_verts, mGlobalIdTag, coord_arrays);
  if (MB_SUCCESS != result) return result;

    // initialize file
  std::ofstream ofile(file_name);
  ofile << "gmvinput ascii" << std::endl;

  int i, j;
  
    //========================================
    // WRITE COORDINATE DATA TO FILE HERE

  ofile << "nodev " << num_verts << std::endl;
  for (i = 0; i < num_verts; i++) 
    ofile << xcoord[i] << " " << ycoord[i] << " " << zcoord[i] << std::endl;
  

    //========================================

  delete [] xcoord;
  delete [] ycoord;
  delete [] zcoord;
  
    // iterate over types in selected dimension

  MBRange sub_range;
  std::vector<int> connect;
  std::vector<MBEntityHandle> connecth;
  for (MBEntityType otype = MBCN::TypeDimensionMap[user_dimension].first;
       otype <= MBCN::TypeDimensionMap[user_dimension].second; otype++) {

    if (otype == MBPOLYGON || otype == MBPOLYHEDRON) continue;
      
      // get the first element of this type in the range, and one past the last
    MBRange::iterator lower = std::lower_bound(elements.begin(), elements.end(), 
                                               CREATE_HANDLE(otype, MB_START_ID, i)),
      upper = std::lower_bound(elements.begin(), elements.end(), 
                               CREATE_HANDLE(otype+1, MB_START_ID, i));

    if (lower == upper) continue;
    
      // copy these elements into a subrange
    sub_range.clear();
    std::copy(lower, upper, mb_range_inserter(sub_range));

      // make sure the connectivity array is big enough
    int verts_per = MBCN::VerticesPerEntity(otype);
    if (connect.size() < verts_per*sub_range.size())
      connect.reserve(verts_per*sub_range.size());
    
      // get the connectivity
    result = mWriteIface->get_element_array(sub_range.size(),
                                            verts_per,
                                            mGlobalIdTag, sub_range,
                                            mGlobalIdTag, 1, &connect[0]);
    if (MB_SUCCESS != result) return result;

      //========================================
      // WRITE CONNECTIVITY DATA TO FILE HERE

    for (i = 0; i < (int) sub_range.size(); i++) {
      ofile << gmvTypeNames[otype] << " " << verts_per << std::endl;
      for (j = i*verts_per; j < (int) (i+1)*verts_per; j++)
        ofile << connect[j] << " ";
      ofile << std::endl;
    }

      //========================================
  }

    // write polygons/hedra, if any
  MBRange polygons, polyhedra;
  result = mbImpl->get_entities_by_type(output_set, MBPOLYGON, polygons, true);
  if (MB_SUCCESS != result) return result;
  
  result = mbImpl->get_entities_by_type(output_set, MBPOLYHEDRON, polyhedra, true);
  if (MB_SUCCESS != result) return result;

  if (polygons.size() == 0) return result;
  
    // mark polyhedra with global ids
  result = mWriteIface->assign_ids(polyhedra, mGlobalIdTag, 1);
  if (MB_SUCCESS != result) return result;

  ofile << "faces " << polygons.size() << " " << polyhedra.size() << std::endl;

  for (MBRange::iterator rit = polygons.begin(); rit != polygons.end(); rit++) {
      // get the vertices
    connecth.clear();
    result = mbImpl->get_connectivity(&(*rit), 1, connecth, true);
    if (MB_SUCCESS != result) return result;

    if (0 == connecth.size()) continue;
    
      // get the polyhedra, if any
    if (user_dimension == 3) {
      polyhedra.clear();
      result = mbImpl->get_adjacencies(MBRange(*rit, *rit), 3, false, polyhedra);
      if (MB_SUCCESS != result) return result;
    
        // put them in the connect array
      connecth.push_back((polyhedra.size() > 0 ? *polyhedra.begin() : 0));
      connecth.push_back((polyhedra.size() > 1 ? *polyhedra.rbegin() : 0));
    }
    
      // replace handles with ids
    connect.reserve(connecth.size());

      // pre-set polyhedra ids in case there aren't any
    connect[connecth.size()] = 0;
    connect[connecth.size()+1] = 0;
    result = mbImpl->tag_get_data(mGlobalIdTag, &connecth[0], 
                                  connecth.size()-2+polyhedra.size(),
                                  &connect[0]);
    if (MB_SUCCESS != result) return result;
    
      // write the data
    ofile << connecth.size()-2;
    
    for (i = 0; i < (int)connecth.size(); i++)
      ofile << " " << connect[i];

    ofile << std::endl;
  }

  ofile.close();
  
  return MB_SUCCESS;
}

