/** \brief This test shows how to perform local point-in-element searches with MOAB's new tree searching functionality.  
 *
 * MOAB's SpatialLocator functionality performs point-in-element searches over a local or parallel mesh.
 * SpatialLocator is flexible as to what kind of tree is used and what kind of element basis functions are 
 * used to localize elements and interpolate local fields.
 */

#include <iostream>
#include <sstream>

#include "moab/Core.hpp"
#ifdef MOAB_HAVE_MPI
#  include "moab/ParallelComm.hpp"
#endif
#include "moab/Range.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/CN.hpp"
#include "moab/SpatialLocator.hpp"

#include "nanoflann.hpp"

using namespace moab;
using namespace std;

using namespace std;
using namespace nanoflann;

void dump_mem_usage();


// And this is the "dataset to kd-tree" adaptor class:
template <typename Derived, typename coord_t>
struct PointCloudAdaptor
{
  //typedef typename Derived::coord_t coord_t;

  const Derived &obj; //!< A const ref to the data set origin

  /// The constructor that sets the data set source
  PointCloudAdaptor(const Derived &obj_) : obj(obj_) { }

  /// CRTP helper method
  inline const Derived& derived() const { return obj; }

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const { return obj.size()/3; }

  // Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
  inline coord_t kdtree_distance_euclidean(const coord_t *p1, const size_t idx_p2,size_t /*size*/) const
  {
    assert(idx_p2 < obj.size());
    size_t offset = idx_p2*3;
    const coord_t d0=p1[0]-obj[offset];
    const coord_t d1=p1[1]-obj[offset+1];
    const coord_t d2=p1[2]-obj[offset+2];
    return d0*d0+d1*d1+d2*d2;
  }

  // Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
  inline coord_t kdtree_distance_angle(const coord_t *p1, const size_t idx_p2,size_t /*size*/) const
  {
    assert(idx_p2 < obj.size());
    size_t offset = idx_p2*3;
    const coord_t d0=p1[0]-obj[offset];
    const coord_t d1=p1[1]-obj[offset+1];
    const coord_t d2=p1[2]-obj[offset+2];
    return atan2((d0*d0+d1*d1+d2*d2),(std::sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2])*std::sqrt(obj[offset]*obj[offset]+obj[offset+1]*obj[offset+1]+obj[offset+2]*obj[offset+2])))*180/PI;
  }

  inline coord_t kdtree_distance(const coord_t *p1, const size_t idx_p2, size_t size) const
  {
    return kdtree_distance_euclidean(p1,idx_p2,size);
    // return kdtree_distance_angle(p1,idx_p2,size);
  }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an immediate value, the
  //  "if/else's" are actually solved at compile time.
  inline coord_t kdtree_get_pt(const size_t idx, int dim) const
  {
    assert(idx < obj.size());
    if (dim==0) return obj[idx*3];
    else if (dim==1) return obj[idx*3+1];
    else return obj[idx*3+2];
  }

  // Optional bounding-box computation: return false to default to a standard bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
  //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }

}; // end of PointCloudAdaptor


template <typename T>
void generateRandomPointCloud(std::vector<T> &point, const size_t N, const T max_range = 10)
{
  std::cout << "Generating "<< N << " point cloud...";
  point.resize(N*3);
  for (size_t i=0;i<N*3;i+=3)
  {
    point[i]   = max_range * (rand() % 1000) / T(1000);
    point[i+1] = max_range * (rand() % 1000) / T(1000);
    point[i+2] = max_range * (rand() % 1000) / T(1000);
  }

  std::cout << "done\n";
}

template <typename num_t>
ErrorCode moab_get_coordinates(Core* mb, std::vector<num_t>& cloud)
{
  ErrorCode rval;
  // Get all 3d elements in the file
  Range verts;
  std::vector<double> coords;
  rval = mb->get_entities_by_dimension(0, 3, verts);MB_CHK_SET_ERR(rval, "Error getting 3d elements");
  //rval = mb->get_entities_by_dimension(0, 0, verts);MB_CHK_SET_ERR(rval, "Error getting vertices");
  cloud.resize(verts.size()*3);
  rval = mb->get_coords(verts, &cloud[0]);MB_CHK_SET_ERR(rval, "Error getting coordinates");
  return MB_SUCCESS;
}

template <>
ErrorCode moab_get_coordinates<float>(Core* mb, std::vector<float>& cloud)
{
  ErrorCode rval;
  // Get all 3d elements in the file
  Range verts;
  std::vector<double> coords;
  rval = mb->get_entities_by_dimension(0, 0, verts);MB_CHK_SET_ERR(rval, "Error getting vertices");
  cloud.resize(verts.size()*3);
  coords.resize(verts.size()*3);
  rval = mb->get_coords(verts, &coords[0]);MB_CHK_SET_ERR(rval, "Error getting coordinates");
  for (size_t i=0;i<verts.size()*3;i++)
    cloud[i]   = static_cast<float>(coords[i]);
  return MB_SUCCESS;
}

template <typename num_t>
ErrorCode kdtree_demo(const size_t N, Core* mb=NULL, moab::ParallelComm* pcomm=NULL, bool use_reference=false)
{
  std::vector<num_t> cloud;
  ErrorCode rval;

  // Generate points:
  if (mb == NULL) {
    generateRandomPointCloud(cloud, N);
  }
  else {
    // Get all 3d elements in the file
    rval = moab_get_coordinates(mb, cloud);MB_CHK_SET_ERR(rval, "Error getting coordinates");
  }
  
  int rank = pcomm->rank();
  num_t query_pt[3] = { 0.5, 0.5, 0.5};
  //num_t query_pt[3] = { 1.45, 1.25, 1.65};

  typedef PointCloudAdaptor<std::vector<num_t>, num_t > PC2KD;
  if (use_reference)
  {
    Range elems;
    rval = mb->get_entities_by_dimension(0, 3, elems);MB_CHK_SET_ERR(rval, "Error getting vertices");
    // Code for MB KDtree here
    // Create a tree to use for the location service
    AdaptiveKDTree tree(mb);

    // Build the SpatialLocator
    SpatialLocator sl(mb, elems, &tree);
    
    // Get the box extents
    CartVect box_extents;
    BoundBox box = sl.local_box();
    box_extents = box.bMax - box.bMin;

    // Query at random places in the tree
    CartVect params;
    int is_inside = 1;
    int num_inside = 0;
    EntityHandle elem;
    bool multiple_leaves = true;
    const double search_radius = 0.25;
    std::vector<EntityHandle> ptelems;
    std::vector<double> distout;
    for (int i = 0; i < 1/*N*/; i++) {
      // pos = box.bMin + CartVect(box_extents[0] * .01 * (rand() % 100), box_extents[1] * .01 * (rand() % 100),
      //     box_extents[2] * .01 * (rand() % 100));
      // CartVect pos(query_pt);
      // rval = sl.locate_point(pos.array(), elem, params.array(), &is_inside, 0.0, 0.0);MB_CHK_ERR(rval);
      // if (is_inside) num_inside++;
      rval = tree.point_search(query_pt, elem, 1e-10, 1e-8, &multiple_leaves);MB_CHK_ERR(rval);
      mb->list_entity(elem);
    MPI_Barrier( MPI_COMM_WORLD ) ;
      cout << "["<<rank<<"] MB::knnSearch(): num_results=1\n";
      cout << "["<<rank<<"] ret_index["<< elems.index(elem) << "]" <<  " point = (" << cloud[elems.index(elem)*3] << ", " << cloud[elems.index(elem)*3+1] << ", " << cloud[elems.index(elem)*3+2] << endl;
      cout << "\n";
      rval = tree.distance_search(query_pt, search_radius, ptelems, 1e-10, 1e-8, &distout);MB_CHK_ERR(rval);
      cout << "["<<rank<<"] MB::knnSearch(): num_results=" << ptelems.size() << "\n";
      for (size_t j=0;j<ptelems.size();j++)
        cout << "["<<rank<<"] idx["<< j << "]=" << elems.index(ptelems[j]) << " point = (" << cloud[elems.index(ptelems[j])*3] << ", " << cloud[elems.index(ptelems[j])*3+1] << ", " << cloud[elems.index(ptelems[j])*3+2] << ") dist["<< j << "]=" << distout[j] << endl;
      cout << "\n";
    }

    MPI_Barrier( MPI_COMM_WORLD ) ;
    cout << "["<<rank<<"] Mesh contains " << elems.size() << " elements of type "
              << CN::EntityTypeName(mb->type_from_handle(*elems.begin())) << endl;
    cout << "["<<rank<<"] Bounding box min-max = (" << box.bMin[0] << "," << box.bMin[1] << "," << box.bMin[2] << ")-("
              << box.bMax[0] << "," << box.bMax[1] << "," << box.bMax[2] << ")" << endl;
    cout << "["<<rank<<"] Queries inside box = " << num_inside << "/" << N << " = "
              << 100.0*((double)num_inside) / N << "%" << endl;
  }
  else {
    const PC2KD  pc2kd(cloud); // The adaptor

    // construct a kd-tree index:
    typedef KDTreeSingleIndexAdaptor<
      L2_Simple_Adaptor<num_t, PC2KD > ,
      PC2KD,
      3 /* dim */
      > my_kd_tree_t;

      dump_mem_usage();

      my_kd_tree_t   index(3 /*dim*/, pc2kd, KDTreeSingleIndexAdaptorParams(30 /* max leaf */) );
      index.buildIndex();
      dump_mem_usage();

    // ----------------------------------------------------------------
    // findNeighbors():  Perform a distance-based neighbor search
    // ----------------------------------------------------------------
    {
      const size_t num_results = 1;
      size_t ret_index;
      num_t out_dist_sqr;
      nanoflann::KNNResultSet<num_t> resultSet(num_results);
      resultSet.init(&ret_index, &out_dist_sqr );
      index.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10)); 

      MPI_Barrier( MPI_COMM_WORLD ) ;
      if (!rank) {
        cout << "["<<rank<<"] knnSearch(nn="<<num_results<<"): \n";
        cout << "["<<rank<<"] ret_index["<< ret_index << "]" <<  " point = (" << cloud[ret_index*3] << ", " << cloud[ret_index*3+1] << ", " << cloud[ret_index*3+2] << ") dist=" << out_dist_sqr << endl;
      }
      MPI_Barrier( MPI_COMM_WORLD ) ;
      if (rank) {
        cout << "["<<rank<<"] knnSearch(nn="<<num_results<<"): \n";
        cout << "["<<rank<<"] ret_index["<< ret_index << "]" <<  " point = (" << cloud[ret_index*3] << ", " << cloud[ret_index*3+1] << ", " << cloud[ret_index*3+2] << ") dist=" << out_dist_sqr << endl;
      }
    }

    // ----------------------------------------------------------------
    // knnSearch():  Perform a search for the N closest points
    // ----------------------------------------------------------------
    {
      const size_t num_results = 5;
      std::vector<size_t>   ret_index(num_results);
      std::vector<num_t> out_dist_sqr(num_results);
      index.knnSearch(&query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);

      MPI_Barrier( MPI_COMM_WORLD ) ;
      if (!rank) {
        cout << "["<<rank<<"] knnSearch(): num_results=" << num_results << "\n";
        for (size_t i=0;i<num_results;i++)
          cout << "["<<rank<<"] idx["<< i << "]=" << ret_index[i] << " point = (" << cloud[ret_index[i]*3] << ", " << cloud[ret_index[i]*3+1] << ", " << cloud[ret_index[i]*3+2] << ") dist["<< i << "]=" << out_dist_sqr[i] << endl;
        cout << "\n";
      }
      MPI_Barrier( MPI_COMM_WORLD ) ;
      if (rank) {
        cout << "["<<rank<<"] knnSearch(): num_results=" << num_results << "\n";
        for (size_t i=0;i<num_results;i++)
          cout << "["<<rank<<"] idx["<< i << "]=" << ret_index[i] << " point = (" << cloud[ret_index[i]*3] << ", " << cloud[ret_index[i]*3+1] << ", " << cloud[ret_index[i]*3+2] << ") dist["<< i << "]=" << out_dist_sqr[i] << endl;
        cout << "\n";
      }
    }

    // ----------------------------------------------------------------
    // radiusSearch():  Perform a search for the N closest points
    // ----------------------------------------------------------------
    {
      const num_t search_radius = static_cast<num_t>(0.25);
      std::vector<std::pair<size_t,num_t> >   ret_matches;

      nanoflann::SearchParams params;
      params.sorted = true;

      const size_t nMatches = index.radiusSearch(&query_pt[0], search_radius, ret_matches, params);

      MPI_Barrier( MPI_COMM_WORLD ) ;
      if (!rank) {
        cout << "["<<rank<<"] radiusSearch(): radius=" << search_radius << " -> " << nMatches << " matches\n";
        for (size_t i=0;i<nMatches;i++)
          cout << "["<<rank<<"] idx["<< i << "]=" << ret_matches[i].first << " point = (" << cloud[ret_matches[i].first*3] << ", " << cloud[ret_matches[i].first*3+1] << ", " << cloud[ret_matches[i].first*3+2] << ") dist["<< i << "]=" << ret_matches[i].second << endl;
        cout << "\n";
      }
      MPI_Barrier( MPI_COMM_WORLD ) ;
      if (rank) {
        cout << "["<<rank<<"] radiusSearch(): radius=" << search_radius << " -> " << nMatches << " matches\n";
        for (size_t i=0;i<nMatches;i++)
          cout << "["<<rank<<"] idx["<< i << "]=" << ret_matches[i].first << " point = (" << cloud[ret_matches[i].first*3] << ", " << cloud[ret_matches[i].first*3+1] << ", " << cloud[ret_matches[i].first*3+2] << ") dist["<< i << "]=" << ret_matches[i].second << endl;
        cout << "\n";
      }
    }

    if (mb)
    {
      // ADding a new point
      double newvert[3] = {0.55, 0.55, 0.55};
      EntityHandle newverthandle;
      mb->create_vertex(newvert, newverthandle);
      cloud.push_back(newvert[0]);
      cloud.push_back(newvert[1]);
      cloud.push_back(newvert[2]);
      index.buildIndex();
      const num_t search_radius = static_cast<num_t>(0.25);
      std::vector<std::pair<size_t,num_t> >   ret_matches;

      nanoflann::SearchParams params;
      params.sorted = true;

      const size_t nMatches = index.radiusSearch(&query_pt[0], search_radius, ret_matches, params);

      MPI_Barrier( MPI_COMM_WORLD ) ;
      if (!rank) {
        cout << "["<<rank<<"] radiusSearch(): radius=" << search_radius << " -> " << nMatches << " matches\n";
        for (size_t i=0;i<nMatches;i++)
          cout << "["<<rank<<"] idx["<< i << "]=" << ret_matches[i].first << " point = (" << cloud[ret_matches[i].first*3] << ", " << cloud[ret_matches[i].first*3+1] << ", " << cloud[ret_matches[i].first*3+2] << ") dist["<< i << "]=" << ret_matches[i].second << endl;
        cout << "\n";
      }
      MPI_Barrier( MPI_COMM_WORLD ) ;
      if (rank) {
        cout << "["<<rank<<"] radiusSearch(): radius=" << search_radius << " -> " << nMatches << " matches\n";
        for (size_t i=0;i<nMatches;i++)
          cout << "["<<rank<<"] idx["<< i << "]=" << ret_matches[i].first << " point = (" << cloud[ret_matches[i].first*3] << ", " << cloud[ret_matches[i].first*3+1] << ", " << cloud[ret_matches[i].first*3+2] << ") dist["<< i << "]=" << ret_matches[i].second << endl;
        cout << "\n";
      }
    }

  }
  

  return MB_SUCCESS;
}

void dump_mem_usage()
{
  FILE* f=fopen("/proc/self/statm","rt");
  if (!f) return;
  char str[300];
  size_t n=fread(str,1,200,f);
  str[n]=0;
  printf("MEM: %s\n",str);
  fclose(f);
}

#ifdef MOAB_HAVE_HDF5
string test_file_name = string(MESH_DIR) + string("/64bricks_512hex_256part.h5m");
#else
string test_file_name = string(MESH_DIR) + string("/mbtest1.vtk");
#endif
int main(int argc, char **argv)
{
  ErrorCode rval;
  const int nghostrings = 1;
  const int num_dim = 3;
  int num_queries = 100;
  std::stringstream sstr;

  if (argc > 3) {
    cout << "Usage: " << argv[0] << " <filename> [num_queries]" << endl;
    return 0;
  }
  else if (argc == 3) {
    test_file_name = argv[1];
    num_queries = atoi(argv[2]);
  }
  else if (argc == 2) {
    num_queries = atoi(argv[1]);
  }
  else {
    num_queries = 100;
  }

  // Instantiate
  Core mb;

  // get the ParallelComm instance
  ParallelComm *pcomm = new ParallelComm(&mb, MPI_COMM_WORLD);
  
  string roptions,woptions;
  if (pcomm->size() > 1) { // if reading in parallel, need to tell it how
    sstr.str("");
    sstr << "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=" 
         << num_dim << ".0." << nghostrings
         << ";DEBUG_IO=0;DEBUG_PIO=0";
    roptions = sstr.str();
    woptions = "PARALLEL=WRITE_PART";
  }

  moab::EntityHandle fileset;
  rval = mb.create_meshset(MESHSET_SET, fileset); MB_CHK_SET_ERR(rval, "Error creating meshset");
  // Load the file
  rval = mb.load_file(test_file_name.c_str(), &fileset, roptions.c_str());MB_CHK_SET_ERR(rval, "Error loading file");

  // NanoFLANN based queries  
  //rval = kdtree_demo<float>(num_queries);MB_CHK_SET_ERR(rval, "Error with Kdtree demo with float type");
  //rval = kdtree_demo<double>(num_queries);MB_CHK_SET_ERR(rval, "Error with Kdtree demo with double type");

  rval = kdtree_demo<double>(num_queries, &mb, pcomm, false);MB_CHK_SET_ERR(rval, "Error with Kdtree demo with MOAB Instance");
  // MOAB AdaptiveKDTree based queries
  //rval = kdtree_demo<double>(num_queries, &mb, true);MB_CHK_SET_ERR(rval, "Error with Kdtree demo with MOAB Instance");

  return 0;
}
