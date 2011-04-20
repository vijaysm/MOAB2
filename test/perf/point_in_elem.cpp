#include "moab/Core.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/Range.hpp"
#include "moab/CartVect.hpp"
#include "moab/GeomUtil.hpp"
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <assert.h>

#define CHK(ErrorCode) do { if (MB_SUCCESS != (ErrorCode)) fail( (ErrorCode), __FILE__, __LINE__ ); } while(false)

using namespace moab;

static void fail( ErrorCode error_code, const char* file_name, int line_number );

enum TreeType {  UseKDTree, UseNoTree, UseDefaultTree = UseKDTree };

const char* default_str = "(default)";
const char* empty_str = "";
inline const char* is_default_tree( TreeType type ) {
  return type == UseDefaultTree ? default_str : empty_str;
}

const long DEFAULT_NUM_TEST = 100000;
const long HARD_MAX_UNIQUE_POINTS = 100000;
const long HARD_MIN_UNIQUE_POINTS = 1000;
const long FRACTION_UNIQUE_POINTS = 100;

static void usage( char* argv0, bool help = false )
{
  const char* usage_str = "[-k|-v] [-n <N>] [-d <N>] [-e <N>] <input_mesh>";
  if (!help) {
    std::cerr << "Usage: " << argv0 << " " << usage_str << std::endl;
    std::cerr << "     : " << argv0 << " -h" << std::endl;
    exit(1);
  }
  
  std::cout << "Usage: " << argv0 << " " << usage_str << std::endl
            << "  -k : Use kD Tree " << is_default_tree(UseKDTree) << std::endl
            << "  -v : Use no tree" << is_default_tree(UseNoTree) << std::endl
            << "  -n : Specify number of test points (default = " << DEFAULT_NUM_TEST << ")" << std::endl
            << "  -d : Specify maximum tree depth" << std::endl
            << "  -e : Specify maximum elements per leaf" << std::endl
            << " <input_mesh> : Mesh to build and query." << std::endl
            << std::endl;
  exit(0);
}

static void generate_random_points( Interface& mesh, size_t num_points,
                                    std::vector<CartVect>& points,
                                    std::vector<EntityHandle>& point_elems );

static void do_kdtree_test( Interface& mesh, int tree_depth, int elem_per_leaf,
                            long num_test, const std::vector<CartVect>& points,
                            std::vector<EntityHandle>& point_elems,
                            clock_t& build_time, clock_t& test_time, size_t& depth );

static void do_linear_test( Interface& mesh, int tree_depth, int elem_per_leaf,
                            long num_test, const std::vector<CartVect>& points,
                            std::vector<EntityHandle>& point_elems,
                            clock_t& build_time, clock_t& test_time, size_t& depth );

int main( int argc, char* argv[] )
{
  const char* input_file = 0;
  long tree_depth = -1;
  long elem_per_leaf = -1;
  long num_points = DEFAULT_NUM_TEST;
  TreeType type = UseDefaultTree;
  char* endptr;
  
    // PARSE ARGUMENTS
    
  long* valptr = 0;
  for (int i = 1; i < argc; ++i) {
    if (valptr) {
      *valptr = strtol( argv[i], &endptr, 0 );
      if (*valptr < 1 || *endptr) {
        std::cerr << "Invalid value for '" << argv[i-1] << "' flag: " << argv[i] << std::endl;
        exit(1);
      }
      valptr = 0;
    }
    else if (argv[i][0] == '-' && argv[i][1] && !argv[i][2]) {
      switch (argv[i][1]) {
        case 'h': usage(argv[0],true); break;
        case 'k': type = UseKDTree; break;
        case 'v': type = UseNoTree; break;
        case 'd': valptr = &tree_depth; break;
        case 'e': valptr = &elem_per_leaf; break;
        case 'n': valptr = &num_points; break; 
        default:
          std::cerr << "Invalid flag: " << argv[i] << std::endl;
          usage(argv[0]);
          break;
      }
    }
    else if (!input_file) {
      input_file = argv[i];
    }
    else {
      std::cerr << "Unexpected argument: " << argv[i] << std::endl;
      usage(argv[0]);
    }
  }
  if (valptr) {
    std::cerr << "Expected value following '" << argv[argc-1] << "' flag" << std::endl;
    usage(argv[0]);
  }
  if (!input_file) {
    std::cerr << "No input file specified." << std::endl;
    usage(argv[0]);
  }
  
    // LOAD MESH
    
  Core moab;
  Interface& mb = moab;
  ErrorCode rval;
  std::string init_msg, msg;
  mb.get_last_error( init_msg );
  rval = mb.load_file( input_file );
  if (MB_SUCCESS != rval) {
    std::cerr << input_file << " : failed to read file." << std::endl;
    mb.get_last_error( msg );
    if (msg != init_msg) 
      std::cerr << "message : " << msg << std::endl;
  }
  
    // GENERATE TEST POINTS
  
  int num_unique = num_points / FRACTION_UNIQUE_POINTS;
  if (num_unique > HARD_MAX_UNIQUE_POINTS)
    num_unique = HARD_MAX_UNIQUE_POINTS;
  else if (num_unique < HARD_MIN_UNIQUE_POINTS)
    num_unique = num_points;
  std::vector<CartVect> points;
  std::vector<EntityHandle> elems;
  generate_random_points( mb, num_unique, points, elems );
  
    // GET MEMORY USE BEFORE BUILDING TREE
  
  unsigned long init_total_storage;
  mb.estimated_memory_use( 0, 0, &init_total_storage );
  
    // RUN TIMING TEST
  clock_t build_time, test_time;
  size_t actual_depth;
  std::vector<EntityHandle> results(points.size());
  switch (type) {
    case UseKDTree: 
      do_kdtree_test( mb, tree_depth, elem_per_leaf, num_points, points, results, build_time, test_time, actual_depth );
      break;
    case UseNoTree: 
      do_linear_test( mb, tree_depth, elem_per_leaf, num_points, points, results, build_time, test_time, actual_depth );
      break;
  }
  
  unsigned long fini_total_storage;
  mb.estimated_memory_use( 0, 0, &fini_total_storage );
  
    // VALIDATE RESULTS
  int fail = 0;
  if (results != elems) {
    std::cout << "WARNING:  Test produced invalid results" << std::endl;
    fail = 1;
  }
  
    // SUMMARIZE RESULTS
  std::cout << "Number of test queries: " << num_points << std::endl;
  std::cout << "Tree build time: " << ((double)build_time)/CLOCKS_PER_SEC << " seconds" << std::endl;
  std::cout << "Total query time: " << ((double)test_time)/CLOCKS_PER_SEC << " seconds" << std::endl;
  std::cout << "Time per query: " << ((double)test_time)/CLOCKS_PER_SEC/num_points << " seconds" << std::endl;
  std::cout << "Tree depth: " << actual_depth << std::endl;
  if (actual_depth)
    std::cout << "Total query time/depth: " << ((double)test_time)/CLOCKS_PER_SEC/actual_depth << " seconds" << std::endl;
  std::cout << std::endl;
  std::cout << "Bytes before tree construction: " << init_total_storage << std::endl;
  std::cout << "Bytes after tree construction: " << fini_total_storage <<std::endl;
  std::cout << "Difference: " << fini_total_storage - init_total_storage << " bytes" << std::endl;
  return fail;
}


void fail( ErrorCode error_code, const char* file, int line )
{
  std::cerr << "Internal error (error code " << error_code << ") at " << file << ":" << line << std::endl;
  abort();
}

const int HexSign[8][3] = { { -1, -1, -1 },
                            {  1, -1, -1 },
                            {  1,  1, -1 },
                            { -1,  1, -1 },
                            { -1, -1,  1 },
                            {  1, -1,  1 },
                            {  1,  1,  1 },
                            { -1,  1,  1 } };

static CartVect random_point_in_hex( Interface& mb, EntityHandle hex )
{
  const double f = RAND_MAX/2;
  CartVect xi( ((double)rand() - f)/f, 
                 ((double)rand() - f)/f, 
                 ((double)rand() - f)/f );
  CartVect coords[8];
  const EntityHandle* conn;
  int len;
  ErrorCode rval = mb.get_connectivity( hex, conn, len, true );
  if (len != 8 && MB_SUCCESS != rval) {
    std::cerr << "Invalid element" << std::endl;
    assert(false);
    abort();
  }
  rval = mb.get_coords( conn, 8, reinterpret_cast<double*>(coords) );
  CHK(rval);
  
  CartVect point(0,0,0);
  for (unsigned i = 0; i < 8; ++i) {
    double coeff = 0.125;
    for (unsigned j = 0; j < 3; ++j)
      coeff *= 1 + HexSign[i][j] * xi[j];
    point += coeff * coords[i];
  }
  
  return point;
}

void generate_random_points( Interface& mb, size_t num_points,
                             std::vector<CartVect>& points,
                             std::vector<EntityHandle>& point_elems )
{
  Range elems;
  ErrorCode rval;
  rval = mb.get_entities_by_dimension( 0, 3, elems );
  CHK(rval);
  if (!elems.all_of_type(MBHEX)) {
    std::cerr << "Warning: ignoring non-hexahedral elements." << std::endl;
    std::pair< Range::iterator, Range::iterator > p = elems.equal_range(MBHEX);
    elems.erase( p.second, elems.end() );
    elems.erase( elems.begin(), p.first );
  }
  if (elems.empty()) {
    std::cerr << "Input file contains no hexahedral elements." << std::endl;
    exit(1);
  }
  
  points.resize( num_points );
  point_elems.resize( num_points );
  const size_t num_elem = elems.size();
  for (size_t i = 0; i < num_points; ++i) {
    size_t offset = 0;
    for (size_t x = num_elem; x > 0; x /= RAND_MAX)
      offset += rand();
    offset %= num_elem;
    point_elems[i] = elems[offset];
    points[i] = random_point_in_hex( mb, point_elems[i] );
  }
}

void do_kdtree_test( Interface& mb, int tree_depth, int elem_per_leaf,
                     long num_test, const std::vector<CartVect>& points,
                     std::vector<EntityHandle>& point_elems,
                     clock_t& build_time, clock_t& test_time, size_t& depth )
{
  ErrorCode rval;
  clock_t init = clock();
  AdaptiveKDTree tool( &mb );
  EntityHandle root;
  AdaptiveKDTree::Settings settings;
  if (tree_depth > 0)
    settings.maxTreeDepth = tree_depth;
  if (elem_per_leaf > 0)
    settings.maxEntPerLeaf = elem_per_leaf;
  Range all_hexes;
  rval = mb.get_entities_by_type( 0, MBHEX, all_hexes );
  CHK(rval);
  rval = tool.build_tree( all_hexes, root, &settings );
  CHK(rval);
  all_hexes.clear();
  build_time = clock() - init;
  
  EntityHandle leaf;
  std::vector<EntityHandle> hexes;
  std::vector<EntityHandle>::iterator j;
  CartVect coords[8];
  const EntityHandle* conn;
  int len;
  for (long i = 0; i < num_test; ++i) {
    const size_t idx = (size_t)i % points.size();
    rval = tool.leaf_containing_point( root, points[idx].array(), leaf ); CHK(rval);
    hexes.clear();
    rval = mb.get_entities_by_handle( leaf, hexes ); CHK(rval);
    for (j = hexes.begin(); j != hexes.end(); ++j) {
      rval = mb.get_connectivity( *j, conn, len, true ); CHK(rval);
      rval = mb.get_coords( conn, 8, reinterpret_cast<double*>(coords) ); CHK(rval);
      if (GeomUtil::point_in_trilinear_hex( coords, points[idx], 1e-12 )) {
        point_elems[idx] = *j;
        break;
      }
    }
    if (j == hexes.end())
      point_elems[idx] = 0;
  }
  
  test_time = clock() - build_time;
  unsigned min_d, max_d;
  tool.depth( root, min_d, max_d );
  depth = max_d;
}

void do_linear_test( Interface& mb, int , int ,
                     long num_test, const std::vector<CartVect>& points,
                     std::vector<EntityHandle>& point_elems,
                     clock_t& build_time, clock_t& test_time, size_t& depth )
{
  clock_t init = clock();
  Range hexes;
  ErrorCode rval = mb.get_entities_by_type( 0, MBHEX, hexes );
  CHK(rval);
  depth = 0;
  point_elems.resize( points.size() );
  build_time = clock() - init;
  
  CartVect coords[8];
  const EntityHandle* conn;
  int len;
  for (long i = 0; i < num_test; ++i) {
    const size_t idx = (size_t)i % points.size();
    for (Range::iterator h = hexes.begin(); h != hexes.end(); ++h) {
      rval = mb.get_connectivity( *h, conn, len, true ); CHK(rval);
      rval = mb.get_coords( conn, 8, reinterpret_cast<double*>(coords) ); CHK(rval);
      if (GeomUtil::point_in_trilinear_hex( coords, points[idx], 1e-12 )) {
        point_elems[idx] = *h;
        break;
      }
    }
  }
  
  test_time = clock() - build_time;
}


  
