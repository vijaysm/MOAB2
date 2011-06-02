#include "moab/Core.hpp"
#include "moab/CartVect.hpp"
#include "OrientedBox.hpp"
#include "moab/OrientedBoxTreeTool.hpp"
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <time.h>
#include <signal.h>
#include <assert.h>

const int NUM_RAYS = 40000;
const int NUM_XSCT = 20000;

using namespace moab;

static void usage( )
{
  std::cerr << "obb_time [-r <int>] [-i <int>] <filename>" << std::endl
      << "  -r - Specify total rays to fire." << std::endl
      << "       Zero implies unbounded. Default: " << NUM_RAYS << std::endl
      << "  -i - Specify total intersecting rays to fire." << std::endl
      << "       Zero implies unbounded. Default: " << NUM_XSCT << std::endl
      << "  -s - Use set-based tree." << std::endl
      << "  -p - Measure and report traversal performance statistics" << std::endl
      << "  The input file should be generated using the '-s'" << std::endl
      << "  option with 'obb_test'" << std::endl;
  exit(1);
}

void generate_ray( const CartVect& sphere_center,
                   double sphere_radius,
                   CartVect& point,
                   CartVect& dir )
{
  const int H = RAND_MAX/2;
  point[0] = (double)rand()/H - 1;
  point[1] = (double)rand()/H - 1;
  point[2] = (double)rand()/H - 1;
  point *= sphere_radius;
  
  dir[0] = (double)rand() * -point[0];
  dir[1] = (double)rand() * -point[1];
  dir[2] = (double)rand() * -point[2];
  dir.normalize();

  point += sphere_center;
}

ErrorCode read_tree( Interface* instance,
                       const char* filename,
                       EntityHandle& tree_root_out )
{
  ErrorCode rval = instance->load_mesh( filename );
  if (MB_SUCCESS != rval)
    return rval;
  
  Tag tag;
  rval = instance->tag_get_handle( "OBB_ROOT", 1, MB_TYPE_HANDLE, tag );
  if (MB_SUCCESS != rval)
    return rval;
  
  const EntityHandle root = 0;
  return instance->tag_get_data( tag, &root, 1, &tree_root_out );
}

// global variables for CLI options
int num_rays = NUM_RAYS;
int num_xsct = NUM_XSCT;
const char* filename = 0;
bool do_sets = false;
bool do_trv_stats = false;

// global to make accessable to signal handler
int rays = 0, xsct = 0, gen = 0;
clock_t t;

extern "C" {
  void signal_handler( int ) {
    t = clock() - t;
    std::cout << filename << ":" << std::endl
              << rays << " of " << num_rays << " ray fires done" << std::endl
              << xsct << " of " << num_xsct << " intersecting fires" << std::endl
              << gen  << " unique rays used" << std::endl
              << (double)t/CLOCKS_PER_SEC << " seconds" << std::endl;
    exit(1);
  }
}

int main( int argc, char* argv[] )
{
  signal( SIGINT, &signal_handler );
  
  for (int i = 1; i < argc; ++i)
  {
    if (!strcmp( argv[i], "-r")) {
      ++i;
      if (i == argc || !argv[i][0]) {
        std::cerr << "Expected value following '-r'" << std::endl;
        usage();
      }
      char* end;
      long t = strtol( argv[i], &end, 0 );
      num_rays = (int)t;
      if (*end || t < 0 || num_rays != t) {
        std::cerr << "Expected positive integer following '-r'" << std::endl;
        usage();
      }
    }
    else if (!strcmp( argv[i], "-i")) {
      ++i;
      if (i == argc || !argv[i][0]) {
        std::cerr << "Expected value following '-i'" << std::endl;
        usage();
      }
      char* end;
      long t = strtol( argv[i], &end, 0 );
      num_xsct = (int)t;
      if (*end || t < 0 || num_xsct != t) {
        std::cerr << "Expected positive integer following '-i'" << std::endl;
        usage();
      }
    }
    else if (!strcmp( argv[i], "-s")) {
      do_sets = true;
    }
    else if (!strcmp( argv[i], "-p")) {
      do_trv_stats = true;
    }
    else if (filename) {
      std::cerr << "Invalid options or multiple file names specified." << std::endl;
        usage();
    }
    else {
      filename = argv[i];
    }
  }
  if (!filename) {
    std::cerr << "No file name specified." << std::endl;
    usage();
  }
    
  Core instance;
  Interface* iface = &instance;
  EntityHandle root;
  ErrorCode rval = read_tree( iface, filename, root );
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to read \"" <<filename<<'"'<<std::endl;
    return 2;
  }
  
  OrientedBoxTreeTool tool(iface);
  OrientedBox box;
  rval = tool.box( root, box );
  if (MB_SUCCESS != rval) {
    std::cerr << "Corrupt tree.  Cannot get box for root node." << std::endl;
    return 3;
  }

  OrientedBoxTreeTool::TrvStats* stats = NULL;
  if( do_trv_stats ){
    stats = new OrientedBoxTreeTool::TrvStats;
  }
  
  const unsigned cached = 1000;
  std::vector<double> intersections;
  std::vector<EntityHandle> sets, facets;
  CartVect point, dir;
  std::vector<CartVect> randrays;
  randrays.reserve( cached );
  int cached_idx = 0;
  
  t = clock();
  for (;;) {
    if (!num_rays) {
      if (xsct >= num_xsct)
        break;
    }
    else if (!num_xsct) {
      if (rays >= num_rays)
        break;
    }
    else if (rays >= num_rays && xsct >= num_xsct)
      break;
    
    ++rays;
    CartVect point, dir;
    if (randrays.size() < cached) {
      generate_ray( box.center, box.outer_radius(), point, dir );
      ++gen;
    }
    else {
      point = randrays[cached_idx++];
      dir   = randrays[cached_idx++];
      cached_idx = cached_idx % randrays.size();
    }
    
    intersections.clear();
    if (do_sets) {
      sets.clear();
      facets.clear();
      rval = tool.ray_intersect_sets( intersections, sets, facets, root, 1e-6, 1, point.array(), dir.array(), 0,  stats );
    }
    else {
      rval = tool.ray_intersect_triangles( intersections, facets, root, 1e-6, point.array(), dir.array(), 0, stats );
    }
    if (MB_SUCCESS != rval) {
      std::cerr << "Rayfire #" << rays << " failed." << std::endl;
      return 4;
    }
    
    if (!intersections.empty()) {
      ++xsct;
    }
    
    if (randrays.size() < cached && 
        (!intersections.empty() || !num_xsct || xsct >= num_xsct)) {
      randrays.push_back( point );
      randrays.push_back( dir );
    }
  }

  t = clock() - t;
  std::cout << rays << " ray fires done" << std::endl
            << gen  << " unique rays used" << std::endl
            << xsct << " intersecting fires" << std::endl
            << (double)t/CLOCKS_PER_SEC << " seconds" << std::endl;
  
  if( do_trv_stats ){
    std::cout << "Traversal statistics: " << std::endl;
    stats->print( std::cout );
  }

  return 0;
}




  

