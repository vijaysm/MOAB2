#include "MBCore.hpp"
#include "MBCartVect.hpp"
#include "MBAdaptiveKDTree.hpp"
#include "MBGeomUtil.hpp"
#include "MBInternals.hpp"
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>

const int MAX_TAG_VALUE = 32; // maximum value to use when tagging entities with tree cell number
const char* TAG_NAME = "TREE_CELL";

std::string clock_to_string( clock_t t );
std::string mem_to_string( unsigned long mem );

void build_tree( MBInterface* interface, 
                 MBAdaptiveKDTree::Settings settings, 
                 unsigned meshset_flags );
void build_grid( MBInterface* iface,
                 MBAdaptiveKDTree::Settings settings,
                 unsigned meshset_flags,
                 const double dims[3] );
void delete_existing_tree( MBInterface* interface );
void print_stats( MBInterface* interface );
void tag_triangles( MBInterface* interface );
void tag_vertices( MBInterface* interface );
void write_tree_blocks( MBInterface* interface, const char* file );

static const char* ds( MBAdaptiveKDTree::CandidatePlaneSet scheme )
{
  static const char def[] = " (default)";
  const static char non[] = "";
  MBAdaptiveKDTree::Settings st;
  return scheme == st.candidatePlaneSet ? def : non;
}
  
static void usage( bool err = true )
{
  std::ostream& s = err ? std::cerr : std::cout;
  s << "kd_tree_tool [-s|-S] [-d <int>] [-n <int>] [-u|-p|-m|-v] [-N <int>] <input file> <output file>" << std::endl
    << "kd_tree_tool [-s|-S] [-d <int>] -G <dims> <output file>" << std::endl
    << "kd_tree_tool [-h]" << std::endl;
  if (!err) {
    MBAdaptiveKDTree::Settings st;
    s << "Tool to build adaptive kd-Tree from triangles" << std::endl;
    s << "  -s        Use range-based sets for tree nodes" << std::endl
      << "  -S        Use vector-based sets for tree nodes" << std::endl
      << "  -d <int>  Specify maximum depth for tree. Default: " << st.maxTreeDepth << std::endl
      << "  -n <int>  Specify maximum entities per leaf. Default: " << st.maxEntPerLeaf << std::endl
      << "  -u        Use 'SUBDIVISION' scheme for tree construction" << ds(MBAdaptiveKDTree::SUBDIVISION) << std::endl
      << "  -p        Use 'SUBDIVISION_SNAP' tree construction algorithm." << ds(MBAdaptiveKDTree::SUBDIVISION_SNAP) << std::endl
      << "  -m        Use 'VERTEX_MEDIAN' tree construction algorithm." << ds(MBAdaptiveKDTree::VERTEX_MEDIAN) << std::endl
      << "  -v        Use 'VERTEX_SAMPLE' tree construction algorithm." << ds(MBAdaptiveKDTree::VERTEX_SAMPLE) << std::endl
      << "  -N <int>  Specify candidate split planes per axis.  Default: " << st.candidateSplitsPerDir << std::endl
      << "  -t        Tag triangles will tree cell number." << std::endl
      << "  -T        Write tree boxes to file." << std::endl
      << "  -G <dims> Generate grid - no input triangles.  Dims must be " << std::endl
      << "            HxWxD or H." << std::endl
      << std::endl;
  }
  exit( err );
}

static void memory_use( unsigned long& vsize, unsigned long& rss )
{
  char buffer[512];
  int filp = open( "/proc/self/stat", O_RDONLY );
  ssize_t r = read( filp, buffer, sizeof(buffer)-1 );
  close( filp );
  if (r < 0) r = 0;
  vsize = rss = 0;
  buffer[r] = '\0';
  sscanf( buffer, "%*d %*s %*c "         // pid command state
                  "%*d %*d "             // ppid pgrp
                  "%*d %*d %*d "         // session tty_nr tpgid
                  "%*u "                // flags
                  "%*u %*u %*u %*u " // minflt cminflt majflt cmajflt
                  "%*u %*u %*d %*d " // utime stime cutime cstime
                  "%*d %*d %*d "      // priority nice (unused)
                  "%*d %*u "           // itrealval starttime
                  "%lu %lu",             &vsize, &rss );
  rss *= getpagesize();
}

static int parseint( int& i, int argc, char* argv[] )
{
  char* end;
  ++i;
  if (i == argc) {
    std::cerr << "Expected value following '" << argv[i-1] << "'" << std::endl;
    usage();
  }
  
  int result = strtol( argv[i], &end, 0 );
  if (result < 0 || *end) {
    std::cerr << "Expected positive integer following '" << argv[i-1] << "'" << std::endl;
    usage();
  }
  
  return result;
}

static void parsedims( int& i, int argc, char* argv[], double dims[3] )
{
  char* end;
  ++i;
  if (i == argc) {
    std::cerr << "Expected value following '" << argv[i-1] << "'" << std::endl;
    usage();
  }
  
  dims[0] = strtod( argv[i], &end );
  if (*end == '\0') {
    dims[2] = dims[1] = dims[0];
    return;
  }
  else if (*end != 'x' && *end != 'X')
    goto error;
  
  ++end;
  dims[1] = strtod( end, &end );
  if (*end == '\0') {
    dims[2] = dims[1];
    return;
  }
  else if (*end != 'x' && *end != 'X')
    goto error;
  
  ++end;
  dims[2] = strtod( end, &end );
  if (*end != '\0')
    goto error;
  
  return;
error:
  std::cerr << "Invaild dimension specification." << std::endl;
  usage();
}
    

int main( int argc, char* argv[] )
{
  const char* input_file = 0;
  const char* output_file = 0;
  const char* tree_file = 0;
  MBAdaptiveKDTree::Settings settings;
  unsigned meshset_flags = MESHSET_SET;
  bool tag_tris = false;
  clock_t load_time, build_time, stat_time, tag_time, write_time, block_time;
  bool make_grid = false;
  double dims[3];
  
  for (int i = 1; i < argc; ++i) {
    if (argv[i][0] != '-') {
      if (!input_file)
        input_file = argv[i];
      else if (!output_file)
        output_file = argv[i];
      else 
        usage();
      continue;
    }
    
    if (!argv[i][1] || argv[i][2])
      usage();
    
    switch (argv[i][1]) {
      case 's': meshset_flags = MESHSET_SET;                        break;
      case 'S': meshset_flags = MESHSET_ORDERED;                    break;
      case 'd': settings.maxTreeDepth  = parseint( i, argc, argv ); break;
      case 'n': settings.maxEntPerLeaf = parseint( i, argc, argv ); break;
      case 'u': settings.candidatePlaneSet = MBAdaptiveKDTree::SUBDIVISION;      break;
      case 'p': settings.candidatePlaneSet = MBAdaptiveKDTree::SUBDIVISION_SNAP; break;
      case 'm': settings.candidatePlaneSet = MBAdaptiveKDTree::VERTEX_MEDIAN;    break;
      case 'v': settings.candidatePlaneSet = MBAdaptiveKDTree::VERTEX_SAMPLE;    break;
      case 'N': settings.candidateSplitsPerDir = parseint( i, argc, argv );      break;
      case 't': tag_tris = true;                                    break;
      case 'T': tree_file = argv[++i];                              break;
      case 'G': make_grid = true; parsedims( i, argc, argv, dims ); break;
      case 'h': usage(false);
      default: usage();
    }
  }
  
  if (make_grid != !output_file)
    usage();
  
  MBErrorCode rval;
  MBCore moab_core;
  MBInterface* interface = &moab_core;
  
  if (make_grid) {
    load_time = 0;
    output_file = input_file;
    input_file = 0;
    build_time = clock();
    build_grid( interface, settings, meshset_flags, dims );
    build_time = clock() - build_time;
  }
  else {
    load_time = clock();
    rval = interface->load_mesh( input_file );
    if (MB_SUCCESS != rval) {
      std::cerr << "Error reading file: " << input_file << std::endl;
      exit(2);
    }
    load_time = clock() - load_time;


    delete_existing_tree( interface );

    std::cout << "Building tree..." << std::endl;
    build_time = clock();
    build_tree( interface, settings, meshset_flags );
    build_time = clock() - build_time;
  }
  
  std::cout << "Calculating stats..." << std::endl;
  print_stats( interface );
  stat_time = clock() - build_time;
  
  if (tag_tris) {
    std::cout << "Tagging tree..." << std::endl;
    tag_triangles( interface );
    tag_vertices( interface );
  }
  tag_time = clock() - stat_time;
  
  std::cout << "Writing file... "; std::cout.flush();
  rval = interface->write_mesh( output_file );
  if (MB_SUCCESS != rval) {
    std::cerr << "Error writing file: " << output_file << std::endl;
    exit(3);
  }
  write_time = clock() - tag_time; 
  std::cout << "Wrote " << output_file << std::endl;
  
  if (tree_file) {
    std::cout << "Writing tree block rep..."; std::cout.flush();
    write_tree_blocks( interface, tree_file );
    std::cout << "Wrote " << tree_file << std::endl;
  }
  block_time = clock() - write_time;
  
  std::cout   << "Times:  " 
              << "    Load" 
              << "   Build"
              << "   Stats"
              << "   Write";
  if (tag_tris)
    std::cout << "Tag Sets";
  if (tree_file)
    std::cout << "Block   ";
  std::cout   << std::endl;
  
  std::cout   << "        "
              << std::setw(8) << clock_to_string(load_time)
              << std::setw(8) << clock_to_string(build_time)
              << std::setw(8) << clock_to_string(stat_time)
              << std::setw(8) << clock_to_string(write_time);
  if (tag_tris)
    std::cout << std::setw(8) << clock_to_string(tag_time);
  if (tree_file)
    std::cout << std::setw(8) << clock_to_string(block_time);
  std::cout   << std::endl;

  return 0;
}

  
void delete_existing_tree( MBInterface* interface )
{
  MBRange trees;
  MBAdaptiveKDTree tool(interface);

  tool.find_all_trees( trees );
  if (!trees.empty())
    std::cout << "Destroying existing tree(s) contained in file" << std::endl;
  for (MBRange::iterator i = trees.begin(); i != trees.end(); ++i)
    tool.delete_tree( *i );
  
  trees.clear();
  tool.find_all_trees( trees );
  if (!trees.empty()) {
    std::cerr << "Failed to destroy existing trees. Aborting" << std::endl;
    exit( 5 );
  }
}
  
void build_tree( MBInterface* interface, 
                 MBAdaptiveKDTree::Settings settings, 
                 unsigned meshset_flags )
{
  MBErrorCode rval;
  MBEntityHandle root = 0;
  MBRange triangles;
  
  rval = interface->get_entities_by_type( 0, MBTRI, triangles );
  if (MB_SUCCESS != rval || triangles.empty()) {
    std::cerr << "No triangles from which to build tree." << std::endl;
    exit(4);
  }
  
  MBAdaptiveKDTree tool(interface, 0, meshset_flags);
  rval = tool.build_tree( triangles, root, &settings );
  if (MB_SUCCESS != rval || !root) {
    std::cerr << "Tree construction failed." << std::endl;
    exit(4);
  }
}  
  
void build_grid( MBInterface* interface, 
                 MBAdaptiveKDTree::Settings settings, 
                 unsigned meshset_flags,
                 const double dims[3] )
{
  MBErrorCode rval;
  MBEntityHandle root = 0;
  MBAdaptiveKDTree tool(interface, 0, meshset_flags);
  MBAdaptiveKDTreeIter iter;
  MBAdaptiveKDTree::Plane plane;

  double min[3] = { -0.5*dims[0], -0.5*dims[1], -0.5*dims[2] };
  double max[3] = {  0.5*dims[0],  0.5*dims[1],  0.5*dims[2] };
  rval = tool.create_tree( min, max, root );
  if (MB_SUCCESS != rval || !root) {
    std::cerr << "Failed to create tree root." << std::endl;
    exit(4);
  }
  
  rval = tool.get_tree_iterator( root, iter );
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to get tree iterator." << std::endl;
  }
  
  do {
    while (iter.depth() < settings.maxTreeDepth) {
      plane.norm = iter.depth() % 3;
      plane.coord = 0.5 * (iter.box_min()[plane.norm] + iter.box_max()[plane.norm]);
      rval = tool.split_leaf( iter, plane );
      if (MB_SUCCESS != rval) {
        std::cerr << "Failed to split tree leaf at depth " << iter.depth() << std::endl;
        exit(4);
      }
    }
  } while ((rval = iter.step()) == MB_SUCCESS);
  
  if (rval != MB_ENTITY_NOT_FOUND) {
    std::cerr << "Error stepping KDTree iterator." << std::endl;
    exit(4);
  }
}  

std::string clock_to_string( clock_t t )
{
  char unit[5] = "s";
  char buffer[256];
  double dt = t / (double)CLOCKS_PER_SEC;
  //if (dt > 300) {
  //  dt /= 60;
  //  strcpy( unit, "min" );
  //}
  //if (dt > 300) {
  //  dt /= 60;
  //  strcpy( unit, "hr" );
  //}
  //if (dt > 100) {
  //  dt /= 24;
  //  strcpy( unit, "days" );
  //}
  sprintf(buffer,"%0.2f%s",dt,unit);
  return buffer;
}

std::string mem_to_string( unsigned long mem )
{
  char unit[3] = "B";
  if (mem > 9*1024) {
    mem = (mem + 512) / 1024;
    strcpy( unit, "kB" );
  }
  if (mem > 9*1024) {
    mem = (mem + 512) / 1024;
    strcpy( unit, "MB" );
  }
  if (mem > 9*1024) {
    mem = (mem + 512) / 1024;
    strcpy( unit, "GB" );
  }
  char buffer[256];
  sprintf(buffer, "%lu %s", mem, unit );
  return buffer;
}

template <typename T> 
struct SimpleStat 
{
  T min, max, sum, sqr;
  size_t count;
  SimpleStat();
  void add( T value );
  double avg() const { return (double)sum / count; }
  double rms() const { return sqrt( (double)sqr / count ); }
  double dev() const { return sqrt( (count * (double)sqr - (double)sum * (double)sum) / ((double)count * (count - 1) ) ); }
};

template <typename T> SimpleStat<T>::SimpleStat()
  : min(  std::numeric_limits<T>::max() ),
    max(  std::numeric_limits<T>::min() ),
    sum( 0 ), sqr( 0 ), count( 0 )
  {}

template <typename T> void SimpleStat<T>::add( T value )
{
  if (value < min)
    min = value;
  if (value > max)
    max = value;
  sum += value;
  sqr += value*value;
  ++count;
}
  
void print_stats( MBInterface* interface )
{
  MBEntityHandle root;
  MBRange range;
  MBAdaptiveKDTree tool(interface);
  
  tool.find_all_trees( range );
  if (range.size() != 1) {
    if (range.empty())
      std::cerr << "Internal error: Failed to retreive tree." << std::endl;
    else
      std::cerr << "Internal error: Multiple tree roots." << std::endl;
    exit(5);
  }
  
  root = *range.begin();
  range.clear();

  MBRange tree_sets, triangles, verts;
  //interface->get_child_meshsets( root, tree_sets, 0 );
  interface->get_entities_by_type( 0, MBENTITYSET, tree_sets );
  tree_sets.erase( tree_sets.begin(), MBRange::lower_bound( tree_sets.begin(), tree_sets.end(), root ) );
  interface->get_entities_by_type( 0, MBTRI, triangles );
  interface->get_entities_by_type( 0, MBVERTEX, verts );
  triangles.merge( verts );
  tree_sets.insert( root );
  unsigned long set_used, set_amortized, set_store_used, set_store_amortized,
                set_tag_used, set_tag_amortized, tri_used, tri_amortized;
  interface->estimated_memory_use( tree_sets, 
                                   &set_used, &set_amortized, 
                                   &set_store_used, &set_store_amortized,
                                   0, 0, 0, 0,
                                   &set_tag_used, &set_tag_amortized );
  interface->estimated_memory_use( triangles,  &tri_used, &tri_amortized );
  
  int num_tri = 0;
  interface->get_number_entities_by_type( 0, MBTRI, num_tri );
  
  double min[3] = {0,0,0}, max[3] = {0,0,0};
  tool.get_tree_box( root, min, max );
  double diff[3] = { max[0]-min[0], max[1]-min[1], max[2] - min[2] };
  double tree_vol = diff[0]*diff[1]*diff[2];
  double tree_surf_area = 2*(diff[0]*diff[1] + diff[1]*diff[2] + diff[2]*diff[0]);
  
  SimpleStat<unsigned> depth, size;
  SimpleStat<double> vol, surf;
  
  MBAdaptiveKDTreeIter iter;
  tool.get_tree_iterator( root, iter );
  do {
    depth.add( iter.depth() );
    
    int num_leaf_tri;
    interface->get_number_entities_by_type( iter.handle(), MBTRI, num_leaf_tri );
    size.add(num_leaf_tri);
    
    const double* n = iter.box_min();
    const double* x = iter.box_max();
    double dims[3] = {x[0]-n[0], x[1]-n[1], x[2]-n[2]};
    
    double leaf_vol = dims[0]*dims[1]*dims[2];
    vol.add(leaf_vol);
    
    double area = 2.0*(dims[0]*dims[1] + dims[1]*dims[2] + dims[2]*dims[0]);
    surf.add(area);
    
  } while (MB_SUCCESS == iter.step());
  
  unsigned long real_rss, real_vsize;
  memory_use( real_vsize, real_rss );
  
  printf("------------------------------------------------------------------\n");
  printf("tree volume:      %f\n", tree_vol );
  printf("total triangles:  %d\n", num_tri );
  printf("number of leaves: %lu\n", (unsigned long)depth.count );
  printf("number of nodes:  %lu\n", (unsigned long)tree_sets.size() );
  printf("volume ratio:     %0.2f%%\n", 100*(vol.sum / tree_vol));
  printf("surface ratio:    %0.2f%%\n", 100*(surf.sum / tree_surf_area));
  printf("\nmemory:           used  amortized\n");
  printf("            ---------- ----------\n");
  printf("triangles   %10s %10s\n",mem_to_string(tri_used).c_str(), mem_to_string(tri_amortized).c_str());
  printf("sets (total)%10s %10s\n",mem_to_string(set_used).c_str(), mem_to_string(set_amortized).c_str());
  printf("sets        %10s %10s\n",mem_to_string(set_store_used).c_str(), mem_to_string(set_store_amortized).c_str());
  printf("set tags    %10s %10s\n",mem_to_string(set_tag_used).c_str(), mem_to_string(set_tag_amortized).c_str());
  printf("total real  %10s %10s\n",mem_to_string(real_rss).c_str(), mem_to_string(real_vsize).c_str());
  printf("\nleaf stats:        min        avg        rms        max    std.dev\n");
  printf("            ---------- ---------- ---------- ---------- ----------\n");
  printf("depth       %10u %10.1f %10.1f %10u %10.2f\n", depth.min, depth.avg(), depth.rms(), depth.max, depth.dev() );
  printf("triangles   %10u %10.1f %10.1f %10u %10.2f\n", size.min, size.avg(), size.rms(), size.max, size.dev() );
  printf("volume      %10.2g %10.2g %10.2g %10.2g %10.2g\n", vol.min, vol.avg(), vol.rms(), vol.max, vol.dev() );
  printf("surf. area  %10.2g %10.2g %10.2g %10.2g %10.2g\n", surf.min, surf.avg(), surf.rms(), surf.max, surf.dev() );
  printf("------------------------------------------------------------------\n");
}


static int hash_handle( MBEntityHandle handle )
{
  MBEntityID h = ID_FROM_HANDLE(handle);
  return (int)((h * 13 + 7) % MAX_TAG_VALUE) + 1;
}   

void tag_triangles( MBInterface* moab )
{
  MBEntityHandle root;
  MBRange range;
  MBAdaptiveKDTree tool(moab);
  
  tool.find_all_trees( range );
  if (range.size() != 1) {
    if (range.empty())
      std::cerr << "Internal error: Failed to retreive tree." << std::endl;
    else
      std::cerr << "Internal error: Multiple tree roots." << std::endl;
    exit(5);
  }
  
  root = *range.begin();
  range.clear();
  
  MBTag tag;
  int zero = 0;
  moab->tag_create( TAG_NAME, sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER, tag, &zero, true );
  
  MBAdaptiveKDTreeIter iter;
  tool.get_tree_iterator( root, iter );
  
  std::vector<int> tag_vals;
  do {
    range.clear();
    moab->get_entities_by_type( iter.handle(), MBTRI, range );
    tag_vals.clear();
    tag_vals.resize( range.size(), hash_handle(iter.handle()) );
    moab->tag_set_data( tag, range, &tag_vals[0] );
  } while (MB_SUCCESS == iter.step());
}

void tag_vertices( MBInterface* moab )
{
  MBEntityHandle root;
  MBRange range;
  MBAdaptiveKDTree tool(moab);
  
  tool.find_all_trees( range );
  if (range.size() != 1) {
    if (range.empty())
      std::cerr << "Internal error: Failed to retreive tree." << std::endl;
    else
      std::cerr << "Internal error: Multiple tree roots." << std::endl;
    exit(5);
  }
  
  root = *range.begin();
  range.clear();
  
  MBTag tag;
  int zero = 0;
  moab->tag_create( TAG_NAME, sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER, tag, &zero, true );
  
  MBAdaptiveKDTreeIter iter;
  tool.get_tree_iterator( root, iter );
  
  do {
    range.clear();
    moab->get_entities_by_type( iter.handle(), MBTRI, range );
    
    int tag_val = hash_handle(iter.handle());
    MBRange verts;
    moab->get_adjacencies( range, 0, false, verts, MBInterface::UNION );
    for (MBRange::iterator i = verts.begin(); i != verts.end(); ++i) {
      MBCartVect coords;
      moab->get_coords( &*i, 1, coords.array() );
      if (MBGeomUtil::box_point_overlap( MBCartVect(iter.box_min()),
                                         MBCartVect(iter.box_max()),
                                         coords, 1e-7 )) 
        moab->tag_set_data( tag, &*i, 1, &tag_val );
    }
  } while (MB_SUCCESS == iter.step());
}

void write_tree_blocks( MBInterface* interface, const char* file )
{
  MBEntityHandle root;
  MBRange range;
  MBAdaptiveKDTree tool(interface);
  
  tool.find_all_trees( range );
  if (range.size() != 1) {
    if (range.empty())
      std::cerr << "Internal error: Failed to retreive tree." << std::endl;
    else
      std::cerr << "Internal error: Multiple tree roots." << std::endl;
    exit(5);
  }
  
  root = *range.begin();
  range.clear();
  
  MBCore moab2;
  MBTag tag;
  int zero = 0;
  moab2.tag_create( TAG_NAME, sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER, tag, &zero, true );
  
  
  MBAdaptiveKDTreeIter iter;
  tool.get_tree_iterator( root, iter );
  
  do {
    double x1 = iter.box_min()[0];
    double y1 = iter.box_min()[1];
    double z1 = iter.box_min()[2];
    double x2 = iter.box_max()[0];
    double y2 = iter.box_max()[1];
    double z2 = iter.box_max()[2];
    double coords[24] = { x1, y1, z1,
                          x2, y1, z1,
                          x2, y2, z1,
                          x1, y2, z1,
                          x1, y1, z2,
                          x2, y1, z2,
                          x2, y2, z2,
                          x1, y2, z2 };
    MBEntityHandle verts[8], elem;
    for (int i = 0; i < 8; ++i)
      moab2.create_vertex( coords + 3*i, verts[i] );
    moab2.create_element( MBHEX, verts, 8, elem );
    int tag_val = hash_handle(iter.handle());
    moab2.tag_set_data( tag, &elem, 1, &tag_val );
  } while (MB_SUCCESS == iter.step());

  moab2.write_mesh( file );
}
 
