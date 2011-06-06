#define IS_BUILDING_MB
#include "moab/Core.hpp"
#include "moab/CartVect.hpp"
#include "moab/OrientedBoxTreeTool.hpp"
#include "OrientedBox.hpp"
#include "Internals.hpp"
#include "moab/Range.hpp"

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <limits>
#include <stdlib.h>
#include <time.h>
#if !defined(_MSC_VER) && !defined(__MINGW32__)
#  include <unistd.h>
#  include <sys/stat.h>
#  include <fcntl.h>
#endif

using namespace moab;

std::string clock_to_string( clock_t t );
std::string mem_to_string( unsigned long mem );

const int MAX_TAG_VALUE = 20;
const char* const TAG_NAME = "OBB_ID";
const char* const TREE_TAG = "OBB_ROOT";
const char* root_tag = TREE_TAG;

ErrorCode get_root( Interface* moab, EntityHandle& root );
EntityHandle build_tree( Interface* interface, 
                           OrientedBoxTreeTool::Settings settings );
void delete_existing_tree( Interface* interface );
void print_stats( Interface* interface );
void tag_triangles( Interface* interface );
void tag_vertices( Interface* interface );
void write_tree_blocks( Interface* interface, const char* file );

static void usage( bool err = true )
{
  std::ostream& s = err ? std::cerr : std::cout;
  s << "obb_tree_tool [-s|-S] [-d <int>] [-n <int>] <input file> <output file>" << std::endl
    << "obb_tree_tool [-h]" << std::endl;
  if (!err) {
    OrientedBoxTreeTool::Settings st;
    s << "Tool to build adaptive kd-Tree from triangles" << std::endl;
    s << "  -s        Use range-based sets for tree nodes" << std::endl
      << "  -S        Use vector-based sets for tree nodes" << std::endl
      << "  -d <int>  Specify maximum depth for tree. Default: " << st.max_depth << std::endl
      << "  -n <int>  Specify maximum entities per leaf. Default: " << st.max_leaf_entities << std::endl
      << "  -m <real> Specify worst split ratio. Default: " << st.worst_split_ratio << std::endl
      << "  -M <real> Specify best split ratio. Default: " << st.best_split_ratio << std::endl
      << "  -t        Tag triangles will tree cell number." << std::endl
      << "  -T        Write tree boxes to file." << std::endl
      << "  -N        Specify mesh tag containing tree root. Default: \"" << TREE_TAG << '"' << std::endl
      << std::endl;
  }
  exit( err );
}

#if defined(_MSC_VER)  || defined(__MINGW32__)
static void memory_use( unsigned long& vsize, unsigned long& rss )
  { vsize = rss = 0; }
#else
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
#endif

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

static double parsedouble( int& i, int argc, char* argv[] )
{
  char* end;
  ++i;
  if (i == argc) {
    std::cerr << "Expected value following '" << argv[i-1] << "'" << std::endl;
    usage();
  }
  
  double result = strtod( argv[i], &end );
  if (result < 0 || *end) {
    std::cerr << "Expected positive real number following '" << argv[i-1] << "'" << std::endl;
    usage();
  }
  
  return result;
}
    

int main( int argc, char* argv[] )
{
  const char* input_file = 0;
  const char* output_file = 0;
  const char* tree_file = 0;
  OrientedBoxTreeTool::Settings settings;
  bool tag_tris = false;
  clock_t load_time, build_time, stat_time, tag_time, write_time, block_time;
  
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
      case 's': settings.set_options = MESHSET_SET;                        break;
      case 'S': settings.set_options = MESHSET_ORDERED;                    break;
      case 'd': settings.max_depth  = parseint( i, argc, argv );           break;
      case 'n': settings.max_leaf_entities = parseint( i, argc, argv );    break;
      case 'm': settings.worst_split_ratio = parsedouble( i, argc, argv ); break;
      case 'M': settings.best_split_ratio = parsedouble( i, argc, argv );  break;
      case 't': tag_tris = true;                                           break;
      case 'T': if (++i == argc) usage(); tree_file = argv[i];             break;
      case 'N': if (++i == argc) usage(); root_tag = argv[i];              break;
      case 'h': usage(false);
      default: usage();
    }
  }
  
  if (!output_file)
    usage();
  
  ErrorCode rval;
  Core moab_core;
  Interface* interface = &moab_core;
  
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
  build_tree( interface, settings );
  build_time = clock() - build_time;
  
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

ErrorCode get_root( Interface* moab, EntityHandle& root )
{
  Tag tag;
  ErrorCode rval;

  rval = moab->tag_get_handle( root_tag, 1, MB_TYPE_HANDLE, tag );
  if (MB_SUCCESS != rval)
    return rval;
 
  const EntityHandle mesh = 0;
  return moab->tag_get_data( tag, &mesh, 1, &root );
}

  
void delete_existing_tree( Interface* interface )
{
  EntityHandle root;
  ErrorCode rval = get_root(interface, root);
  if (MB_SUCCESS == rval) {
    OrientedBoxTreeTool tool(interface);
    rval = tool.delete_tree( root );
    if (MB_SUCCESS != rval) {
      std::cerr << "Failed to destroy existing trees. Aborting" << std::endl;
      exit( 5 );
    }
  }
}
  
EntityHandle build_tree( Interface* interface, OrientedBoxTreeTool::Settings settings )
{
  ErrorCode rval;
  EntityHandle root = 0;
  Range triangles;
  
  rval = interface->get_entities_by_type( 0, MBTRI, triangles );
  if (MB_SUCCESS != rval || triangles.empty()) {
    std::cerr << "No triangles from which to build tree." << std::endl;
    exit(4);
  }
  
  OrientedBoxTreeTool tool( interface );
  rval = tool.build( triangles, root, &settings );
  if (MB_SUCCESS != rval || !root) {
    std::cerr << "Tree construction failed." << std::endl;
    exit(4);
  }
  
    // store tree root
  Tag roottag;
  rval = interface->tag_get_handle( root_tag, 1, MB_TYPE_HANDLE, roottag, MB_TAG_CREAT|MB_TAG_SPARSE );
  if (MB_SUCCESS != rval) {
    std::cout << "Failed to create root tag: \"" << root_tag << '"' << std::endl;
    exit(2);
  }
  const EntityHandle mesh = 0;
  rval = interface->tag_set_data( roottag, &mesh, 1, &root );
  if (MB_SUCCESS != rval) {
    std::cout << "Failed to set root tag: \"" << root_tag << '"' << std::endl;
    exit(2);
  }
  
  return root;
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

void print_stats( Interface* interface )
{
  EntityHandle root;
  Range range;
  get_root( interface, root );
  OrientedBoxTreeTool tool(interface);

  Range tree_sets, triangles, verts;
  //interface->get_child_meshsets( root, tree_sets, 0 );
  interface->get_entities_by_type( 0, MBENTITYSET, tree_sets );
  tree_sets.erase( tree_sets.begin(), Range::lower_bound( tree_sets.begin(), tree_sets.end(), root ) );
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
  
  tool.stats( root, std::cout );
  
  unsigned long real_rss, real_vsize;
  memory_use( real_vsize, real_rss );
  
  printf("------------------------------------------------------------------\n");
  printf("\nmemory:           used  amortized\n");
  printf("            ---------- ----------\n");
  printf("triangles   %10s %10s\n",mem_to_string(tri_used).c_str(), mem_to_string(tri_amortized).c_str());
  printf("sets (total)%10s %10s\n",mem_to_string(set_used).c_str(), mem_to_string(set_amortized).c_str());
  printf("sets        %10s %10s\n",mem_to_string(set_store_used).c_str(), mem_to_string(set_store_amortized).c_str());
  printf("set tags    %10s %10s\n",mem_to_string(set_tag_used).c_str(), mem_to_string(set_tag_amortized).c_str());
  printf("total real  %10s %10s\n",mem_to_string(real_rss).c_str(), mem_to_string(real_vsize).c_str());
  printf("------------------------------------------------------------------\n");
}


static int hash_handle( EntityHandle handle )
{
  EntityID h = ID_FROM_HANDLE(handle);
  return (int)((h * 13 + 7) % MAX_TAG_VALUE) + 1;
}   

class TriTagger : public OrientedBoxTreeTool::Op
{
private:
  Interface* mMB;
  Tag mTag;
  std::vector<EntityHandle> mHandles;
  std::vector<int> mTagData;
public:
  TriTagger( Tag tag, Interface* moab )
    : mMB(moab), mTag(tag) {}

  ErrorCode visit( EntityHandle, int, bool& descent )
    { descent = true; return MB_SUCCESS; }
  
  ErrorCode leaf( EntityHandle node ) {
    mHandles.clear();
    mMB->get_entities_by_handle( node, mHandles );
    mTagData.clear();
    mTagData.resize( mHandles.size(), hash_handle( node ) );
    mMB->tag_set_data( mTag, &mHandles[0], mHandles.size(), &mTagData[0] );
    return MB_SUCCESS;
  }
};
    

void tag_triangles( Interface* moab )
{
  EntityHandle root;
  ErrorCode rval = get_root( moab, root );
  if (MB_SUCCESS != rval) {
    std::cerr << "Internal error: Failed to retreive tree." << std::endl;
    exit(5);
  }

  Tag tag;
  int zero = 0;
  moab->tag_get_handle( TAG_NAME, 1, MB_TYPE_INTEGER, tag, MB_TAG_DENSE|MB_TAG_CREAT, &zero );
  TriTagger op( tag, moab );
  
  OrientedBoxTreeTool tool(moab);
  rval = tool.preorder_traverse( root, op );
  if (MB_SUCCESS != rval) {
    std::cerr << "Internal error tagging triangles" << std::endl;
    exit(5);
  }
}


class VtxTagger : public OrientedBoxTreeTool::Op
{
private:
  Interface* mMB;
  Tag mTag;
  std::vector<EntityHandle> mHandles;
  std::vector<EntityHandle> mConn;
  std::vector<int> mTagData;
public:
  VtxTagger( Tag tag, Interface* moab )
    : mMB(moab), mTag(tag) {}

  ErrorCode visit( EntityHandle, int, bool& descent )
    { descent = true; return MB_SUCCESS; }
  
  ErrorCode leaf( EntityHandle node ) {
    mHandles.clear();
    mMB->get_entities_by_handle( node, mHandles );
    mConn.clear();
    mMB->get_connectivity( &mHandles[0], mHandles.size(), mConn );
    mTagData.clear();
    mTagData.resize( mConn.size(), hash_handle( node ) );
    mMB->tag_set_data( mTag, &mConn[0], mConn.size(), &mTagData[0] );
    return MB_SUCCESS;
  }
};
    

void tag_vertices( Interface* moab )
{
  EntityHandle root;
  ErrorCode rval = get_root( moab, root );
  if (MB_SUCCESS != rval) {
    std::cerr << "Internal error: Failed to retreive tree." << std::endl;
    exit(5);
  }

  Tag tag;
  int zero = 0;
  moab->tag_get_handle( TAG_NAME, 1, MB_TYPE_INTEGER, tag, MB_TAG_DENSE|MB_TAG_CREAT, &zero );
  VtxTagger op( tag, moab );
  
  OrientedBoxTreeTool tool(moab);
  rval = tool.preorder_traverse( root, op );
  if (MB_SUCCESS != rval) {
    std::cerr << "Internal error tagging vertices" << std::endl;
    exit(5);
  }
}


class LeafHexer : public OrientedBoxTreeTool::Op
{
private:
  OrientedBoxTreeTool* mTool;
  Interface* mOut;
  Tag mTag;
  std::vector<EntityHandle> mHandles;
  std::vector<EntityHandle> mConn;
  std::vector<int> mTagData;
public:
  LeafHexer( OrientedBoxTreeTool* tool, Interface* mb2, Tag tag )
    : mTool(tool), mOut( mb2 ), mTag(tag) {}

  ErrorCode visit( EntityHandle, int, bool& descent )
    { descent = true; return MB_SUCCESS; }
  
  ErrorCode leaf( EntityHandle node ) {
    OrientedBox box;
    ErrorCode rval = mTool->box( node, box );
    EntityHandle h;
    rval = box.make_hex( h, mOut );
    if (MB_SUCCESS !=rval) return rval;
    int i = hash_handle( node );
    return mOut->tag_set_data( mTag, &h, 1, &i );
  }
};

void write_tree_blocks( Interface* interface, const char* file )
{
  EntityHandle root;
  ErrorCode rval = get_root( interface, root );
  if (MB_SUCCESS != rval) {
    std::cerr << "Internal error: Failed to retreive tree." << std::endl;
    exit(5);
  }
  
  Core moab2;
  Tag tag;
  int zero = 0;
  moab2.tag_get_handle( TAG_NAME, 1, MB_TYPE_INTEGER, tag, MB_TAG_DENSE|MB_TAG_CREAT, &zero );

  OrientedBoxTreeTool tool(interface);
  LeafHexer op( &tool, &moab2, tag );
  rval = tool.preorder_traverse( root, op );
  if (MB_SUCCESS != rval) {
    std::cerr << "Internal error: failed to construct leaf hexes" << std::endl;
    exit(5);
  }

  moab2.write_mesh( file );
}
 
