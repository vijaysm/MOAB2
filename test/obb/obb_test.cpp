#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBOrientedBoxTreeTool.hpp"
#include "MBOrientedBox.hpp"
#include "MBTagConventions.hpp"

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <limits>
#include <set>
#include <map>

#include "testdir.h"
const char* NAME = "obb_test";
const char* DEFAULT_FILES[] = { "../../" TEST_DIR "/3k-tri-sphere.vtk",
                                //"../../" TEST_DIR "/4k-tri-plane.vtk",
#ifdef HDF5_FILE
                                "../../" TEST_DIR "/3k-tri-cube.h5m",
#endif
                                0 };

static void usage( const char* error, const char* opt )
{
  const char* default_message = "Invalid option";
  if (opt && !error)
    error = default_message;

  std::ostream& str = error ? std::cerr : std::cout;
  if (error) {
    str << error;
    if (opt)
      str << ": " << opt;
    str << std::endl;
  }

  str << "Usage: "<<NAME<<" [output opts.] [settings] [file ...]" << std::endl;
  str << "       "<<NAME<<" -h" << std::endl;
  if (!error)
    str << " If no file is specified the defautl test files will be used" << std::endl
        << " -h  print help text. " << std::endl
        << " -v  verbose output (may be specified multiple times) " << std::endl
        << " -q  quite (minimal output) " << std::endl
        << " -f <x>:<y>:<z>:<i>:<j>:<k> Do ray fire" << std::endl
        << " -c  write box geometry to Cubit journal file." << std::endl
        << " -k  write leaf contents to vtk files." << std::endl
        << " -K  write contents of leaf boxes intersected by rays to vtk file." << std::endl
        << " -o <name> Save file containing tree and triangles.  Mesh tag \"OBB_ROOT\"." << std::endl
        << " -t <real> specify tolerance" << std::endl
        << " -n <int>  specify max entities per leaf node " << std::endl
        << " -l <int>  specify max tree levels" << std::endl
        << " -r <real> specify worst cell split ratio" << std::endl
        << " -R <real> specify best cell split ratio" << std::endl
        << " -s force construction of surface tree" << std::endl
        << " -S do not build surface tree." << std::endl
        << "    (Default: surface tree if file contains multiple surfaces" << std::endl
        << " -u  use unorderd (MBRange) meshsets for tree nodes" << std::endl
        << " -U  use ordered (vector) meshsets for tree nodes" << std::endl
#if MB_OOB_SPLIT_BY_NON_INTERSECTING
        << " -I <real> specify weight for intersection ratio" << std::endl
#endif
        << " Verbosity (-q sets to 0, each -v increments, default is 1):" << std::endl
        << "  0 - no output" << std::endl
        << "  1 - status messages and error summary" << std::endl
        << "  2 - print tree statistics " << std::endl
        << "  3 - output errors for each node" << std::endl
        << "  4 - print tree" << std::endl
        << "  5 - print tree w/ contents of each node" << std::endl
        << " See documentation for MMOrientedBoxTreeTool::Settings for " << std::endl
        << " a description of tree generation settings." << std::endl
      ;
  exit( !!error );
}

static const char* get_option( int& i, int argc, char* argv[] ) {
  ++i;
  if (i == argc)
    usage( "Expected argument following option", argv[i-1] );
  return argv[i];
}

static int get_int_option( int& i, int argc, char* argv[] ) {
  const char* str = get_option( i, argc, argv );
  char* end_ptr;
  long val = strtol( str, &end_ptr, 0 );
  if (!*str || *end_ptr) 
    usage( "Expected integer following option", argv[i-1] );
  return val;
}

static double get_double_option( int& i, int argc, char* argv[] ) {
  const char* str = get_option( i, argc, argv );
  char* end_ptr;
  double val = strtod( str, &end_ptr );
  if (!*str || *end_ptr) 
    usage( "Expected real number following option", argv[i-1] );
  return val;
}

static bool do_file( const char* filename );
  

enum TriOption { DISABLE, ENABLE, AUTO };

int verbosity = 1;
MBOrientedBoxTreeTool::Settings settings;
double tolerance = 1e-6;
bool write_cubit = false;
bool write_vtk = false;
bool write_ray_vtk = false;
std::vector<MBCartVect> rays;
const char* save_file_name = 0;
TriOption surfTree = AUTO;

void parse_ray( int& i, int argc, char* argv[] );

int main( int argc, char* argv[] )
{
  std::vector<const char*> file_names;
  bool flags = true;
  for (int i = 1; i < argc; ++i) {
    if (flags && argv[i][0] =='-') {
      if (!argv[i][1] || argv[i][2])
        usage(0,argv[i]);
      switch (argv[i][1]) {
        default:  usage( 0, argv[i] );  break;
        case '-': flags = false;        break;
        case 'v': ++verbosity;          break;
        case 'q': verbosity = 0;        break;
        case 'h': usage( 0, 0 );        break;
        case 'c': write_cubit = true;   break;
        case 'k': write_vtk = true;     break;
        case 'K': write_ray_vtk = true; break;
        case 's': surfTree = ENABLE;    break;
        case 'S': surfTree = DISABLE;   break;
        case 'u': settings.set_options = MESHSET_SET; break;
        case 'U': settings.set_options = MESHSET_ORDERED; break;
        case 'o':
          save_file_name = get_option( i, argc, argv );
          DEFAULT_FILES[1] = 0; // only one file can be saved, so by default only do one.
          break;
        case 'n': 
          settings.max_leaf_entities = get_int_option( i, argc, argv );
          break;
        case 'l':
          settings.max_depth = get_int_option( i, argc, argv );
          break;
        case 'r':
          settings.worst_split_ratio = get_double_option( i, argc, argv );
          break;
        case 'R':
          settings.best_split_ratio = get_double_option( i, argc, argv );
          break;
#if MB_OOB_SPLIT_BY_NON_INTERSECTING
        case 'I':
          settings.intersect_ratio_factor = get_double_option( i, argc, argv );
          break;
#endif
        case 't':
          tolerance = get_double_option( i, argc, argv );
          break;
        case 'f':
          parse_ray( i, argc, argv );
          break;
      }
    }
    else {
      file_names.push_back( argv[i] );
    }
  }
  
  if (verbosity) {
    MBCore core;
    std::string version;
    core.impl_version( & version );
    std::cout << version << std::endl;
    if (verbosity > 1) 
      std::cout << "max_leaf_entities:      " << settings.max_leaf_entities      << std::endl
                << "max_depth:              " << settings.max_depth              << std::endl
                << "worst_split_ratio:      " << settings.worst_split_ratio      << std::endl
                << "best_split_ratio:       " << settings.best_split_ratio       << std::endl
#if MB_OOB_SPLIT_BY_NON_INTERSECTING
                << "intersect ratio factor: " << settings.intersect_ratio_factor << std::endl
#endif
                << "tolerance:              " << tolerance                       << std::endl
                << "set type:               " << ((settings.set_options&MESHSET_ORDERED) ? "ordered" : "set") << std::endl
                << std::endl;
  }
  
  if (!settings.valid() || tolerance < 0.0) {
    std::cerr << "Invalid settings specified." << std::endl;
    return 2;
  }
  
  if (file_names.empty()) {
    std::cerr << "No file(s) specified." << std::endl;
    for (int i = 0; DEFAULT_FILES[i]; ++i) {
      std::cerr << "Using default file \"" << DEFAULT_FILES[i] << '"' << std::endl;
      file_names.push_back( DEFAULT_FILES[i] );
    }
  }
  
  if (save_file_name && file_names.size() != 1) {
    std::cerr << "Only one input file allowed if \"-o\" option is specified." << std::endl;
    return 1;
  }
  
  int exit_val = 0;
  for (unsigned j = 0; j < file_names.size(); ++j)
    if (!do_file( file_names[j] ))
      ++exit_val;
  
  return exit_val ? exit_val + 2 : 0;
}


void parse_ray( int& i, int argc, char* argv[] )
{
  MBCartVect point, direction;
  if (6 != sscanf( get_option( i, argc, argv ), "%lf:%lf:%lf:%lf:%lf:%lf",
                   &point[0], &point[1], &point[2],
                   &direction[0], &direction[1], &direction[2] ))
    usage( "Expected ray specified as <x>:<y>:<z>:<i>:<j>:<k>", 0 );
  direction.normalize();
  rays.push_back( point );
  rays.push_back( direction );
}
 

class TreeValidator : public MBOrientedBoxTreeTool::Op
{
  private:
    MBInterface* const instance;
    MBOrientedBoxTreeTool* const tool;
    const bool printing;
    const double epsilon;
    bool surfaces;
    std::ostream& stream;
    MBOrientedBoxTreeTool::Settings settings;
    
    void print( MBEntityHandle handle, const char* string ) {
      if (printing)
        stream << instance->id_from_handle(handle) << ": "
               << string << std::endl;
    }
    
    MBErrorCode error( MBEntityHandle handle, const char* string ) {
      ++error_count;
      print( handle, string );
      return MB_SUCCESS;
    }
    
  public:
    
    unsigned entity_count;
  
    unsigned loose_box_count;
    unsigned child_outside_count;
    unsigned entity_outside_count;
    unsigned num_entities_outside;
    unsigned non_ortho_count;
    unsigned error_count;
    unsigned empty_leaf_count;
    unsigned non_empty_non_leaf_count;
    unsigned entity_invalid_count;
    unsigned unsorted_axis_count;
    unsigned non_unit_count;
    unsigned duplicate_entity_count;
    unsigned bad_outer_radius_count;
    unsigned missing_surface_count;
    unsigned multiple_surface_count;
    std::set<MBEntityHandle> seen;
    int surface_depth;
    MBEntityHandle surface_handle;

    TreeValidator( MBInterface* instance_ptr, 
                   MBOrientedBoxTreeTool* tool_ptr,
                   bool print_errors,
                   std::ostream& str,
                   double tol,
                   bool surfs,
                   MBOrientedBoxTreeTool::Settings s )
      : instance(instance_ptr),
        tool(tool_ptr),
        printing(print_errors), 
        epsilon(tol), 
        surfaces(surfs),
        stream(str),
        settings(s),
        entity_count(0),
        loose_box_count(0),
        child_outside_count(0),
        entity_outside_count(0),
        num_entities_outside(0),
        non_ortho_count(0),
        error_count(0),
        empty_leaf_count(0),
        non_empty_non_leaf_count(0),
        entity_invalid_count(0),
        unsorted_axis_count(0),
        non_unit_count(0),
        duplicate_entity_count(0),
        bad_outer_radius_count(0),
        missing_surface_count(0),
        multiple_surface_count(0),
        surface_depth(-1)
      {}
    
    bool is_valid() const 
      { return 0 == loose_box_count+child_outside_count+entity_outside_count+
                    num_entities_outside+non_ortho_count+error_count+
                    empty_leaf_count+non_empty_non_leaf_count+entity_invalid_count
                    +unsorted_axis_count+non_unit_count+duplicate_entity_count
                    +bad_outer_radius_count+missing_surface_count+multiple_surface_count; }
    
    virtual MBErrorCode visit( MBEntityHandle node,
                               int depth,
                               bool& descend );
                               
    virtual MBErrorCode leaf( MBEntityHandle ) { return MB_SUCCESS; }
};


MBErrorCode TreeValidator::visit( MBEntityHandle node,
                                  int depth,
                                  bool& descend )
{
  MBErrorCode rval;
  descend = true;
  
  MBRange contents;
  rval = instance->get_entities_by_handle( node, contents );
  if (MB_SUCCESS != rval) 
    return error(node, "Error getting contents of tree node.  Corrupt tree?");
  entity_count += contents.size();
  
  if (surfaces) {
      // if no longer in subtree for previous surface, clear
    if (depth <= surface_depth) 
      surface_depth = -1;
      
    MBEntityHandle surface = 0;
    MBRange::iterator surf_iter = contents.lower_bound( MBENTITYSET );
    if (surf_iter != contents.end()) {
      surface = *surf_iter;
      contents.erase( surf_iter );
    }
    
    if (surface) {
      if (surface_depth >=0) {
        ++multiple_surface_count;
        print( node, "Multiple surfaces in encountered in same subtree." );
      }
      else {
        surface_depth = depth;
        surface_handle = surface;
      }
    }
  }
  
  std::vector<MBEntityHandle> children;
  rval = tool->get_moab_instance()->get_child_meshsets( node, children );
  if (MB_SUCCESS != rval || (!children.empty() && children.size() != 2)) 
    return error(node, "Error getting children.  Corrupt tree?");
  
  MBOrientedBox box;
  rval = tool->box( node, box );
  if (MB_SUCCESS != rval) 
    return error(node, "Error getting oriented box from tree node.  Corrupt tree?");
  
  if (children.empty() && contents.empty()) {
    ++empty_leaf_count;
    print( node, "Empty leaf node.\n" );
  }
  else if (!children.empty() && !contents.empty()) {
    ++non_empty_non_leaf_count;
    print( node, "Non-leaf node is not empty." );
  }
  
  if (surfaces && children.empty() && surface_depth < 0) {
    ++missing_surface_count;
    print( node, "Reached leaf node w/out encountering any surface set.");
  }
  
  double dot_epsilon = epsilon*(box.axis[0]+box.axis[1]+box.axis[2]).length();
  if (box.axis[0] % box.axis[1] > dot_epsilon ||
      box.axis[0] % box.axis[2] > dot_epsilon ||
      box.axis[1] % box.axis[2] > dot_epsilon ) {
    ++non_ortho_count;
    print (node, "Box axes are not orthogonal");
  }
  
  if (!children.empty()) {
    for (int i = 0; i < 2; ++i) {
      MBOrientedBox other_box;
      rval = tool->box( children[i], other_box );
      if (MB_SUCCESS != rval) 
        return error( children[i], " Error getting oriented box from tree node.  Corrupt tree?" );
//      else if (!box.contained( other_box, epsilon )) {
//        ++child_outside_count;
//        print( children[i], "Parent box does not contain child box." );
//        char string[64];
//        sprintf(string, "     Volume ratio is %f", other_box.volume()/box.volume() );
//        print( children [i], string );
//      }
        else {
          double vol_ratio = other_box.volume()/box.volume();
          if (vol_ratio > 2.0) {
            char string[64];
            sprintf(string, "child/parent volume ratio is %f", vol_ratio );
            print( children[i], string );
            sprintf(string, "   child/parent area ratio is %f", other_box.area()/box.area() );
            print( children[i], string );
          }
       }
    }
  }
  
  bool bad_element = false;
  bool bad_element_handle = false;
  bool bad_element_conn = false;
  bool duplicate_element = false;
  int num_outside = 0;
  bool boundary[6] = { false, false, false, false, false, false };
  for (MBRange::iterator it = contents.begin(); it != contents.end(); ++it) {
    MBEntityType type = instance->type_from_handle( *it );
    int dim = MBCN::Dimension( type );
    if (dim != 2) {
      bad_element = true;
      continue;
    }
    
    const MBEntityHandle* conn;
    int conn_len;
    rval = instance->get_connectivity( *it, conn, conn_len );
    if (MB_SUCCESS != rval) {
      bad_element_handle = true;
      continue;
    }
    
    std::vector<MBCartVect> coords(conn_len);
    rval = instance->get_coords( conn, conn_len, coords[0].array() );
    if (MB_SUCCESS != rval) {
      bad_element_conn = true;
      continue;
    }
    
    bool outside = false;
    for (std::vector<MBCartVect>::iterator j = coords.begin(); j != coords.end(); ++j) {
      if (!box.contained( *j, epsilon ))
        outside = true;
      else for (int d = 0; d < 3; ++d) {
#if MB_ORIENTED_BOX_UNIT_VECTORS
        double n = box.axis[d] % (*j - box.center);
        if (fabs(n - box.length[d]) <= epsilon)
          boundary[2*d] = true;
        if (fabs(n + box.length[d]) <= epsilon)
          boundary[2*d+1] = true;
#else
        double ln = box.axis[d].length();
        MBCartVect v1 = *j - box.center - box.axis[d];
        MBCartVect v2 = *j - box.center + box.axis[d];
        if (fabs(v1 % box.axis[d]) <= ln * epsilon)
          boundary[2*d] = true;
        if (fabs(v2 % box.axis[d]) <= ln * epsilon)
          boundary[2*d+1] = true;
#endif
      }
    }
    if (outside)
      ++num_outside;
      
    if (!seen.insert(*it).second) {
      duplicate_element = true;
      ++duplicate_entity_count;
    }
  }
  
  MBCartVect alength( box.axis[0].length(), box.axis[1].length(), box.axis[2].length() );
#if MB_ORIENTED_BOX_UNIT_VECTORS
  MBCartVect length = box.length;
#else
  MBCartVect length = alength;
#endif
  
  if (length[0] > length[1] || length[0] > length[2] || length[1] > length[2]) {
    ++unsorted_axis_count;
    print( node, "Box axes are not ordered from shortest to longest." );
  }
  
#if MB_ORIENTED_BOX_UNIT_VECTORS
  if (fabs(alength[0] - 1.0) > epsilon ||
      fabs(alength[1] - 1.0) > epsilon ||
      fabs(alength[2] - 1.0) > epsilon) {
    ++non_unit_count;
    print( node, "Box axes are not unit vectors.");
  }
#endif

#if MB_ORIENTED_BOX_OUTER_RADIUS
  if (fabs(length.length() - box.radius) > tolerance) {
    ++bad_outer_radius_count;
    print( node, "Box has incorrect outer radius.");
  }
#endif

  if (depth+1 < settings.max_depth 
      && contents.size() > (unsigned)(4*settings.max_leaf_entities))
  {
    char string[64];
    sprintf(string, "leaf at depth %d with %u entities", depth, (unsigned)contents.size() );
    print( node, string );
  }
    
      
  bool all_boundaries = true;
  for (int f = 0; f < 6; ++f)
    all_boundaries = all_boundaries && boundary[f];
  
  if (bad_element) {
    ++entity_invalid_count;
    print( node, "Set contained an entity with an inappropriate dimension." );
  }
  if (bad_element_handle) {
    ++error_count;
    print( node, "Error querying face contained in set.");
  }
  if (bad_element_conn) {
    ++error_count;
    print( node, "Error querying connectivity of element.");
  }
  if (duplicate_element) {
    print( node, "Elements occur in multiple leaves of tree.");
  }
  if (num_outside > 0) {
    ++entity_outside_count;
    num_entities_outside += num_outside;
    if (printing)
      stream << instance->id_from_handle( node ) << ": "
             << num_outside << " elements outside box." << std::endl;
  }
  else if (!all_boundaries && !contents.empty()) {
    ++loose_box_count;
    print( node, "Box does not fit contained elements tightly." );
  }

  return MB_SUCCESS;
}

class CubitWriter : public MBOrientedBoxTreeTool::Op
{
  public:
    CubitWriter( FILE* file_ptr, 
                 MBOrientedBoxTreeTool* tool_ptr )
      : file(file_ptr), tool(tool_ptr) {}
    
    MBErrorCode visit ( MBEntityHandle node,
                        int depth,
                        bool& descend );
    MBErrorCode leaf( MBEntityHandle ) { return MB_SUCCESS; }
    
  private:
    FILE* file;
    MBOrientedBoxTreeTool* tool;
};

MBErrorCode CubitWriter::visit( MBEntityHandle node,
                                int ,
                                bool& descend )
{
  descend = true;
  MBOrientedBox box;
  MBErrorCode rval = tool->box( node, box );
  if (rval != MB_SUCCESS)
    return rval;

  double sign[] = {-1, 1};
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k) {
#if MB_ORIENTED_BOX_UNIT_VECTORS
        MBCartVect corner = box.center + box.length[0] * sign[i] * box.axis[0] +
                                         box.length[1] * sign[j] * box.axis[1] +
                                         box.length[2] * sign[k] * box.axis[2];
#else
        MBCartVect corner = box.center + sign[i] * box.axis[0] +
                                         sign[j] * box.axis[1] +
                                         sign[k] * box.axis[2];
#endif
        fprintf( file, "create vertex %f %f %f\n", corner[0],corner[1],corner[2] );
      }
  fprintf( file, "#{i=Id(\"vertex\")-7}\n" );
  fprintf( file, "create surface vertex {i  } {i+1} {i+3} {i+2}\n");
  fprintf( file, "create surface vertex {i+4} {i+5} {i+7} {i+6}\n");
  fprintf( file, "create surface vertex {i+1} {i+0} {i+4} {i+5}\n");
  fprintf( file, "create surface vertex {i+3} {i+2} {i+6} {i+7}\n");
  fprintf( file, "create surface vertex {i+2} {i+0} {i+4} {i+6}\n");
  fprintf( file, "create surface vertex {i+3} {i+1} {i+5} {i+7}\n");
  fprintf( file, "delete vertex {i}-{i+7}\n");
  fprintf( file, "#{s=Id(\"surface\")-5}\n" );
  fprintf( file, "create volume surface {s} {s+1} {s+2} {s+3} {s+4} {s+5} noheal\n" );
  int id = tool->get_moab_instance()->id_from_handle( node );
  fprintf( file, "volume {Id(\"volume\")} name \"cell%d\"\n", id );
  
  return MB_SUCCESS;
}

class VtkWriter : public MBOrientedBoxTreeTool::Op
{
   public:
    VtkWriter( std::string base_name, 
               MBInterface* interface )
      : baseName(base_name), instance(interface) {}
    
    MBErrorCode visit ( MBEntityHandle node,
                        int depth,
                        bool& descend );
                        
    MBErrorCode leaf( MBEntityHandle node ) { return MB_SUCCESS; }
    
  private:
    std::string baseName;
    MBInterface* instance;
};

MBErrorCode VtkWriter::visit( MBEntityHandle node,
                              int ,
                              bool& descend )
{
  descend = true;
  
  MBErrorCode rval;
  int count;
  rval = instance->get_number_entities_by_handle( node, count );
  if (MB_SUCCESS != rval || 0 == count)
    return rval;
  
  int id = instance->id_from_handle( node );
  char id_str[32];
  sprintf( id_str, "%d", id );
  
  std::string file_name( baseName );
  file_name += ".";
  file_name += id_str;
  file_name += ".vtk";
  
  return instance->write_mesh( file_name.c_str(), &node, 1 );
}

  
static bool do_ray_fire_test( MBOrientedBoxTreeTool& tool, 
                              MBEntityHandle root_set,
                              const char* filename,
                              bool have_surface_tree );

static MBErrorCode save_tree( MBInterface* instance,
                              const char* filename,
                              MBEntityHandle tree_root );
  
static bool do_file( const char* filename )
{
  MBErrorCode rval;
  MBCore instance;
  MBInterface* const iface = &instance;
  MBOrientedBoxTreeTool tool( iface );
  bool haveSurfTree = false;
  
  if (verbosity) 
    std::cout << filename << std::endl
              << "------" << std::endl;
  
  rval = iface->load_mesh( filename );
  if (MB_SUCCESS != rval) {
    if (verbosity)
      std::cout << "Failed to read file: \"" << filename << '"' << std::endl;
     return false;
  }
  
    // IF building from surfaces, get surfaces.
    // If AUTO and less than two surfaces, don't build from surfaces.
  MBRange surfaces;
  if (surfTree != DISABLE) {
    MBTag surftag;
    rval = iface->tag_get_handle( GEOM_DIMENSION_TAG_NAME, surftag );
    if (MB_SUCCESS == rval) {
      int dim = 2;
      const void* tagvalues[] = {&dim};
      rval = iface->get_entities_by_type_and_tag( 0, MBENTITYSET,
                                &surftag, tagvalues, 1, surfaces );
      if (MB_SUCCESS != rval && MB_ENTITY_NOT_FOUND != rval)
        return false;
    }
    else if (MB_TAG_NOT_FOUND != rval) 
      return false;
    
    if (ENABLE == surfTree && surfaces.empty()) {
      std::cerr << "No Surfaces found." << std::endl;
      return false;
    }
    
    haveSurfTree = (ENABLE == surfTree) || (surfaces.size() > 1);
  }
  
  MBEntityHandle root;
  MBRange entities;
  if (!haveSurfTree) {
    rval = iface->get_entities_by_dimension( 0, 2, entities );
    if (MB_SUCCESS != rval) {
      std::cerr << "get_entities_by_dimension( 2 ) failed." << std::endl;
      return false;
    }
  
    if (entities.empty()) {
      if (verbosity)
        std::cout << "File \"" << filename << "\" contains no 2D elements" << std::endl;
      return false;
    }
  
    if (verbosity) 
      std::cout << "Building tree from " << entities.size() << " 2D elements" << std::endl;

    rval = tool.build( entities, root, &settings );
    if (MB_SUCCESS != rval) {
      if (verbosity)
        std::cout << "Failed to build tree." << std::endl;
      return false;
    }
  }
  else {
      // Get list of entities so can verify total count later
    for (MBRange::iterator s = surfaces.begin(); s != surfaces.end(); ++s) {
      rval= iface->get_entities_by_dimension( *s, 2, entities );
      if (MB_SUCCESS != rval)
        return false;
    }

    if (verbosity)
      std::cout << "Building tree from " << surfaces.size() << " surfaces" 
                << " (" << entities.size() << " elements)" << std::endl;
    entities.merge( surfaces );
    
    rval = tool.set_build( surfaces, root, &settings );
    if (MB_SUCCESS != rval) {
      if (verbosity)
        std::cout << "Failed to build tree." << std::endl;
      return false;
    }
   
      // Get list of entities so can verify total count later
    for (MBRange::iterator s = surfaces.begin(); s != surfaces.end(); ++s) {
      rval= iface->get_entities_by_dimension( *s, 2, entities );
      if (MB_SUCCESS != rval)
        return false;
    }
  }    

  if (write_cubit) {
    std::string name = filename;
    name += ".boxes.jou";
    FILE* ptr = fopen( name.c_str(), "w+" );
    if (!ptr) {
      perror( name.c_str() );
    }
    else {
      if (verbosity)
        std::cout << "Writing: " << name << std::endl;
      fprintf(ptr,"graphics off\n");
      CubitWriter op( ptr, &tool );
      tool.preorder_traverse( root, op );
      fprintf(ptr,"graphics on\n");
      fclose( ptr );
    }
  }
  
  if (write_vtk) {
    VtkWriter op( filename, iface );
    if (verbosity)
      std::cout << "Writing leaf contents as : " << filename 
                << ".xxx.vtk where 'xxx' is the set id" << std::endl;
    tool.preorder_traverse( root, op );
  }  

  bool print_errors = false, print_contents = false;
  switch (verbosity) {
    default:
      print_contents = true;
    case 4:
      tool.print( root, std::cout, print_contents );
    case 3:
      print_errors = true;
    case 2:
      rval = tool.stats( root, std::cout );
      if (MB_SUCCESS != rval)
        std::cout << "****** Failed to get tree statistics ******" << std::endl;
    case 1:
    case 0:
      ;
  }  
  
  TreeValidator op( iface, &tool, print_errors, std::cout, tolerance, haveSurfTree, settings ); 
  rval = tool.preorder_traverse( root, op );
  bool result = op.is_valid();
  if (MB_SUCCESS != rval) {
    result = false;
    if (verbosity)
      std::cout << "Errors traversing tree.  Corrupt tree?" << std::endl;
  }
  
  bool missing = (op.entity_count != entities.size());
  if (missing)
    result = false;
  
  if (verbosity) {
    if (result)
      std::cout << std::endl << "No errors detected." << std::endl;
    else
      std::cout << std::endl << "*********************** ERROR SUMMARY **********************" << std::endl;
    if (op.child_outside_count)
      std::cout << "* " << op.child_outside_count << " child boxes not contained in parent." << std::endl;
    if (op.entity_outside_count)
      std::cout << "* " << op.entity_outside_count << " nodes containing entities outside of box." << std::endl;
    if (op.num_entities_outside)
      std::cout << "* " << op.num_entities_outside << " entities outside boxes." << std::endl;
    if (op.empty_leaf_count)
      std::cout << "* " << op.empty_leaf_count << " empty leaf nodes." << std::endl;
    if (op.non_empty_non_leaf_count)
      std::cout << "* " << op.non_empty_non_leaf_count << " non-leaf nodes containing entities." << std::endl;
    if (op.duplicate_entity_count)
      std::cout << "* " << op.duplicate_entity_count << " duplicate entities in leaves." << std::endl;
    if (op.missing_surface_count)
      std::cout << "* " << op.missing_surface_count << " leaves outside surface subtrees." << std::endl;
    if (op.multiple_surface_count)
      std::cout << "* " << op.multiple_surface_count << " surfaces within other surface subtrees." << std::endl;
    if (op.non_ortho_count)
      std::cout << "* " << op.non_ortho_count << " boxes with non-orthononal axes." << std::endl;
    if (op.non_unit_count)
      std::cout << "* " << op.non_unit_count << " boxes with non-unit axes." << std::endl;
    if (op.bad_outer_radius_count)
      std::cout << "* " << op.bad_outer_radius_count << " boxes incorrect outer radius." << std::endl;
    if (op.unsorted_axis_count)
      std::cout << "* " << op.unsorted_axis_count << " boxes axes in unsorted order." << std::endl;
    if (op.loose_box_count)
      std::cout << "* " << op.loose_box_count << " boxes that do not fit the contained entities tightly." << std::endl;
    if (op.error_count + op.entity_invalid_count)
      std::cout << "* " << op.error_count + op.entity_invalid_count
                << " other errors while traversing tree." << std::endl;
    if (missing)
      std::cout << "* tree built from " << entities.size() << " entities contains "
                << op.entity_count << " entities." << std::endl;
    if (!result)
      std::cout << "************************************************************" << std::endl;
  }
  
  if (result && save_file_name) {
    if (MB_SUCCESS == save_tree( iface, save_file_name, root ))
      std::cerr << "Wrote '" << save_file_name << "'" << std::endl;
    else
      std::cerr << "FAILED TO WRITE '" << save_file_name << "'" << std::endl;
  }
  
  if (!do_ray_fire_test( tool, root, filename, haveSurfTree )) {
    if (verbosity)
      std::cout << "Ray fire test failed." << std::endl;
    result = false;
  }

  rval = tool.delete_tree( root );
  if (MB_SUCCESS != rval) {
    if (verbosity)
      std::cout << "delete_tree failed." << std::endl;
    result = false;
  }
  
  return result;
}

struct RayTest {
  const char* description;
  unsigned expected_hits;
  MBCartVect point, direction;
};

static bool do_ray_fire_test( MBOrientedBoxTreeTool& tool, 
                              MBEntityHandle root_set,
                              const char* filename,
                              bool haveSurfTree )
{
  if (verbosity > 1)
    std::cout << "beginning ray fire tests" << std::endl;
 
  MBOrientedBox box;
  MBErrorCode rval = tool.box( root_set, box );
  if (MB_SUCCESS != rval) {
    if (verbosity)
      std::cerr << "Error getting box for tree root set" << std::endl;
    return false;
  }
  
  /* Do standard ray fire tests */
  
  RayTest tests[] = { 
   { "half-diagonal from center", 1, box.center,                            1.5 * box.dimensions() },
   { "large axis through box",    2, box.center - 1.2 * box.scaled_axis(2), box.axis[2] },
   { "small axis through box",    2, box.center - 1.2 * box.scaled_axis(0), box.axis[0] },
   { "parallel miss",             0, box.center + 2.0 * box.scaled_axis(1), box.axis[2] },
   { "skew miss",                 0, box.center + box.dimensions(),          box.dimensions() * box.axis[2] }
   };
  
  bool result = true;
  const size_t num_test = sizeof(tests)/sizeof(tests[0]);
  for (size_t i = 0; i < num_test; ++i) {
    tests[i].direction.normalize();
    if (verbosity > 2) 
      std::cout << "  " << tests[i].description << " " << tests[i].point << " " << tests[i].direction << std::endl;
    
    std::vector<double> intersections;
    rval = tool.ray_intersect_triangles( intersections, root_set, tolerance, tests[i].point.array(), tests[i].direction.array(), 0 );
    if (MB_SUCCESS != rval) {
      if (verbosity)
        std::cout << "  Call to MBOrientedBoxTreeTool::ray_intersect_triangles failed." << std::endl;
      result = false;
      continue;
    }
    
    if (intersections.size() != tests[i].expected_hits) {
      if (verbosity > 2)
        std::cout << "  Expected " << tests[i].expected_hits << " and got "
                  << intersections.size() << " hits for ray fire of " 
                  << tests[i].description << std::endl;
      if (verbosity > 3) {
        for (unsigned j = 0; j < intersections.size(); ++j)
          std::cout << "  " << intersections[j];
        std::cout << std::endl;
      }
      result = false;
    }
    
    if (!haveSurfTree)
      continue;
    
    const int NUM_NON_TOL_INT = 1;
    std::vector<double> intersections2;
    std::vector<MBEntityHandle> surf_handles;
    rval = tool.ray_intersect_sets( intersections2, surf_handles, root_set, tolerance, NUM_NON_TOL_INT, tests[i].point.array(), tests[i].direction.array(), 0 );
    if (MB_SUCCESS != rval) {
      if (verbosity)
        std::cout << "  Call to MBOrientedBoxTreeTool::ray_intersect_sets failed." << std::endl;
      result = false;
      continue;
    }

    if (surf_handles.size() != intersections2.size()) {
      if (verbosity)
        std::cout << "  ray_intersect_sets did not return sets for all intersections." << std::endl;
      result = false;
    }

    double non_tol_dist2 = std::numeric_limits<double>::max();
    int non_tol_count2 = 0;
    for (size_t i = 0; i < intersections2.size(); ++i) {
      if (intersections2[i] > tolerance) {
        ++non_tol_count2;
        if (intersections2[i] < non_tol_dist2)
          non_tol_dist2 = intersections2[i];
      }
    }

    if (non_tol_count2 > NUM_NON_TOL_INT) {
      if (verbosity)
        std::cout << "  Requested " << NUM_NON_TOL_INT << "intersections "
                  << "  beyond tolernace.  Got " << non_tol_count2 << std::endl;
      result = false;
    }

    double non_tol_dist = std::numeric_limits<double>::max();
    int non_tol_count = 0;
    for (size_t i = 0; i < intersections.size(); ++i) {
      if (intersections[i] > tolerance) {
        ++non_tol_count;
        if (intersections[i] < non_tol_dist)
          non_tol_dist = intersections[i];
      }
    }

    if (!NUM_NON_TOL_INT)
      continue;

    if (!non_tol_count && non_tol_count2) {
      if (verbosity)
        std::cout << "  ray_intersect_sets returned intersection not found by ray_intersect_triangles" << std::endl;
      result = false;
      continue;
    }
    else if (non_tol_count && !non_tol_count2) {
      if (verbosity)
        std::cout << "  ray_intersect_sets missed intersection found by ray_intersect_triangles" << std::endl;
      result = false;
      continue;
    }
    else if (non_tol_count && non_tol_count2 && 
             fabs(non_tol_dist - non_tol_dist2) > tolerance) {
      if (verbosity)
        std::cout << "  ray_intersect_sets and ray_intersect_triangles did not find same closest intersection" << std::endl;
      result = false;
    }
  }
  
  /* Do ray fire for any user-specified rays */
  
  for (size_t i = 0; i < rays.size(); i += 2) {
    std::cout << rays[i] << "+" << rays[i+1] << " : ";
    
    if (!haveSurfTree) {
      MBRange leaves;
      std::vector<double> intersections;
      rval = tool.ray_intersect_boxes( leaves, root_set, tolerance, rays[i].array(), rays[i+1].array(), 0 );
      if (MB_SUCCESS != rval) {
        std::cout << "FAILED" << std::endl;
        result = false;
        continue;
      }

      if (!leaves.empty() && write_ray_vtk) {
        std::string num, name(filename);
        std::stringstream s;
        s << (i/2);
        s >> num;
        name += "-ray";
        name += num;
        name += ".vtk";

        std::vector<MBEntityHandle> sets(leaves.size());
        std::copy( leaves.begin(), leaves.end(), sets.begin() );
        tool.get_moab_instance()->write_mesh( name.c_str(), &sets[0], sets.size() );
        if (verbosity)
          std::cout << "(Wrote " << name << ") ";
      }      

      rval = tool.ray_intersect_triangles( intersections, leaves, tolerance, rays[i].array(), rays[i+1].array(), 0 );
      if (MB_SUCCESS != rval) {
        std::cout << "FAILED" << std::endl;
        result = false;
        continue;
      }

      if (intersections.empty()) {
        std::cout << "(none)" << std::endl;
        continue;
      }

      std::cout << intersections[0];
      for (unsigned j = 1; j < intersections.size(); ++j)
        std::cout << ", " << intersections[j];
      std::cout << std::endl;

      if (!leaves.empty() && write_cubit && verbosity > 2) {
        std::cout << "  intersected boxes:";
        for (MBRange::iterator i= leaves.begin(); i!= leaves.end(); ++i)
          std::cout << " " << tool.get_moab_instance()->id_from_handle(*i);
        std::cout << std::endl;
      }
    }
    else {
      std::vector<double> intersections;
      std::vector<MBEntityHandle> surfaces;

      rval = tool.ray_intersect_sets( intersections, surfaces, root_set, tolerance, 1000, rays[i].array(), rays[i+1].array(), 0 );
      if (MB_SUCCESS != rval) {
        std::cout << "FAILED" << std::endl;
        result = false;
        continue;
      }

      if (!surfaces.empty() && write_ray_vtk) {
        std::string num, name(filename);
        std::stringstream s;
        s << (i/2);
        s >> num;
        name += "-ray";
        name += num;
        name += ".vtk";

        tool.get_moab_instance()->write_mesh( name.c_str(), &surfaces[0], surfaces.size() );
        if (verbosity)
          std::cout << "(Wrote " << name << ") ";
      }
      
      if (intersections.size() != surfaces.size()) {
        std::cout << "Mismatched output lists." << std::endl;
        result = false;
        continue;
      }
      
      if (intersections.empty()) {
        std::cout << "(none)" << std::endl;
        continue;
      }
      
      
      MBTag idtag;
      rval = tool.get_moab_instance()->tag_get_handle( GLOBAL_ID_TAG_NAME, idtag );
      if (MB_SUCCESS != rval) {
        std::cout << "NO GLOBAL_ID TAG." << std::endl;
        continue;
      }
      std::vector<int> ids(surfaces.size());
      rval = tool.get_moab_instance()->tag_get_data( idtag, &surfaces[0], surfaces.size(), &ids[0] );
      if (MB_SUCCESS != rval) {
        std::cout << "NO GLOBAL_ID TAG ON SURFACE." << std::endl;
        continue;
      }
      
        // group by surfaces
      std::map<int,double> intmap;
      for (unsigned j = 0; j < intersections.size(); ++j)
        intmap[ids[j]] = intersections[j];
        
      std::map<int,double>::iterator it = intmap.begin();
      int prevsurf = it->first;
      std::cout << "Surf" << it->first << " " << it->second;
      for (++it; it != intmap.end(); ++it) {
        std::cout << ", ";
        if (it->first != prevsurf) {
          prevsurf = it->first;
          std::cout << "Surf" << it->first << " ";
        }
        std::cout << it->second;
      }
      std::cout << std::endl;
    }
  }
  
  return result;
}


MBErrorCode save_tree( MBInterface* instance,
                       const char* filename,
                       MBEntityHandle tree_root )
{
  MBErrorCode rval;
  MBTag tag;
  
  rval = instance->tag_get_handle( "OBB_ROOT", tag );
  if (MB_SUCCESS == rval) {
    int size;
    MBDataType type;
    rval = instance->tag_get_size( tag, size );
    if (MB_SUCCESS != rval)
      return rval;
    rval = instance->tag_get_data_type( tag, type );
    if (MB_SUCCESS != rval)
      return rval;

    if (size != sizeof(MBEntityHandle) || type != MB_TYPE_HANDLE)
      return MB_FAILURE;
  }
  else {
    rval = instance->tag_create( "OBB_ROOT", sizeof(MBEntityHandle), MB_TAG_SPARSE, MB_TYPE_HANDLE, tag, 0 );
    if (MB_SUCCESS != rval)
      return rval;
  }
  
  rval = instance->tag_set_data( tag, 0, 0, &tree_root );
  if (MB_SUCCESS != rval)
    return rval;
  
  return instance->write_mesh( filename );
}


  

