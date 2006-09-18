#include "AcisQueryEngine.hpp"
#include "AcisModifyEngine.hpp"
#include "CGMApp.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryQueryEngine.hpp"
#include "ModelQueryEngine.hpp"
#include "RefEntityName.hpp"
#include "GMem.hpp"

#include "RefGroup.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"

#include "SenseEntity.hpp"
#include "Surface.hpp"
#include "Curve.hpp"

#include "MBCore.hpp"
#include "MBInterface.hpp"
#include "MBRange.hpp"
#include "MBTagConventions.hpp"

#include <stdio.h>

#include "cubfile.h"
#include "geom_file_type.h"

// Tag name used for saving sense of faces in volumes.
// We assume that the surface occurs in at most two volumes.
// Code will error out if more than two volumes per surface.
// The tag data is a pair of tag handles, representing the
// forward and reverse volumes, respectively.  If a surface
// is non-manifold in a single volume, the same volume will
// be listed for both the forward and reverse slots.
const char GEOM_SENSE_TAG_NAME[] = "GEOM_SENSE_2";

const char NAME_TAG_NAME[] = "NAME";
const int NAME_TAG_SIZE = 32;

// default settings
const char DEFAULT_NAME[] = "cgm2moab";
const double DEFAULT_DISTANCE = 0.001;
const int DEFAULT_NORMAL = 5;

// settings
const char* name = DEFAULT_NAME;
double dist_tol = DEFAULT_DISTANCE;
int norm_tol = DEFAULT_NORMAL;
double len_tol = 0.0;
bool have_dist_tol = true;
bool have_norm_tol = true;
bool have_len_tol = false;
bool actuate_attribs = true;
bool have_actuate_attribs = false;
const char* file_type = 0;

void print_usage()
{
  std::cout << name << 
      " [-d <tol>|-D] [-n <tol>|-N] [-l <tol>|-L] [-A|-a] [-t <type>] <input_file> <outupt_file>" <<
      std::endl;
  std::cout << name << " -h" << std::endl;
}

void usage()
{
  print_usage();
  exit(1);
}

void help()
{
  print_usage();
  std::cout << "\t-d  max. distance deviation between facet and geometry" << std::endl
            << "\t-D  no distance tolerance" << std::endl
            << "\t    (default:" << DEFAULT_DISTANCE << ")" <<std::endl
            << "\t-n  max. normal angle deviation (degrees) between facet and geometry" << std::endl
            << "\t-N  no normal tolerance" << std::endl
            << "\t-l  max. facet edge length" << std::endl
            << "\t-L  no facet edge length maximum" << std::endl
            << "\t    (default:none)" << std::endl
            << "\t-a  force actuation of all CGM attributes" << std::endl
            << "\t-A  disable all CGM attributes" << std::endl
            << "\t-t  specify input file type (default is autodetect)" << std::endl
            << std::endl;
  exit(0);
}


double parse_tol( char* argv[], int argc, int& i );
bool copy_geometry( MBInterface* to_this );

int main( int argc, char* argv[] )
{
  name = argv[0];
  
  const char* input_name = 0;
  const char* output_name = 0;
  
    // Process CL args
  bool process_options = true;
  for (int i = 1; i < argc; ++i) {
    if (!process_options || argv[i][0] != '-') {
      if (output_name) {
        std::cerr << "Unexpected argument: " << argv[i] << std::endl;
        usage();
      }
      else if (input_name)
        output_name = argv[i];
      else 
        input_name = argv[i];
      continue;
    }
    
    if (!argv[i][1] || argv[i][2]) {  // two chars long
      std::cerr << "Invalid option: " << argv[i] << std::endl;
      usage();
    }
    
    switch(argv[i][1]) {
      default:
        std::cerr << "Invalid option: " << argv[i] << std::endl;
        usage();
      case 'd':
        dist_tol = parse_tol( argv, argc, i );
        have_dist_tol = true;
        break;
      case 'D':
        have_dist_tol = false;
        break;
      case 'n':
        norm_tol = (int)round(parse_tol( argv, argc, i ));
        have_norm_tol = true;
        break;
      case 'N':
        have_norm_tol = false;
        break;
      case 'l':
        len_tol = parse_tol( argv, argc, i );
        have_len_tol = true;
        break;
      case 'L':
        have_len_tol = false;
        break;
      case 'a':
        actuate_attribs = true;
        have_actuate_attribs = true;
        break;
      case 'A':
        actuate_attribs = false;
        have_actuate_attribs = true;
        break;
      case 't':
        ++i;
        if (i == argc) {
          std::cerr << "Expected argument following '-t'" << std::endl;
          usage();
        }
        file_type = argv[i];
        break;
      case 'h':
        help();
      case '-':
        process_options = false;
        break;
    }
  }
  
  if (!output_name)
    usage();
        
  
    // Initialize CGM
  CGMApp::instance()->startup( 0, NULL );
  AcisQueryEngine::instance();
  AcisModifyEngine::instance();
  if (have_actuate_attribs) {
    CGMApp::instance()->attrib_manager()->set_all_auto_read_flags( actuate_attribs );
    CGMApp::instance()->attrib_manager()->set_all_auto_actuate_flags( actuate_attribs );
  }
  
  
    // Intitalize MOAB
  MBCore moab;
  MBInterface* iface = &moab;
  
    // Get CGM file type
  if (!file_type) {
    file_type = get_geom_file_type( input_name );
    if (!file_type) {
      std::cerr << input_name << " : unknown file type, try '-t'" << std::endl;
      exit(1);
    }
  }
  
    // If CUB file, extract ACIS geometry
  CubitStatus s;
  if (!strcmp( file_type, "CUBIT" )) {
    char* temp_name = tempnam( 0, "cgm" );
    if (!temp_name) {
      perror(name);
      exit(2);
    }
    
    FILE* cub_file = fopen( input_name, "r" );
    if (!cub_file) {
      free(temp_name);
      perror(input_name);
      exit(2);
    }
    
    FILE* tmp_file = fopen( temp_name, "w" );
    if (!tmp_file) {
      free(temp_name);
      perror(temp_name);
      exit(2);
    }
    
    int rval = cub_file_type( cub_file, tmp_file, CUB_FILE_ACIS );
    fclose( cub_file );
    fclose( tmp_file );
    if (rval) {
      remove( temp_name );
      free( temp_name );
      exit(2);
    }
    
    s = GeometryQueryTool::instance()->import_solid_model( temp_name, "ACIS_SAT" );
    remove( temp_name );
    free( temp_name );
  }
  else {
    s = GeometryQueryTool::instance()->import_solid_model( input_name, file_type );
  }
  if (CUBIT_SUCCESS != s) {
    std::cerr << "Failed to read '" << input_name << "' of type '" << file_type << "'" << std::endl;
    exit(2);
  }
  
    // copy geometry facets into mesh database
  if (!copy_geometry( iface )) {
    std::cerr << "Internal error copying geometry" << std::endl;
    exit(5);
  }
  
    // write mesh database
  MBErrorCode rval = iface->write_mesh( output_name );
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to write '" << output_name << "'" << std::endl;
    exit(2);
  }
  
  return 0;
}

// parse double CL arg
double parse_tol( char* argv[], int argc, int& i )
{
  ++i;
  if (i == argc) {
    std::cerr << "Expected value following option '" << argv[i-1] << "'" << std::endl;
    usage();
  }
  
  char* endptr = 0;
  double val = strtod( argv[i], &endptr );
  if (!endptr || *endptr || val < DBL_EPSILON) {
    std::cerr << "Invalid tolerance value for '" << argv[i-1] <<"': " << argv[i] << std::endl;
    usage();
  }
  
  return val;
}

// copy geometry into mesh database
bool copy_geometry( MBInterface* iface )
{
  MBErrorCode rval;

    // get some tag handles
  MBTag geom_tag, id_tag, sense_tag, name_tag;
  rval = iface->tag_create( GEOM_DIMENSION_TAG_NAME, sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER, geom_tag, 0, true );
  assert(!rval);
  rval = iface->tag_create( GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER, id_tag, 0, true );
  assert(!rval);
  rval = iface->tag_create( GEOM_SENSE_TAG_NAME, 2*sizeof(MBEntityHandle), MB_TAG_SPARSE, MB_TYPE_HANDLE, sense_tag, 0, true );
  assert(!rval);
  rval = iface->tag_create( NAME_TAG_NAME, NAME_TAG_SIZE, MB_TAG_SPARSE, MB_TYPE_OPAQUE, name_tag, 0, true );
  assert(!rval);
  
    // CGM data
  std::map<RefEntity*,MBEntityHandle> entmap[5]; // one for each dim, and one for groups
  std::map<RefEntity*,MBEntityHandle>::iterator ci;
  const char* const names[] = { "Vertex", "Curve", "Surface", "Volume" };
  DLIList<RefEntity*> entlist;
  DLIList<ModelEntity*> me_list;

    // create entity sets for all geometric entities
  for (int dim = 0; dim < 4; ++dim) {
    entlist.clean_out();
    GeometryQueryTool::instance()->ref_entity_list( names[dim], entlist, true );
    
    entlist.reset();
    for (int i = entlist.size(); i--; ) {
      RefEntity* ent = entlist.get_and_step();
      MBEntityHandle handle;
      rval = iface->create_meshset( dim == 1 ? MESHSET_ORDERED : MESHSET_SET, handle );
      if (MB_SUCCESS != rval)
        return false;
    
      entmap[dim][ent] = handle;
      
      rval = iface->tag_set_data( geom_tag, &handle, 1, &dim );
      if (MB_SUCCESS != rval)
        return rval;
      int id = ent->id();
      rval = iface->tag_set_data( id_tag, &handle, 1, &id );
      if (MB_SUCCESS != rval)
        return rval;
    }
  }
  
    // create topology for all geometric entities
  for (int dim = 1; dim < 4; ++dim) {
    for (ci = entmap[dim].begin(); ci != entmap[dim].end(); ++ci) {
      entlist.clean_out();
      ci->first->get_child_ref_entities( entlist );
    
      entlist.reset();
      for (int i = entlist.size(); i--; ) {
        RefEntity* ent = entlist.get_and_step();
        MBEntityHandle h = entmap[dim-1][ent];
        rval = iface->add_parent_child( ci->second, h );
        if (MB_SUCCESS != rval)
          return false;
      }
    }
  }
  
    // store CoFace senses
  for (ci = entmap[2].begin(); ci != entmap[2].end(); ++ci) {
    RefFace* face = (RefFace*)(ci->first);
    BasicTopologyEntity *forward = 0, *reverse = 0;
    for (SenseEntity* cf = face->get_first_sense_entity_ptr();
         cf; cf = cf->next_on_bte()) {
      BasicTopologyEntity* vol = cf->get_parent_basic_topology_entity_ptr();
      if (cf->get_sense() == CUBIT_UNKNOWN || 
          cf->get_sense() != face->get_surface_ptr()->bridge_sense()) {
        if (reverse) {
          std::cout << "Surface " << face->id() << " has reverse senes " <<
                       "with multiple volume " << reverse->id() << " and " <<
                       "volume " << vol->id() << std::endl;
          return false;
        }
        reverse = vol;
      }
      if (cf->get_sense() == CUBIT_UNKNOWN || 
          cf->get_sense() == face->get_surface_ptr()->bridge_sense()) {
        if (forward) {
          std::cout << "Surface " << face->id() << " has forward senes " <<
                       "with multiple volume " << forward->id() << " and " <<
                       "volume " << vol->id() << std::endl;
          return false;
        }
        forward = vol;
      }
    }
    
    if (forward || reverse) {
      MBEntityHandle tag_data[2] = {0,0};
      if (forward)
        tag_data[0] = entmap[3][forward];
      if (reverse)
        tag_data[1] = entmap[3][reverse];
      MBEntityHandle h = ci->second;
      rval = iface->tag_set_data( sense_tag, &h, 1, tag_data );
      if (MB_SUCCESS != rval)
        return false;
    }
  }
  
    // create entity sets for all ref groups
  DLIList<CubitString*> name_list;
  entlist.clean_out();
  GeometryQueryTool::instance()->ref_entity_list( "group", entlist );
  entlist.reset();
  for (int i = entlist.size(); i--; ) {
    RefEntity* grp = entlist.get_and_step();
    name_list.clean_out();
    RefEntityName::instance()->get_refentity_name( grp, name_list, true );
    if (name_list.size() == 0)
      continue;
    CubitString name = *name_list.get();
    
    MBEntityHandle h;
    rval = iface->create_meshset( MESHSET_SET, h );
    if (MB_SUCCESS != rval)
      return false;
    
    char namebuf[NAME_TAG_SIZE];
    memset( namebuf, '\0', NAME_TAG_SIZE );
    strncpy( namebuf, name.c_str(), NAME_TAG_SIZE - 1 );
    if (name.length() >= (unsigned)NAME_TAG_SIZE) 
      std::cout << "WARNING: group name '" << name.c_str() 
                << "' truncated to '" << namebuf << "'" << std::endl;
    rval = iface->tag_set_data( name_tag, &h, 1, namebuf );
    if (MB_SUCCESS != rval)
      return false;
      
    int id = grp->id();
    rval = iface->tag_set_data( id_tag, &h, 1, &id );
    if (MB_SUCCESS != rval)
      return false;
      
    entmap[4][grp] = h;
  }
  
    // store contents for each group
  entlist.reset();
  for (ci = entmap[4].begin(); ci != entmap[4].end(); ++ci) {
    RefGroup* grp = (RefGroup*)(ci->first);
    entlist.clean_out();
    grp->get_child_ref_entities( entlist );
    
    MBRange entities;
    while (entlist.size()) {
      RefEntity* ent = entlist.pop();
      int dim = ent->dimension();
      if (dim < 0) {
        if (entmap[4].find(ent) != entmap[4].end())
          entities.insert( entmap[4][ent] );
      }
      else if (dim < 4) {
        if (entmap[dim].find(ent) != entmap[dim].end())
          entities.insert( entmap[dim][ent] );
      }
    }
    
    if (!entities.empty()) {
      rval = iface->add_entities( ci->second, entities );
      if (MB_SUCCESS != rval)
        return false;
    }
  }
  
    // done with volumes and groups
  entmap[3].clear();
  entmap[4].clear();
  
    // create geometry for all vertices and replace 
    // vertex set handles with vertex handles in map
  for (ci = entmap[0].begin(); ci != entmap[0].end(); ++ci) {
    CubitVector pos = dynamic_cast<RefVertex*>(ci->first)->coordinates();
    double coords[3] = {pos.x(), pos.y(), pos.z()};
    MBEntityHandle vh;
    rval = iface->create_vertex( coords, vh );
    if (MB_SUCCESS != rval)
      return false;
    
    rval = iface->add_entities( ci->second, &vh, 1 );
    if (MB_SUCCESS != rval)
      return false;
    
    ci->second = vh;
  }
  
    // create geometry for all curves
  GMem data;
  for (ci = entmap[1].begin(); ci != entmap[1].end(); ++ci) {
    RefEdge* edge = dynamic_cast<RefEdge*>(ci->first);
    Curve* curve = edge->get_curve_ptr();
    GeometryQueryEngine* gqe = curve->get_geometry_query_engine();
    data.clean_out();
    int count;
    CubitStatus s = gqe->get_graphics( curve, count, &data, have_dist_tol ? dist_tol : 0.0 );
    if (CUBIT_SUCCESS != s)
      return false;
      
    std::vector<CubitVector> points;
    for (int i = 0; i < data.pointListCount; ++i)
      points.push_back( CubitVector( data.point_list()[i].x,
                                     data.point_list()[i].y,
                                     data.point_list()[i].z ) );

      // need to reverse data?
    if (curve->bridge_sense() == CUBIT_REVERSED) 
      std::reverse( points.begin(), points.end() );
    
       // check for closed curve
    RefVertex *start_vtx, *end_vtx;
    start_vtx = edge->start_vertex();
    end_vtx = edge->end_vertex();
    
      // Special case for point curve
    if (points.size() < 2) {
      if (start_vtx != end_vtx || curve->measure() > GEOMETRY_RESABS) {
        std::cerr << "Warning: No facetting for curve " << edge->id() << std::endl;
        continue;
      }
      MBEntityHandle h = entmap[0][start_vtx];
      rval = iface->add_entities( ci->second, &h, 1 );
      if (MB_SUCCESS != rval)
        return false;
      continue;
    }
    
    const bool closed = (points.front() - points.back()).length() < GEOMETRY_RESABS;
    if (closed != (start_vtx == end_vtx)) {
      std::cerr << "Warning: topology and geometry inconsistant for possibly closed curve "
                << edge->id() << std::endl;
    }
    
      // check proximity of vertices to end coordinates
    if ((start_vtx->coordinates() - points.front()).length() > GEOMETRY_RESABS
     || (  end_vtx->coordinates() - points.back() ).length() > GEOMETRY_RESABS) {
      std::cerr << "Warning: vertices not at ends of curve " << edge->id() << std::endl;
    }
    
      // create interior points
    std::vector<MBEntityHandle> verts, edges;
    verts.push_back( entmap[0][start_vtx] );
    for (size_t i = 1; i < points.size() - 1; ++i) {
      double coords[] = { points[i].x(), points[i].y(), points[i].z() };
      MBEntityHandle h;
      rval = iface->create_vertex( coords, h );
      if (MB_SUCCESS != rval)
        return false;
      verts.push_back( h );
    }
    verts.push_back( entmap[0][end_vtx] );
    
      // create edges
    for (size_t i = 0; i < verts.size()-1; ++i) {
      MBEntityHandle h;
      rval = iface->create_element( MBEDGE, &verts[i], 2, h );
      if (MB_SUCCESS != rval)
        return false;
      edges.push_back( h );
    }
    
      // if closed, remove duplicate
    if (verts.front() == verts.back())
      verts.pop_back();
    
    rval = iface->add_entities( ci->second, &verts[0], verts.size() );
    if (MB_SUCCESS != rval)
      return false;
    rval = iface->add_entities( ci->second, &edges[0], edges.size() );
    if (MB_SUCCESS != rval)
      return false;
  }
  
    // create geometry for all surfaces
  for (ci = entmap[2].begin(); ci != entmap[2].end(); ++ci) {
    RefFace* face = dynamic_cast<RefFace*>(ci->first);
    Surface* surf = face->get_surface_ptr();
    GeometryQueryEngine* gqe = surf->get_geometry_query_engine();
    data.clean_out();
    int nt, np, nf;
    CubitStatus s = gqe->get_graphics( surf, nt, np, nf, &data, 
                                       have_norm_tol ? norm_tol : 0, 
                                       have_dist_tol ? dist_tol : 0.0,
                                       have_len_tol  ?  len_tol : 0.0  );
    if (CUBIT_SUCCESS != s)
      return false;

      // declare array of all vertex handles
    std::vector<MBEntityHandle> verts( data.pointListCount, 0 );
    
      // get list of vertices in surface
    me_list.clean_out();
    ModelQueryEngine::instance()->query_model( *face, DagType::ref_vertex_type(), me_list );

      // for each vertex, find coincident point in facets
    for (int i = me_list.size(); i--; ) {
      RefVertex* vtx = dynamic_cast<RefVertex*>(me_list.get_and_step());
      CubitVector pos = vtx->coordinates();

      for (int j = 0; j < data.pointListCount; ++j) {
        CubitVector vpos( data.point_list()[j].x,
                          data.point_list()[j].y,
                          data.point_list()[j].z );
        if ((pos - vpos).length_squared() < GEOMETRY_RESABS*GEOMETRY_RESABS) {
          if (verts[j])
            std::cerr << "Warning: Coincident vertices in surface " << face->id() << std::endl;
          verts[j] = entmap[0][vtx];
          break;
        }
      }
    }
    
      // now create vertices for the remaining points in the facetting
    for (int i = 0; i < data.pointListCount; ++i) {
      if (verts[i]) // if a geometric vertex
        continue;
      double coords[] = { data.point_list()[i].x,
                          data.point_list()[i].y,
                          data.point_list()[i].z };
      rval = iface->create_vertex( coords, verts[i] );
      if (MB_SUCCESS != rval)
        return rval;
    }
    
      // now create facets
    MBRange facets;
    std::vector<MBEntityHandle> corners;
    for (int i = 0; i < data.fListCount; i += data.facet_list()[i]+1) {
      int* facet = data.facet_list() + i;
      corners.resize( *facet );
      for (int j = 1; j <= *facet; ++j) 
        corners[j-1] = verts[facet[j]];
      MBEntityType type;
      if (*facet == 3)
        type = MBTRI;
      else {
        std::cerr << "Warning: non-triangle facet in surface " << face->id() << std::endl;
        if (*facet == 4)
          type = MBQUAD;
        else
          type = MBPOLYGON;
      }
      
      if (surf->bridge_sense() == CUBIT_REVERSED)
        std::reverse( corners.begin(), corners.end() );
      
      MBEntityHandle h;
      rval = iface->create_element( type, &corners[0], corners.size(), h );
      if (MB_SUCCESS != rval)
        return false;
        
      facets.insert( h );
    }
    
      // add vertices and facets to surface set
    rval = iface->add_entities( ci->second, &verts[0], verts.size() );
    if (MB_SUCCESS != rval)
      return false;
    rval = iface->add_entities( ci->second, facets );
    if (MB_SUCCESS != rval)
      return false;
  }
  
  return true;
}

    
