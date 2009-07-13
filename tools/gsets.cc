#include <iostream>
#include <stdlib.h>
#include "MBRange.hpp"
#include "MBCore.hpp"
#include "MBTagConventions.hpp"
#include "GeomTopoTool.hpp"

MBTag geomTag = 0;
MBTag blockTag = 0;
MBTag sideTag = 0;
MBTag nodeTag = 0;
MBTag nameTag = 0;
MBTag idTag = 0;
bool printAnonSets = false;
bool printSVSense = false;

MBCore mb;
GeomTopoTool geomTool(&mb);

void usage( const char* name, bool brief = true ) {
  std::ostream& str = brief ? std::cerr : std::cout;
  if (!brief)
   str << name << ": A tool to export entity set parent/child relations" << std::endl
       << "      for use as input to graphviz" << std::endl;
  str << "Usage: " << name << " [-a | [-g] [-m] [-n] ] <input_file>" << std::endl
      << "       " << name << " -h" << std::endl;
  if (brief)
    exit(1);
  str << "  The default behavior is equivalent to \"-gmn\"." << std::endl
      << "  If any of the following options are used to specify which " << std::endl
      << "  sets to output, then there are no defaults.  Only the " << std::endl
      << "  indicated sets will be output." << std::endl
      << "  -a  : write all sets (default is only geom, mesh, and named)" << std::endl
      << "  -g  : write geometric topology sets" << std::endl
      << "  -m  : write material sets and boundary condition sets" << std::endl
      << "  -n  : write named sets" << std::endl
      << "  -s  : label surface-volume links with sense" << std::endl
      << "  The default link behavior is to both child links" << std::endl
      << "  and containment with solid lines." << std::endl
      << "  -P  : do not write child links" << std::endl
      << "  -p  : write child links with dashed lines" << std::endl
      << "  -C  : do not write containment links" << std::endl
      << "  -c  : write containment links with dashed lines" << std::endl;
  exit(0);
}

enum Link { NONE = 0, SOLID, DASHED };
void write_dot( Link contained, Link children );
void dot_nodes( std::ostream& s, MBRange& sets_out );
void dot_children( std::ostream& s, const MBRange& sets, bool dashed );
void dot_contained( std::ostream& s, const MBRange& sets, bool dashed );

int main( int argc, char* argv[] )
{
  Link children = SOLID, contained = SOLID;
  bool printGeomSets = true;
  bool printMeshSets = true;
  bool printNamedSets = true;
  const char* input_file = 0;
  bool geom_flag = false, mesh_flag = false, name_flag = false, all_flag = false;
  bool no_more_flags = false;
  for (int i = 1; i < argc; ++i) {
    if (no_more_flags || argv[i][0] != '-') {
      if (input_file) 
        usage(argv[0]);
      input_file = argv[i];
      continue;
    }
    for (int j = 1; argv[i][j]; ++j) {
      switch (argv[i][j]) {
        case 'a': all_flag = true; break;
        case 'g': geom_flag = true; break;
        case 'm': mesh_flag = true; break;
        case 'n': name_flag = true; break;
        case 's': printSVSense = true; break;
        case 'P': children = NONE; break;
        case 'p': children = DASHED; break;
        case 'C': contained = NONE; break;
        case 'c': contained = DASHED; break;
        case '-': no_more_flags = true; break;
        case 'h': usage(argv[0], false);
        default:
          std::cerr << "Unknown flag: '" << argv[i][j] << "'" << std::endl;
          usage(argv[0]);
      }
    }
  }
  
  if (!input_file) {
    std::cerr << "No input file specified." << std::endl;
    usage(argv[0]);
  }
  
  if (all_flag) {
    printGeomSets = printMeshSets = printNamedSets = printAnonSets = true;
  }
  else if (geom_flag || mesh_flag || name_flag) {
    printGeomSets = geom_flag;
    printMeshSets = mesh_flag;
    printNamedSets = name_flag;
  }
  
  if (MB_SUCCESS != mb.load_mesh( input_file )) {
    std::cerr << input_file << ": file read failed." << std::endl;
    return 1;
  }

  MBTag t;
  if (printGeomSets) {
    if (MB_SUCCESS == mb.tag_get_handle( GEOM_DIMENSION_TAG_NAME, t )) {
      geomTag = t;
    }
  }
  if (printMeshSets) {
    if (MB_SUCCESS == mb.tag_get_handle( MATERIAL_SET_TAG_NAME, t )) {
      blockTag = t;
    }
    if (MB_SUCCESS == mb.tag_get_handle( DIRICHLET_SET_TAG_NAME, t )) {
      nodeTag = t;
    }
    if (MB_SUCCESS == mb.tag_get_handle( NEUMANN_SET_TAG_NAME, t )) {
      sideTag = t;
    }
  }
  if (printNamedSets) {
    if (MB_SUCCESS == mb.tag_get_handle( NAME_TAG_NAME, t )) {
      nameTag = t;
    }
  }
  if (MB_SUCCESS == mb.tag_get_handle( GLOBAL_ID_TAG_NAME, t ))
    idTag = t;
  
  write_dot( contained, children );
  return 0;
}
      
  
void write_dot( Link contained, Link children )
{
  MBRange sets;
  std::cout << "digraph {" << std::endl;
  dot_nodes( std::cout, sets );
  std::cout << std::endl;
  if (contained)
    dot_contained( std::cout, sets, contained == DASHED );
  if (children)
    dot_children( std::cout, sets, children == DASHED );
  std::cout << "}" << std::endl;
}

void dot_get_sets( MBRange& curr_sets, MBRange& result_sets, 
                   MBTag tag, void* tag_val = 0 )
{
  if (!tag)
    return;
    
  result_sets.clear();
  mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &tag, &tag_val, 1, result_sets );
  result_sets = subtract( result_sets,  curr_sets );
  curr_sets.merge( result_sets );
}

void dot_write_node( std::ostream& s, MBEntityHandle h, 
                     const char* label, int* id = 0 )
{
  s << 's' << mb.id_from_handle(h) << " [label = \"" << label;
  if (id)
    s << ' ' << *id;
  s << "\"];" << std::endl;
}

void dot_write_id_nodes( std::ostream& s,
                         const MBRange& entites,
                         MBTag id_tag,
                         const char* type_name )
{
  int id;
  for (MBRange::iterator i = entites.begin(); i != entites.end(); ++i)
    if (MB_SUCCESS == mb.tag_get_data( id_tag, &*i, 1, &id )) 
      dot_write_node( s, *i, type_name, &id );
}

void dot_nodes( std::ostream& s, MBRange& sets )
{
  MBRange vol_sets, surf_sets, curv_sets, vert_sets;
  MBRange block_sets, side_sets, node_sets;
  MBRange named_sets, other_sets;

  dot_get_sets( sets, named_sets, nameTag );
    
  int dim = 3;
  dot_get_sets( sets, vol_sets, geomTag, &dim );
  dim = 2;
  dot_get_sets( sets, surf_sets, geomTag, &dim );
  dim = 1;
  dot_get_sets( sets, curv_sets, geomTag, &dim );
  dim = 0;
  dot_get_sets( sets, vert_sets, geomTag, &dim );
  
  dot_get_sets( sets, block_sets, blockTag );
  dot_get_sets( sets, side_sets, sideTag );
  dot_get_sets( sets, node_sets, nodeTag );

  if (printAnonSets) {
    mb.get_entities_by_type( 0, MBENTITYSET, other_sets );
    MBRange xsect = subtract( other_sets,  sets );
    sets.swap(other_sets);
    other_sets.swap(xsect);
  }
  
  dot_write_id_nodes( s, vol_sets , idTag, "Volume"  );
  dot_write_id_nodes( s, surf_sets, idTag, "Surface" );
  dot_write_id_nodes( s, curv_sets, idTag, "Curve"   );
  dot_write_id_nodes( s, vert_sets, idTag, "Vertex"  );
  dot_write_id_nodes( s, block_sets, blockTag, "Block" );
  dot_write_id_nodes( s, side_sets, sideTag, "Neumann Set" );
  dot_write_id_nodes( s, node_sets, nodeTag, "Dirichlet Set" );
  
  MBRange::iterator i;
  char name[NAME_TAG_SIZE+1];
  for (i = named_sets.begin(); i != named_sets.end(); ++i) {
    if (MB_SUCCESS == mb.tag_get_data( nameTag, &*i, 1, name )) {
      name[NAME_TAG_SIZE] = '\0';
      dot_write_node( s, *i, name );
    }
  }
  for (i = other_sets.begin(); i != other_sets.end(); ++i) {
    int id = mb.id_from_handle(*i);
    dot_write_node( s, *i, "EntitySet ", &id );
  }
}

void dot_down_link( std::ostream& s,
                    MBEntityHandle parent,
                    MBEntityHandle child,
                    bool dashed,
                    const char* label = 0 )
{
  s << 's' << mb.id_from_handle(parent) << " -> "
    << 's' << mb.id_from_handle(child);
  if (dashed && label) 
    s << " [style = dashed label = \"" << label << "\"]";
  else if (dashed)
    s << " [style = dashed]";
  else if (label)
    s << " [label = \"" << label << "\"]";
  s << ';' << std::endl;
}


void dot_children( std::ostream& s, 
                   const MBRange& sets,
                   bool dashed )
{
  bool forward;
  const char *fstr = "forward", *rstr = "reverse";
  for (MBRange::iterator i = sets.begin(); i != sets.end(); ++i) {
    MBRange parents;
    mb.get_parent_meshsets( *i, parents );
    parents = intersect( parents,  sets );
    
    for (MBRange::iterator j = parents.begin(); j != parents.end(); ++j) {
      const char* linklabel = 0;
      if (printSVSense && MB_SUCCESS == geomTool.get_sense( *i, *j, forward ))
        linklabel = forward ? fstr : rstr;
      dot_down_link( s, *j, *i, dashed, linklabel);
    }
  }
}


void dot_contained( std::ostream& s, 
                    const MBRange& sets,
                    bool dashed )
{
  for (MBRange::iterator i = sets.begin(); i != sets.end(); ++i) {
    MBRange contained;
    mb.get_entities_by_type(*i, MBENTITYSET, contained );
    contained = intersect( contained,  sets );
    
    for (MBRange::iterator j = contained.begin(); j != contained.end(); ++j)
      dot_down_link( s, *i, *j, dashed );
  }
}
