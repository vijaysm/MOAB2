/* Simple example of use of moab::AdaptiveKDTree class.
   
   Given a hexahedral mesh, find the hexahedron containing each
   input position.
 */


#include "moab/Core.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/Range.hpp"
#include "moab/GeomUtil.hpp"

#include <iostream>
#include <string>

const double EPSILON = 1e-6; // tolerance to use in intersection checks

// Help with error handling.  Given ErrorCode, print
// corresponding string and any available message.
void print_error( moab::Interface& mb, moab::ErrorCode err )
{
  std::string message;
  std::string code;
  if (moab::MB_SUCCESS != mb.get_last_error( message ))
    message.clear();
  code = mb.get_error_string(err);
  std::cerr << "Error: " << code << std::endl;
  if (!message.empty())
    std::cerr << "  " << message << std::endl;
}

// Print diagnostic info for unexpected failures.
#define CHKERR(err) do { if (moab::MB_SUCCESS != (err)) { \
  print_error( mb, (err) ); \
  std::cerr << "Unexpected failure at: " << __FILE__ << ":" << __LINE__ << std::endl; \
  return 2; \
  } } while (false)

// Given an entity set and a point, find the hex contained in the
// entity set which in turn contains the specified point.  Returns
// 0 if point is not in any hexahedron.
moab::EntityHandle hex_containing_point( moab::Interface& mb,
                                         moab::EntityHandle set,
                                         const double point[3] );

// Print hex containing point.
void print_hex( moab::Interface& mb, moab::EntityHandle hex );


int main( )
{
    // Ask user for file to read
  std::string filename;
  std::cout << "Hex mesh file name: ";
  std::cin >> filename;
  
    // Read file into MOAB instance
  moab::ErrorCode rval;
  moab::Core moab;
  moab::Interface& mb = moab;
  rval = mb.load_file( filename.c_str() );
  if (moab::MB_SUCCESS != rval) {
    print_error(mb,rval);
    std::cerr << filename << ": file load failed" << std::endl;
    return 1;
  }
  
    // Get all hex elemeents
  moab::Range elems;
  rval = mb.get_entities_by_type( 0, moab::MBHEX, elems ); CHKERR(rval);
  if (elems.empty()) {
    std::cerr << filename << ": file containd no hexahedra" << std::endl;
    return 1;
  }
  
    // Build a kD-tree from hex elements
  moab::EntityHandle tree_root;
  moab::AdaptiveKDTree tool( &mb );
  rval = tool.build_tree( elems, tree_root ); CHKERR(rval);
  
    // Loop forever (or until EOF), asking user for a point
    // to query and printing the hex element containing that
    // point.
  for (;;) {
    double point[3];
    std::cout << "Point coordinates: ";
    if (!(std::cin >> point[0] >> point[1] >> point[2]))
      break;
  
    moab::EntityHandle leaf;
    rval = tool.leaf_containing_point( tree_root, point, leaf ); CHKERR(rval);
    moab::EntityHandle hex = hex_containing_point( mb, leaf, point );
    if (0 == hex) 
      std::cout << "Point is not contained in any hexahedron." << std::endl;
    else
      print_hex( mb, hex );
  }
  
  return 0;
}

moab::EntityHandle hex_containing_point( moab::Interface& mb,
                                         moab::EntityHandle set,
                                         const double point[3] )
{
  moab::ErrorCode rval;
  moab::CartVect pt(point); // input location
  moab::CartVect coords[8]; // coordinates of corners of hexahedron
  const moab::EntityHandle* conn; // hex connectivity
  int conn_len;
  
    // Get hexes in leaf
  std::vector<moab::EntityHandle> hexes;
  rval = mb.get_entities_by_type( set, moab::MBHEX, hexes ); CHKERR(rval);

    // Check which hex the point is in
  std::vector<moab::EntityHandle>::const_iterator i;
  for (i = hexes.begin(); i != hexes.end(); ++i) {
    rval = mb.get_connectivity( *i, conn, conn_len ); CHKERR(rval);
    rval = mb.get_coords( conn, 8, &coords[0][0] ); CHKERR(rval);
    if (moab::GeomUtil::point_in_trilinear_hex( coords, pt, EPSILON ))
      return *i;
  }
  
    // Return 0 if no hex contains point.
  return 0;  
}

void print_hex( moab::Interface& mb, moab::EntityHandle hex )
{
    // Get MOAB's internal ID for hex element
  int id = mb.id_from_handle(hex);
  
    // Get vertex handles for hex corners
  const moab::EntityHandle* conn; // hex connectivity
  int conn_len;
  mb.get_connectivity( hex, conn, conn_len );
  
    // Get coordinates of vertices
  double coords[3*8]; 
  mb.get_coords( conn, 8, coords );
  
    // Print
  std::cout << " Point is in hex " << id << " with corners: " << std::endl;
  for (int i = 0; i < 8; ++i) {
    std::cout << " (" << coords[3*i] 
              << ", " << coords[3*i+1] 
              << ", " << coords[3*i+2] 
              << ")" << std::endl;
  }
}


  
