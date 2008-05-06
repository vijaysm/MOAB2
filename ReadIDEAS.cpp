#include <iostream>
#include <fstream>
#include <vector>
#include "assert.h"

#include "ReadIDEAS.hpp"
#include "MBInterface.hpp"
#include "MBInternals.hpp"
#include "MBReadUtilIface.hpp"
#include "FileTokenizer.hpp"
#include "MBRange.hpp"

MBReaderIface* ReadIDEAS::factory( MBInterface* iface )
  { return new ReadIDEAS( iface ); }

ReadIDEAS::ReadIDEAS(MBInterface* impl)
    : MBI(impl)
{
  impl->query_interface("MBReadUtilIface", reinterpret_cast<void**>(&readMeshIface));
}

MBErrorCode ReadIDEAS::load_file(const char* fname, 
                                 MBEntityHandle& meshset, 
                                 const FileOptions& options,
                                 const int* material_set_list,
                                 int num_material_sets ) {

  file.open( fname );

  MBErrorCode rval;

  rval = MBI->create_meshset(MESHSET_SET, mesh_handle);
  if (MB_SUCCESS != rval) return rval;

  char line[10000];
  file.getline(line, 10000);
  std::string s = line;
  if (s.find("-1") > s.length()) return MB_FAILURE;

  while (! file.eof() ) {
    file.getline(line, 10000);
    s = line;

    unsigned int header_id = (unsigned int) strtol(line, NULL, 10);
    switch (header_id) {
      case VERTEX_LIST :
        create_vertices(); 
      break;
      case MAKE_TETRAHEDRA :
        create_tetrahedral_elements();
      break;
      default:
        skip_header(); 
      break;
    }

  }

  meshset = mesh_handle;
  file.close();
  return MB_SUCCESS;

}

void ReadIDEAS::skip_header() {

  // Go until finding a pair of -1 lines
  char *ctmp;
  char line[10000];
  std::string s;

  int end_of_block = 0;

  long int il;

  while (! file.eof() ) {
    file.getline(line, 10000);

    il = std::strtol(line, &ctmp, 10);
    if (il == -1) {
      s = ctmp+1;
      if (s.empty()) end_of_block++;
    }
    else end_of_block = 0;

    if (end_of_block >= 2) break;

  }

}


void ReadIDEAS::create_vertices() {

  // Read two lines: first has some data, second has coordinates
  double *coords;
  char line1[10000], line2[10000];
  int il1, il2;
  char *ctmp1, *ctmp2;
  std::string s1, s2;

  MBErrorCode rval;

  int top_of_block = file.tellg();
  unsigned int num_verts = 0;

  while (! file.eof() ) {

    file.getline(line1, 10000);
    file.getline(line2, 10000);

    // Check if we are at the end of the block
    il1 = std::strtol(line1, &ctmp1, 10);
    il2 = std::strtol(line2, &ctmp2, 10);
    if ((il1 == -1) && (il2 == -1)) {
      s1 = ctmp1+1;
      s2 = ctmp2+1;
      if ((s1.empty()) && (s2.empty())) break;     
    }
    num_verts++;
  }

  file.seekg( top_of_block );
  coords = new double [ 3*num_verts ];

  for (unsigned int i = 0; i < num_verts; i++) {

    file.getline(line1, 10000);
    file.getline(line2, 10000);

    // Get the doubles out of the 2nd line
    coords[3*i  ] = std::strtod(line2, &ctmp2);
    coords[3*i+1] = std::strtod(ctmp2+1, &ctmp2);
    coords[3*i+2] = std::strtod(ctmp2+1, NULL);

  }

  MBRange verts;
  rval = MBI->create_vertices(coords, num_verts, verts);
  assert( MB_SUCCESS == rval );

  rval = MBI->add_entities( mesh_handle, verts );
  assert( MB_SUCCESS == rval );

  file.getline(line1, 10000);
  file.getline(line2, 10000);

}


void ReadIDEAS::create_tetrahedral_elements() {

  MBEntityHandle connect[4];
  char line1[10000], line2[10000];
  int il1, il2;
  char *ctmp1, *ctmp2;
  std::string s1, s2;

  MBErrorCode rval;
  MBEntityHandle handle;

  MBRange verts;
  rval = MBI->get_entities_by_type( mesh_handle, MBVERTEX, verts, true);
  MBEntityHandle vstart = *(verts.begin());
 
  while (! file.eof() ) {

    file.getline(line1, 10000);
    file.getline(line2, 10000);

    // Check if we are at the end of the block
    il1 = std::strtol(line1, &ctmp1, 10);
    il2 = std::strtol(line2, &ctmp2, 10);
    if ((il1 == -1) && (il2 == -1)) {
      s1 = ctmp1+1;
      s2 = ctmp2+1;
      if ((s1.empty()) && (s2.empty())) break;     
    }

    // Get the connectivity out of the 2nd line
    connect[0] = vstart + strtol( line2, &ctmp2, 10) - 1;
    connect[1] = vstart + strtol( ctmp2+1, &ctmp2, 10) - 1;
    connect[2] = vstart + strtol( ctmp2+1, &ctmp2, 10) - 1;
    connect[3] = vstart + strtol( ctmp2+1, &ctmp2, 10) - 1;

    // Make the element
    rval = MBI->create_element(MBTET, connect, 4, handle);
    assert( MB_SUCCESS == rval );

    rval = MBI->add_entities( mesh_handle, &handle, 1);
    assert( MB_SUCCESS == rval );

  }

}
