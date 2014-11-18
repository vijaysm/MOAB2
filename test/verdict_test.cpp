#include "moab/Core.hpp"
#include "moab/CartVect.hpp"
#include "moab/Range.hpp"
#include "moab/verdict.h"
#include <iostream>
#include <iomanip>
#include <cstdio>
#include "TestUtil.hpp"


std::string TestDir( STRINGIFY(MESHDIR) );

std::string filename = TestDir + "/mbtest1.g";

using namespace moab;

int main( int argc, char* argv[] )
{
  ErrorCode rval;
  Core moab_core;
  Interface* interface = &moab_core;
  if(argc > 1) if (argc > 1) filename = std::string(argv[1]);
  rval = interface->load_mesh( filename.c_str());
  if (MB_SUCCESS != rval) {
      std::cerr << "Error reading file: " << filename.c_str() << std::endl;
      exit(2);
    }
  //std::cout << "loaded mesh file: " << filename.c_str() << std::endl;
  Range elems;
  interface->get_entities_by_dimension( 0, 3, elems );
  //std::cout << "# 3 dim entities = " << elems.size() << std::endl;

  double t_volume = 0, h_volume = 0;
  double edge_ratio_max=0.;
  double edge_ratio_max_tet = 0;
  for (Range::const_iterator ei = elems.begin(); ei != elems.end(); ++ei)
    {
      std::vector<EntityHandle> connect;
      double t_verdict_coords [4][3];
      double h_verdict_coords [8][3];

      rval = interface->get_connectivity( &(*ei), 1, connect );
      if (MB_SUCCESS != rval) return rval;
      for(int i=0; i< (int) connect.size(); i++){
          // get volume of each element and add
          double coords[3];
          rval = interface->get_coords(&connect[i], 1, coords);

          if( connect.size() == 4){
              t_verdict_coords[i][0] = coords[0];
              t_verdict_coords[i][1] = coords[1];
              t_verdict_coords[i][2] = coords[2];
            }
          if( connect.size() == 8){
              h_verdict_coords[i][0] = coords[0];
              h_verdict_coords[i][1] = coords[1];
              h_verdict_coords[i][2] = coords[2];
            }
        }
      if( connect.size() == 4)
      {
        t_volume+=v_tet_volume(4, t_verdict_coords);
        double edge_ratio = v_tet_edge_ratio( 8, t_verdict_coords );
        if (edge_ratio>edge_ratio_max_tet)
          edge_ratio_max_tet = edge_ratio;
      }
      if( connect.size() == 8)
      {
        h_volume+=v_hex_volume (8, h_verdict_coords);
        double edge_ratio = v_hex_edge_ratio( 8, h_verdict_coords );
        if (edge_ratio>edge_ratio_max)
          edge_ratio_max = edge_ratio;
      }
    }
  if (argc == 1 && (fabs(t_volume-1) > 1.e-10 || fabs(h_volume - 1.)> 1.e-10 ) )
  {
    std::cerr << "Test failed, volume of tet or hex elements not equal to 1" << std::endl;
    return 1;
  }
  std::cout << " edge_ratio_max hex: " << edge_ratio_max << "\n";
  std::cout << " edge_ratio_max tet: " << edge_ratio_max_tet << "\n";

  return 0;
}

