#include "MBCore.hpp"
#include "MBCartVect.hpp"

namespace MBElemUtil {

  bool nat_coords_trilinear_hex(const MBCartVect* hex_corners, 
                                const MBCartVect& x, 
                                MBCartVect& xi,
                                double tol);
  bool point_in_trilinear_hex(const MBCartVect *hex_corners, 
                              const MBCartVect& xyz,
                              double etol);
  
  bool point_in_trilinear_hex(const MBCartVect *hex_corners, 
                              const MBCartVect& xyz, 
                              const MBCartVect& box_min, 
                              const MBCartVect& box_max,
                              double etol);

    //wrapper to hex_findpt
  void nat_coords_trilinear_hex2(const MBCartVect* hex_corners, 
                                 const MBCartVect& x, 
                                 MBCartVect& xi,
                                 double til);



  void hex_findpt(double *xm[3],
                  int n,
                  MBCartVect xyz, 
                  MBCartVect& rst,
                  double& dist);
} // namespace MBElemUtil
