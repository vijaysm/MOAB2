#include "MBCore.hpp"
#include "MBCartVect.hpp"

class MBElemUtil {

public:

  static void nat_coords_trilinear_hex(MBCartVect*, 
                                       MBCartVect, 
                                       MBCartVect&,
                                       double);
  static bool point_in_trilinear_hex(MBCartVect *hex, 
                                     MBCartVect xyz,
                                     double etol);
  
  static bool point_in_trilinear_hex(MBCartVect *hex, 
                                     MBCartVect xyz, 
                                     MBCartVect box_min, 
                                     MBCartVect box_max,
                                     double etol);


  static void hex_findpt(double *xm[3],
                         int n,
                         MBCartVect xyz, 
                         MBCartVect& rst,
                         double& dist);
};
