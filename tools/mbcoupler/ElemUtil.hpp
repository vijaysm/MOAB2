#ifndef MOAB_ELEM_UTIL_HPP
#define MOAB_ELEM_UTIL_HPP

#include "moab/Core.hpp"
#include "moab/CartVect.hpp"

namespace moab {
namespace ElemUtil {

  bool nat_coords_trilinear_hex(const CartVect* hex_corners, 
                                const CartVect& x, 
                                CartVect& xi,
                                double tol);
  bool point_in_trilinear_hex(const CartVect *hex_corners, 
                              const CartVect& xyz,
                              double etol);
  
  bool point_in_trilinear_hex(const CartVect *hex_corners, 
                              const CartVect& xyz, 
                              const CartVect& box_min, 
                              const CartVect& box_max,
                              double etol);

    //wrapper to hex_findpt
  void nat_coords_trilinear_hex2(const CartVect* hex_corners, 
                                 const CartVect& x, 
                                 CartVect& xi,
                                 double til);



  void hex_findpt(double *xm[3],
                  int n,
                  CartVect xyz, 
                  CartVect& rst,
                  double& dist);

  void hex_eval(double *field,
		int n,
		CartVect rst,
		double &value);

  bool integrate_trilinear_hex(const CartVect* hex_corners,
                               double *corner_fields,
                               double& field_val,
                               int num_pts);

} // namespace ElemUtil
 
} // namespace moab

#endif /*MOAB_ELEM_UTIL_HPP*/
