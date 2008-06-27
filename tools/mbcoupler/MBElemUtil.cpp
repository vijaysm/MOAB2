#include <iostream>
#include <assert.h>

#include "MBElemUtil.hpp"
#include "MBMatrix3.hpp"
#include "types.h"

namespace MBElemUtil {

/**\brief Class representing a 3-D mapping function (e.g. shape function for volume element) */
class VolMap {
  public:
      /**\brief Return $\vec \xi$ corresponding to logical center of element */
    virtual MBCartVect center_xi() const = 0;
      /**\brief Evaluate mapping function (calculate $\vec x = F($\vec \xi)$ )*/
    virtual MBCartVect evaluate( const MBCartVect& xi ) const = 0;
      /**\brief Evaluate Jacobian of mapping function */
    virtual MBMatrix3 jacobian( const MBCartVect& xi ) const = 0;
      /**\brief Evaluate inverse of mapping function (calculate $\vec \xi = F^-1($\vec x)$ )*/
    bool solve_inverse( const MBCartVect& x, MBCartVect& xi, double tol ) const ;
};

bool VolMap::solve_inverse( const MBCartVect& x, MBCartVect& xi, double tol ) const
{
  const double error_tol_sqr = tol*tol;
  const int max_iterations = 10000;
  int iterations = 0;
  xi = center_xi();
  MBCartVect delta = evaluate(xi) - x;
  MBMatrix3 J;
  while (delta.length_squared() > error_tol_sqr) {
    J = jacobian(xi);
    if (!J.invert() || ++iterations > max_iterations)
      return false;
    xi -= J * delta;
    delta = evaluate( xi ) - x;
  }
  return true;
}

/**\brief Shape function for trilinear hexahedron */
class LinearHexMap : public VolMap {
  public:
    LinearHexMap( const MBCartVect* corner_coords ) : corners(corner_coords) {}
    virtual MBCartVect center_xi() const;
    virtual MBCartVect evaluate( const MBCartVect& xi ) const;
    virtual MBMatrix3 jacobian( const MBCartVect& xi ) const;
  private:
    const MBCartVect* corners;
    static const double corner_xi[8][3];
};

const double LinearHexMap::corner_xi[8][3] = { { -1, -1, -1 },
                                               {  1, -1, -1 },
                                               {  1,  1, -1 },
                                               { -1,  1, -1 },
                                               { -1, -1,  1 },
                                               {  1, -1,  1 },
                                               {  1,  1,  1 },
                                               { -1,  1,  1 } };
MBCartVect LinearHexMap::center_xi() const
  { return MBCartVect(0.0); }

MBCartVect LinearHexMap::evaluate( const MBCartVect& xi ) const
{
  MBCartVect x(0.0);
  for (unsigned i = 0; i < 8; ++i) {
    const double N_i = (1 + xi[0]*corner_xi[i][0])
                     * (1 + xi[1]*corner_xi[i][1])
                     * (1 + xi[2]*corner_xi[i][2]);
    x += N_i * corners[i];
  }
  x *= 0.125;
  return x;
}

MBMatrix3 LinearHexMap::jacobian( const MBCartVect& xi ) const
{
  MBMatrix3 J(0.0);
  for (unsigned i = 0; i < 8; ++i) {
    const double   xi_p = 1 + xi[0]*corner_xi[i][0];
    const double  eta_p = 1 + xi[1]*corner_xi[i][1];
    const double zeta_p = 1 + xi[2]*corner_xi[i][2];
    const double dNi_dxi   = corner_xi[i][0] * eta_p * zeta_p;
    const double dNi_deta  = corner_xi[i][1] *  xi_p * zeta_p;
    const double dNi_dzeta = corner_xi[i][2] *  xi_p *  eta_p;
    J(0,0) += dNi_dxi   * corners[i][0];
    J(1,0) += dNi_dxi   * corners[i][1];
    J(2,0) += dNi_dxi   * corners[i][2];
    J(0,1) += dNi_deta  * corners[i][0];
    J(1,1) += dNi_deta  * corners[i][1];
    J(2,1) += dNi_deta  * corners[i][2];
    J(0,2) += dNi_dzeta * corners[i][0];
    J(1,2) += dNi_dzeta * corners[i][1];
    J(2,2) += dNi_dzeta * corners[i][2];
  }
  return J *= 0.125;
}

bool nat_coords_trilinear_hex( const MBCartVect* corner_coords,
                               const MBCartVect& x,
                               MBCartVect& xi,
                               double tol )
{
  return LinearHexMap( corner_coords ).solve_inverse( x, xi, tol );
}


//
// nat_coords_trilinear_hex2
//  Duplicate functionality of nat_coords_trilinear_hex using hex_findpt
// 
void nat_coords_trilinear_hex2(const MBCartVect hex[8], 
                               const MBCartVect& xyz,
                               MBCartVect &ncoords,
                               double etol)       

{
  const int ndim = 3;
  const int nverts = 8;
  const int vertMap[nverts] = {0,1,3,2, 4,5,7,6}; //Map from nat to lex ordering

  const int n = 2; //linear
  real coords[ndim*nverts]; //buffer

  real *xm[ndim];
  for(int i=0; i<ndim; i++)
    xm[i] = coords + i*nverts;
    
  //stuff hex into coords
  for(int i=0; i<nverts; i++){
    real vcoord[ndim];
    hex[i].get(vcoord);
   
    for(int d=0; d<ndim; d++)
      coords[d*nverts + vertMap[i]] = vcoord[d];
    
  }

  double dist = 0.0;
  MBElemUtil::hex_findpt(xm, n, xyz, ncoords, dist);
  if (3*EPS < dist) {
      // outside element, set extremal values to something outside range
    for (int j = 0; j < 3; j++) {
      if (ncoords[j] < (-1.0-etol) || ncoords[j] > (1.0+etol))
        ncoords[j] *= 10;
    }
  }
  
}

bool point_in_trilinear_hex(const MBCartVect *hex, 
                            const MBCartVect& xyz,
                            double etol) 
{

      const double one = 1.000001;

      MBCartVect  nat(0.);
      nat_coords_trilinear_hex(hex, xyz, nat, etol);
      
      for (unsigned int i = 0; i < 3; i++) {
            if ((nat[i] > one) || (nat[i] < -one)) return false;
      }

      return true;

}


bool point_in_trilinear_hex(const MBCartVect *hex, 
                            const MBCartVect& xyz, 
                            const MBCartVect& box_min, 
                            const MBCartVect& box_max,
                            double etol) 
{

      const double one = 1.000001;

      if ((xyz[0] < box_min[0]) || (xyz[0] > box_max[0])) return false;
      if ((xyz[1] < box_min[1]) || (xyz[1] > box_max[1])) return false;
      if ((xyz[2] < box_min[2]) || (xyz[2] > box_max[2])) return false;

      MBCartVect  nat(0.);
      nat_coords_trilinear_hex(hex, xyz, nat, etol);
      
      for (unsigned int i = 0; i < 3; i++) {
            if ((nat[i] > one) || (nat[i] < -one)) return false;
      }

      return true;

}




// hex_findpt
// Wrapper to James Lottes' findpt routines
// Find the parametric coordinates of an xyz point inside
// a 3d hex spectral element with n nodes per dimension
// xm: coordinates fields, value of x,y,z for each of then n*n*n gauss-lobatto nodes. Nodes are in lexicographical order (x is fastest-changing)
// n: number of nodes per dimension -- n=2 for a linear element
// xyz: input, point to find
// rst: output: parametric coords of xyz inside the element. If xyz is outside the element, rst will be the coords of the closest point
// dist: output: distance between xyz and the point with parametric coords rst
extern "C"{
#include "types.h"
#include "poly.h"
#include "tensor.h"
#include "findpt.h"
#include "extrafindpt.h"
#include "errmem.h"
}

void hex_findpt(real *xm[3],
                int n,
                MBCartVect xyz,
                MBCartVect &rst,
                double &dist)       
{

  //compute stuff that only depends on the order -- could be cached
  real *z[3];
  lagrange_data ld[3];
  opt_data_3 data;

  //triplicates
  for(int d=0; d<3; d++){
    z[d] = tmalloc(real, n);
    lobatto_nodes(z[d], n); 
    lagrange_setup(&ld[d], z[d], n);
  }

  opt_alloc_3(&data, ld);

  //find nearest point
  real x_star[3];
  xyz.get(x_star);

  real r[3] = {0, 0, 0 }; // initial guess for parametric coords
  unsigned c = opt_no_constraints_3;
  dist = opt_findpt_3(&data, (const real **)xm, x_star, r, &c);
  //c tells us if we landed inside the element or exactly on a face, edge, or node

  //copy parametric coords back
  rst = r;

  //Clean-up (move to destructor if we decide to cache)
  opt_free_3(&data);  
  for(int d=0; d<3; ++d) 
    lagrange_free(&ld[d]);
  for(int d=0; d<3; ++d) 
    free(z[d]);
}

} // namespace MBElemUtil
