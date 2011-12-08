#include <iostream>
#include <limits>
#include <assert.h>

#include "ElemUtil.hpp"
#include "Matrix3.hpp"
#include "types.h"

namespace moab { 
namespace ElemUtil {

  bool debug = false;

/**\brief Class representing a 3-D mapping function (e.g. shape function for volume element) */
class VolMap {
  public:
      /**\brief Return $\vec \xi$ corresponding to logical center of element */
    virtual CartVect center_xi() const = 0;
      /**\brief Evaluate mapping function (calculate $\vec x = F($\vec \xi)$ )*/
    virtual CartVect evaluate( const CartVect& xi ) const = 0;
      /**\brief Evaluate Jacobian of mapping function */
    virtual Matrix3 jacobian( const CartVect& xi ) const = 0;
      /**\brief Evaluate inverse of mapping function (calculate $\vec \xi = F^-1($\vec x)$ )*/
    bool solve_inverse( const CartVect& x, CartVect& xi, double tol ) const ;
};


bool VolMap::solve_inverse( const CartVect& x, CartVect& xi, double tol ) const
{
  const double error_tol_sqr = tol*tol;
  double det;
  xi = center_xi();
  CartVect delta = evaluate(xi) - x;
  Matrix3 J;
  while (delta % delta > error_tol_sqr) {
    J = jacobian(xi);
    det = J.determinant();
    if (det < std::numeric_limits<double>::epsilon())
      return false;
    xi -= J.inverse(1.0/det) * delta;
    delta = evaluate( xi ) - x;
  }
  return true;
}

/**\brief Shape function for trilinear hexahedron */
class LinearHexMap : public VolMap {
  public:
    LinearHexMap( const CartVect* corner_coords ) : corners(corner_coords) {}
    virtual CartVect center_xi() const;
    virtual CartVect evaluate( const CartVect& xi ) const;
    virtual double evaluate_scalar_field( const CartVect& xi, const double *f_vals ) const;
    virtual Matrix3 jacobian( const CartVect& xi ) const;
  private:
    const CartVect* corners;
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
CartVect LinearHexMap::center_xi() const
  { return CartVect(0.0); }

CartVect LinearHexMap::evaluate( const CartVect& xi ) const
{
  CartVect x(0.0);
  for (unsigned i = 0; i < 8; ++i) {
    const double N_i = (1 + xi[0]*corner_xi[i][0])
                     * (1 + xi[1]*corner_xi[i][1])
                     * (1 + xi[2]*corner_xi[i][2]);
    x += N_i * corners[i];
  }
  x *= 0.125;
  return x;
}

double LinearHexMap::evaluate_scalar_field( const CartVect& xi, const double *f_vals ) const
{
  double f(0.0);
  for (unsigned i = 0; i < 8; ++i) {
    const double N_i = (1 + xi[0]*corner_xi[i][0])
                     * (1 + xi[1]*corner_xi[i][1])
                     * (1 + xi[2]*corner_xi[i][2]);
    f += N_i * f_vals[i];
  }
  f *= 0.125;
  return f;
}

Matrix3 LinearHexMap::jacobian( const CartVect& xi ) const
{
  Matrix3 J(0.0);
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


bool nat_coords_trilinear_hex( const CartVect* corner_coords,
                               const CartVect& x,
                               CartVect& xi,
                               double tol )
{
  return LinearHexMap( corner_coords ).solve_inverse( x, xi, tol );
}


//
// nat_coords_trilinear_hex2
//  Duplicate functionality of nat_coords_trilinear_hex using hex_findpt
// 
void nat_coords_trilinear_hex2(const CartVect hex[8], 
                               const CartVect& xyz,
                               CartVect &ncoords,
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
  ElemUtil::hex_findpt(xm, n, xyz, ncoords, dist);
  if (3*EPS < dist) {
      // outside element, set extremal values to something outside range
    for (int j = 0; j < 3; j++) {
      if (ncoords[j] < (-1.0-etol) || ncoords[j] > (1.0+etol))
        ncoords[j] *= 10;
    }
  }
  
}

bool point_in_trilinear_hex(const CartVect *hex, 
                            const CartVect& xyz,
                            double etol) 
{
  CartVect xi;
  return nat_coords_trilinear_hex( hex, xyz, xi, etol )
      && fabs(xi[0])-1 < etol 
      && fabs(xi[1])-1 < etol 
      && fabs(xi[2])-1 < etol;
}


bool point_in_trilinear_hex(const CartVect *hex, 
                            const CartVect& xyz, 
                            const CartVect& box_min, 
                            const CartVect& box_max,
                            double etol) 
{
    // all values scaled by 2 (eliminates 3 flops)
  const CartVect mid = box_max + box_min;
  const CartVect dim = box_max - box_min;
  const CartVect pt = 2*xyz - mid;
  return fabs(pt[0]) - dim[0] < etol &&
         fabs(pt[1]) - dim[1] < etol &&
         fabs(pt[2]) - dim[2] < etol &&
         point_in_trilinear_hex( hex, xyz, etol );
}



// Wrapper to James Lottes' findpt routines
// hex_findpt
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
                CartVect xyz,
                CartVect &rst,
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




// hex_eval
// Evaluate a field in a 3d hex spectral element with n nodes per dimension, at some given parametric coordinates
// field: field values for each of then n*n*n gauss-lobatto nodes. Nodes are in lexicographical order (x is fastest-changing)
// n: number of nodes per dimension -- n=2 for a linear element
// rst: input: parametric coords of the point where we want to evaluate the field
// value: output: value of field at rst

void hex_eval(real *field,
	      int n,
	      CartVect rstCartVec,
	      double &value)       
{
  int d;
  real rst[3];
  rstCartVec.get(rst);

  //can cache stuff below
  lagrange_data ld[3]; 
  real *z[3];
  for(d=0;d<3;++d){
    z[d] = tmalloc(real, n);
    lobatto_nodes(z[d], n);
    lagrange_setup(&ld[d], z[d], n);
  } 

  //cut and paste -- see findpt.c
  const unsigned 
    nf = n*n,
    ne = n,
    nw = 2*n*n + 3*n;
  real *od_work = tmalloc(real, 6*nf + 9*ne + nw);

  //piece that we shouldn't want to cache
  for(d=0; d<3; d++){
    lagrange_0(&ld[d], rst[d]);
  }
  
  value = tensor_i3(ld[0].J,ld[0].n,
		    ld[1].J,ld[1].n,
		    ld[2].J,ld[2].n,
		    field,
		    od_work);

  //all this could be cached
  for(d=0; d<3; d++){
    free(z[d]);
    lagrange_free(&ld[d]); 
  }
  free(od_work);
}


// Gaussian quadrature points for a trilinear hex element.
// Five 2d arrays are defined.
// One for the single gaussian point solution, 2 point solution, 
// 3 point solution, 4 point solution and 5 point solution.
// There are 2 columns, one for Weights and the other for Locations
//                                Weight         Location

const double gauss_1[1][2] = { {  2.0,           0.0          } };

const double gauss_2[2][2] = { {  1.0,          -0.5773502691 },
                               {  1.0         ,  0.5773502691 } };

const double gauss_3[3][2] = { {  0.5555555555, -0.7745966692 },
                               {  0.8888888888,  0.0          },
                               {  0.5555555555,  0.7745966692 } };

const double gauss_4[4][2] = { {  0.3478548451, -0.8611363116 },
                               {  0.6521451549, -0.3399810436 },
                               {  0.6521451549,  0.3399810436 },
                               {  0.3478548451,  0.8611363116 } };

const double gauss_5[5][2] = { {  0.2369268851,  -0.9061798459 },
                               {  0.4786286705,  -0.5384693101 },
                               {  0.5688888889,   0.0          },
                               {  0.4786286705,   0.5384693101 },
                               {  0.2369268851,   0.9061798459 } };

// Function to integrate the field defined by field_fn function
// over the volume of the trilinear hex defined by the hex_corners

bool integrate_trilinear_hex(const CartVect* hex_corners,
                             double *corner_fields,
                             double& field_val,
                             int num_pts)
{
  // Create the LinearHexMap object using the hex_corners array of CartVects
  LinearHexMap hex(hex_corners);

  // Use the correct table of points and locations based on the num_pts parameter
  const double (*g_pts)[2] = 0;
  switch (num_pts) {
  case 1: 
    g_pts = gauss_1;
    break;

  case 2:
    g_pts = gauss_2;
    break;

  case 3:
    g_pts = gauss_3;
    break;

  case 4:
    g_pts = gauss_4;
    break;

  case 5:
    g_pts = gauss_5;
    break;

  default:  // value out of range
    return false;
  }

  // Test code - print Gaussian Quadrature data
  if (debug) {
    for (int r=0; r<num_pts; r++)
      for (int c=0; c<2; c++)
        std::cout << "g_pts[" << r << "][" << c << "]=" << g_pts[r][c] << std::endl;
  }
  // End Test code

  double soln = 0.0;

  for (int i=0; i<num_pts; i++) {     // Loop for xi direction
    double w_i  = g_pts[i][0];
    double xi_i = g_pts[i][1];
    for (int j=0; j<num_pts; j++) {   // Loop for eta direction
      double w_j   = g_pts[j][0];
      double eta_j = g_pts[j][1];
      for (int k=0; k<num_pts; k++) { // Loop for zeta direction
        double w_k    = g_pts[k][0];
        double zeta_k = g_pts[k][1];

        // Calculate the "real" space point given the "normal" point
        CartVect normal_pt(xi_i, eta_j, zeta_k);

        // Calculate the value of F(x(xi,eta,zeta),y(xi,eta,zeta),z(xi,eta,zeta)
        double field = hex.evaluate_scalar_field(normal_pt, corner_fields);

        // Calculate the Jacobian for this "normal" point and its determinant
        Matrix3 J = hex.jacobian(normal_pt);
        double det = J.determinant();

        // Calculate integral and update the solution
        soln = soln + (w_i*w_j*w_k*field*det);
      }
    }
  }

  // Set the output parameter
  field_val = soln;

  return true;
}

} // namespace ElemUtil


namespace Element {


  
  inline const std::vector<CartVect>& Map::get_vertices() {
    return this->vertex;
  }
  //
  void Map::set_vertices(const std::vector<CartVect>& v) {
    if(v.size() != this->vertex.size()) {
      throw ArgError();
    }
    this->vertex = v;
  }// Map::set_vertices()
  //
  CartVect Map::ievaluate(const CartVect& x, double tol, const CartVect& x0) const {
    const double error_tol_sqr = tol*tol;
    double det;
    CartVect xi = x0;
    CartVect delta = evaluate(xi) - x;
    Matrix3 J;

    int iters=0;
    while (delta % delta > error_tol_sqr) {
      if(++iters>50)
        throw Map::EvaluationError();

      J = jacobian(xi);
      det = J.determinant();
      if (det < std::numeric_limits<double>::epsilon())
        throw Map::EvaluationError();
      xi -= J.inverse(1.0/det) * delta;
      delta = evaluate( xi ) - x;
    }
    return xi;
  }// Map::ievaluate()




  const double LinearHex::corner[8][3] = { { -1, -1, -1 },
                                           {  1, -1, -1 },
                                           {  1,  1, -1 },
                                           { -1,  1, -1 },
                                           { -1, -1,  1 },
                                           {  1, -1,  1 },
                                           {  1,  1,  1 },
                                           { -1,  1,  1 } };

  LinearHex::LinearHex() : Map(8) {
    for(unsigned int i = 0; i < 8; ++i) {
      this->vertex[i] = CartVect(this->corner[i][0],this->corner[i][1], this->corner[i][2]);
    }
  }// LinearHex::LinearHex()

  /* For each point, its weight and location are stored as an array.
     Hence, the inner dimension is 2, the outer dimension is gauss_count.
     We use a one-point Gaussian quadrature, since it integrates linear functions exactly.
  */
  const double LinearHex::gauss[1][2] = { {  2.0,           0.0          } };

  CartVect LinearHex::evaluate( const CartVect& xi ) const {
    CartVect x(0.0);
    for (unsigned i = 0; i < 8; ++i) {
      const double N_i = 
        (1 + xi[0]*corner[i][0])
      * (1 + xi[1]*corner[i][1])
      * (1 + xi[2]*corner[i][2]);
      x += N_i * this->vertex[i];
    }
    x *= 0.125;
    return x;
  }// LinearHex::evaluate

  Matrix3 LinearHex::jacobian( const CartVect& xi ) const {
    Matrix3 J(0.0);
    for (unsigned i = 0; i < 8; ++i) {
      const double   xi_p = 1 + xi[0]*corner[i][0];
      const double  eta_p = 1 + xi[1]*corner[i][1];
      const double zeta_p = 1 + xi[2]*corner[i][2];
      const double dNi_dxi   = corner[i][0] * eta_p * zeta_p;
      const double dNi_deta  = corner[i][1] *  xi_p * zeta_p;
      const double dNi_dzeta = corner[i][2] *  xi_p *  eta_p;
      J(0,0) += dNi_dxi   * vertex[i][0];
      J(1,0) += dNi_dxi   * vertex[i][1];
      J(2,0) += dNi_dxi   * vertex[i][2];
      J(0,1) += dNi_deta  * vertex[i][0];
      J(1,1) += dNi_deta  * vertex[i][1];
      J(2,1) += dNi_deta  * vertex[i][2];
      J(0,2) += dNi_dzeta * vertex[i][0];
      J(1,2) += dNi_dzeta * vertex[i][1];
      J(2,2) += dNi_dzeta * vertex[i][2];
    }
    return J *= 0.125;
  }// LinearHex::jacobian()

  double LinearHex::evaluate_scalar_field(const CartVect& xi, const double *field_vertex_value) const {
    double f(0.0);
    for (unsigned i = 0; i < 8; ++i) {
      const double N_i = (1 + xi[0]*corner[i][0])
        * (1 + xi[1]*corner[i][1])
        * (1 + xi[2]*corner[i][2]);
      f += N_i * field_vertex_value[i];
    }
    f *= 0.125;
    return f;
  }// LinearHex::evaluate_scalar_field()

  double LinearHex::integrate_scalar_field(const double *field_vertex_values) const {
    double I(0.0);
    for(unsigned int j1 = 0; j1 < this->gauss_count; ++j1) {
      double x1 = this->gauss[j1][1];
      double w1 = this->gauss[j1][0];
      for(unsigned int j2 = 0; j2 < this->gauss_count; ++j2) {
        double x2 = this->gauss[j2][1];
        double w2 = this->gauss[j2][0];
        for(unsigned int j3 = 0; j3 < this->gauss_count; ++j3) {
          double x3 = this->gauss[j3][1];
          double w3 = this->gauss[j3][0];
          CartVect x(x1,x2,x3);
          I += this->evaluate_scalar_field(x,field_vertex_values)*w1*w2*w3*this->det_jacobian(x);
        }
      }
    }
    return I;
  }// LinearHex::integrate_scalar_field()


  const double LinearTet::corner[4][3] = { {0,0,0},
                                           {1,0,0},
                                           {0,1,0},
                                           {0,0,1}};

  LinearTet::LinearTet() : Map(4) {
    for(unsigned int i = 0; i < 4; ++i) {
      this->vertex[i] = CartVect(this->corner[i][0],this->corner[i][1], this->corner[i][2]);
    }
  }// LinearTet::LinearTet()


  void LinearTet::set_vertices(const std::vector<CartVect>& v) {
    this->Map::set_vertices(v);
    this->T = Matrix3(v[1][0]-v[0][0],v[2][0]-v[0][0],v[3][0]-v[0][0],
                      v[1][1]-v[0][1],v[2][1]-v[0][1],v[3][1]-v[0][1],
                      v[1][2]-v[0][2],v[2][2]-v[0][2],v[3][2]-v[0][2]);
    this->T_inverse = this->T.inverse();
    this->det_T = this->T.determinant();
    this->det_T_inverse = (0.0 == this->det_T ? HUGE : 1.0/this->det_T);
  }// LinearTet::set_vertices()


  double LinearTet::evaluate_scalar_field(const CartVect& xi, const double *field_vertex_value) const {
    double f0 = field_vertex_value[0];
    double f = f0;
    for (unsigned i = 1; i < 4; ++i) {
      f += (field_vertex_value[i]-f0)*xi[i-1];
    }
    return f;
  }// LinearTet::evaluate_scalar_field()

  double LinearTet::integrate_scalar_field(const double *field_vertex_values) const {
    double I(0.0);
    for(unsigned int i = 0; i < 4; ++i) {
      I += field_vertex_values[i];
    }
    I *= this->det_T/24.0;
    return I;
  }// LinearTet::integrate_scalar_field()


}// namespace Element

} // namespace moab
