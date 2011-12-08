#ifndef MOAB_ELEM_UTIL_HPP
#define MOAB_ELEM_UTIL_HPP

#include "moab/Core.hpp"
#include "moab/CartVect.hpp"
#include "Matrix3.hpp"

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

  namespace Element {
    /**\brief Class representing a map (diffeomorphism) F parameterizing a 3D element by its canonical preimage.*/
    /*   
         Shape functions on the element can obtained by a pushforward (pullback by the inverse map) 
         of the shape functions on the canonical element. This is done by extending this class.
         
         We further assume that the parameterizing map is defined by the location of n vertices,
         which can be set and retrived on a Map instance.  The number of vertices is fixed at 
         compile time.
    */
    class Map {
    public:
      /**\brief Construct a Map defined by the given std::vector of vertices. */
      Map(const std::vector<CartVect>& v) {this->vertex.resize(v.size()); this->set_vertices(v);};
      /**\brief Construct a Map defined by n vertices. */
      Map(const unsigned int n) {this->vertex = std::vector<CartVect>(n);};
      /**\brief Evaluate the map on \xi (calculate $\vec x = F($\vec \xi)$ )*/
      virtual CartVect evaluate( const CartVect& xi ) const = 0;
      /**\brief Evaluate the inverse map (calculate $\vec \xi = F^-1($\vec x)$ to given tolerance)*/
      CartVect ievaluate( const CartVect& x, double tol, const CartVect& x0 = CartVect(0.0)) const ;
      /* FIX: should evaluate and ievaluate return both the value and the Jacobian (first jet)? */
      /**\brief Evaluate the map's Jacobi matrix. */
      virtual Matrix3 jacobian( const CartVect& xi ) const = 0;
      /* FIX: should this be evaluated in real coordinates and be obtained as part of a Newton solve? */
      /**\brief Evaluate the inverse of the Jacobi matrix. */
      virtual Matrix3 ijacobian( const CartVect& xi ) const {return this->jacobian(xi).inverse();};
      /* det_jacobian and det_ijacobian should be overriden for efficiency. */
      /**\brief Evaluate the determinate of the Jacobi matrix. */
      virtual double  det_jacobian(const CartVect& xi) const {return this->jacobian(xi).determinant();};
      /* FIX: should this be evaluated in real coordinates and be obtained as part of a Newton solve? */
      /**\brief Evaluate the determinate of the inverse Jacobi matrix. */
      virtual double  det_ijacobian(const CartVect& xi) const {return this->jacobian(xi).inverse().determinant();};

      /**\brief Evaluate a scalar field at a point given field values at the vertices. */
      virtual double   evaluate_scalar_field(const CartVect& xi, const double *field_vertex_values) const = 0;
      /**\brief Integrate a scalar field over the element given field values at the vertices. */
      virtual double   integrate_scalar_field(const double *field_vertex_values) const = 0;

      /**\brief Size of the vertices vector. */
      unsigned int size() {return this->vertex.size();}
      /**\brief Retrieve vertices. */
      inline const std::vector<CartVect>& get_vertices();
      /**\brief Set vertices.      */
      virtual void set_vertices(const std::vector<CartVect>& v);
      
      /* Exception thrown when an evaluation fails (e.g., ievaluate fails to converge). */
      class EvaluationError {
      public:
        EvaluationError(){};
      };// class EvaluationError
      
      /* Exception thrown when a bad argument is encountered. */
      class ArgError {
      public:
        ArgError(){};
      };// class ArgError
    protected:
      std::vector<CartVect> vertex;
    };// class Map

    /**\brief Shape function space for trilinear hexahedron, obtained by a pushforward of the canonical linear (affine) functions. */
    class LinearHex : public Map {
    public:
      LinearHex(const std::vector<CartVect>& vertices) : Map(vertices){};
      LinearHex();
      virtual CartVect evaluate( const CartVect& xi ) const;
      virtual Matrix3  jacobian(const CartVect& xi) const;
      double   evaluate_scalar_field(const CartVect& xi, const double *field_vertex_values) const;
      double   integrate_scalar_field(const double *field_vertex_values) const;
    protected:
      /* Preimages of the vertices -- "canonical vertices" -- are known as "corners". */
      static const double corner[8][3];
      static const double gauss[1][2];
      static const unsigned int corner_count = 8;
      static const unsigned int gauss_count  = 1;
      
    };// class LinearHex

    /**\brief Shape function space for a linear tetrahedron, obtained by a pushforward of the canonical affine shape functions. */
    class LinearTet : public Map {
    public:
      LinearTet(const std::vector<CartVect>& vertices) : Map(vertices){};
      LinearTet();
      /* Override the evaluation routines to take advantage of the properties of P1. */
      virtual CartVect evaluate(const CartVect& xi) const {return this->vertex[0] + this->T*xi;};
      virtual CartVect ievaluate(const CartVect& x) const {return this->T_inverse*(x-this->vertex[0]);};
      virtual Matrix3  jacobian(const CartVect& xi)  const {return this->T;};
      virtual Matrix3  ijacobian(const CartVect& xi) const {return this->T_inverse;};
      virtual double   det_jacobian(const CartVect& xi)  const {return this->det_T;};
      virtual double   det_ijacobian(const CartVect& xi) const {return this->det_T_inverse;};
      //
      double   evaluate_scalar_field(const CartVect& xi, const double *field_vertex_values) const;
      double   integrate_scalar_field(const double *field_vertex_values) const;
      //
      /* Override set_vertices so we can precompute the matrices effecting the mapping to and from the canonical simplex. */
      void     set_vertices(const std::vector<CartVect>& v);
    protected:
      static const double corner[4][3];
      Matrix3 T, T_inverse;
      double  det_T, det_T_inverse;
    };// class LinearTet
    
  }// namespace Element
 
} // namespace moab

#endif /*MOAB_ELEM_UTIL_HPP*/
