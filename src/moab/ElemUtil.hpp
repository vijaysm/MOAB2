#ifndef MOAB_ELEM_UTIL_HPP
#define MOAB_ELEM_UTIL_HPP

#include "moab/CartVect.hpp"
#include <vector>
#include "moab/Matrix3.hpp"

// to access data structures for spectral elements
extern "C"{
#include "types.h"
#include "poly.h"
#include "tensor.h"
#include "findpt.h"
#include "extrafindpt.h"
#include "errmem.h"
}

namespace moab {
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
          Map(const CartVect *v, int num_vs) {this->vertex.resize(num_vs); this->set_vertices(v, num_vs);};
            /**\brief Construct a Map defined by n vertices. */
          Map(const unsigned int n) {this->vertex = std::vector<CartVect>(n);};
            /**\brief Evaluate the map on \params (calculate $\vec x = F($\vec \params)$ )*/
          virtual CartVect evaluate( const CartVect& params ) const = 0;
            /**\brief Evaluate the inverse map (calculate $\vec \params = F^-1($\vec x)$ to given tolerance)*/
          bool evaluate_reverse( const CartVect& x, CartVect &params, double tol = 1.0e-10, const CartVect &x0 = CartVect(0.0)) const ;
            /**\brief decide if within the natural param space, with a tolerance*/
          virtual bool is_inside(const CartVect & params, double tol) const = 0;
            /* FIX: should evaluate and ievaluate return both the value and the Jacobian (first jet)? */
            /**\brief Evaluate the map's Jacobi matrix. */
          virtual Matrix3 jacobian( const CartVect& params ) const = 0;
            /* FIX: should this be evaluated in real coordinates and be obtained as part of a Newton solve? */
          virtual Matrix3 jacobian_inverse(const CartVect& params) {return jacobian(params).inverse();}
            /**\brief Evaluate the determinate of the Jacobi matrix. */
          virtual double  det_jacobian(const CartVect& params) const {return this->jacobian(params).determinant();};
            /* FIX: should this be evaluated in real coordinates and be obtained as part of a Newton solve? */
            /**\brief Evaluate the determinate of the inverse Jacobi matrix. */
          virtual double  det_jacobian_inverse(const CartVect& params) const {return this->jacobian(params).inverse().determinant();};
            /**\brief Evaluate a scalar field at a point given field values at the vertices. */
          virtual double evaluate_scalar(const CartVect& params, const double *field_values) const 
              {double tmp; evaluate_vector(params, field_values, 1, &tmp); return tmp;}
            /**\brief Integrate a scalar field over the element given field values at the vertices. */
          virtual double integrate_scalar(const double *field_values) const
              {double tmp; integrate_vector(field_values, 1, &tmp); return tmp;}
            /**\brief Evaluate a vector field of n tuples at a point given interleaved values at the vertices. */
          virtual void evaluate_vector(const CartVect& params, const double *field_values, int num_tuples, double *eval) const = 0;
            /**\brief Integrate a vector field of n tuples over the element given interleaved values at the vertices. */
          virtual void integrate_vector(const double *field_values, int num_tuples, double *integral) const = 0;

            /**\brief Size of the vertices vector. */
          unsigned int num_verts() {return this->vertex.size();}
            /**\brief Retrieve vertices. */
          inline const std::vector<CartVect>& get_vertices();
            /**\brief Set vertices.      */
          virtual void set_vertices(const CartVect *v, int num_vs);
      
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
          LinearHex(const CartVect *vertices, int num_vs) : Map(vertices, num_vs){};
          LinearHex() : Map(0) {}
          virtual CartVect evaluate( const CartVect& params ) const;
          virtual bool is_inside(const CartVect &params, double tol) const;

          virtual Matrix3 jacobian(const CartVect& params) const;
          void evaluate_vector(const CartVect& params, const double *field_vertex_values, int num_tuples, double *eval) const;
          void integrate_vector(const double *field_vertex_values, int num_tuples, double *integral) const;

      protected:
            /* Preimages of the vertices -- "canonical vertices" -- are known as "corners". */
          static const double corner[8][3];
          static const double gauss[1][2];
          static const unsigned int corner_count = 8;
          static const unsigned int gauss_count  = 1;
      
        };// class LinearHex

          /**\brief Shape function space for trilinear hexahedron, obtained by a pushforward of the canonical linear (affine) functions. */
        class QuadraticHex : public Map {
      public:
          QuadraticHex(const CartVect *vertices, int num_vs) : Map(vertices, num_vs){};
          QuadraticHex();
          virtual CartVect evaluate( const CartVect& params ) const;
          virtual bool is_inside(const CartVect &params, double tol) const;

          virtual Matrix3 jacobian(const CartVect& params) const;
          void evaluate_vector(const CartVect& params, const double *field_values, int num_tuples, double *eval) const;
          void integrate_vector(const double *field_values, int num_tuples, double *integral) const;

      protected:
            /* Preimages of the vertices -- "canonical vertices" -- are known as "corners".
             * there are 27 vertices for a tri-quadratic xes*/
          static const int corner[27][3];
          static const double gauss[8][2];// TODO fix me
          static const unsigned int corner_count = 27;
          static const unsigned int gauss_count  = 8; // TODO fix me

        };// class LinearHex
          /**\brief Shape function space for a linear tetrahedron, obtained by a pushforward of the canonical affine shape functions. */
        class LinearTet : public Map {
      public:
          LinearTet(const CartVect *vertices, int num_vs) : Map(vertices, num_vs){};
          LinearTet();
            /* Override the evaluation routines to take advantage of the properties of P1. */
          virtual CartVect evaluate(const CartVect& params) const {return this->vertex[0] + this->T*params;};
          virtual bool evaluate_reverse(const CartVect& x, CartVect &params, double tol = 1.0e-10, const CartVect& = CartVect(0.0)) const {params = this->T_inverse*(x-this->vertex[0]); return is_inside(params, tol);};
          virtual Matrix3  jacobian(const CartVect&)  const {return this->T;};
          virtual Matrix3  jacobian_inverse(const CartVect&) const {return this->T_inverse;};
          virtual double   det_jacobian(const CartVect&)  const {return this->det_T;};
          virtual double   det_jacobian_inverse(const CartVect&) const {return this->det_T_inverse;};
            //
          void evaluate_vector(const CartVect& params, const double *field_values, int num_tuples, double *eval) const;
          void integrate_vector(const double *field_vertex_values, int num_tuples, double *integral) const;
            //
            /* Override set_vertices so we can precompute the matrices effecting the mapping to and from the canonical simplex. */
          void set_vertices(const CartVect *v, int num_vs);
          bool is_inside(const CartVect &params, double tol) const;
      protected:
          static const double corner[4][3];
          Matrix3 T, T_inverse;
          double  det_T, det_T_inverse;
        };// class LinearTet
    
        class SpectralHex : public Map {
      public:
          SpectralHex(const CartVect *vertices, int num_vs) : Map(vertices, num_vs){};
          SpectralHex(int order, double * x, double * y, double *z) ;
          SpectralHex(int order);
          SpectralHex();
          virtual ~SpectralHex();
          void set_gl_points( double * x, double * y, double *z) ;
          virtual CartVect evaluate( const CartVect& params ) const;
          virtual bool evaluate_reverse(const CartVect& x, CartVect &params, double tol = 1.0e-10, const CartVect &init = CartVect(0.0)) const;
          virtual Matrix3  jacobian(const CartVect& params) const;
          void evaluate_vector(const CartVect& params, const double *field_values, int num_tuples, double *eval) const;
          void integrate_vector(const double *field_values, int num_tuples, double *integral) const;
          bool is_inside(const CartVect &params, double tol) const;

            // to compute the values that need to be cached for each element of order n
          void Init(int order);
          void freedata();
      protected:
            /* values that depend only on the order of the element , cached */
            /*  the order in 3 directions */
          static int _n;
          static real *_z[3];
          static lagrange_data _ld[3];
          static opt_data_3 _data;
          static real * _odwork;// work area
// flag for initialization of data
          static bool _init;

          real * _xyz[3];


        };// class SpectralHex

        class SpectralQuad : public Map {
      public:
          SpectralQuad(const CartVect *vertices, int num_vs) : Map(vertices, num_vs){};
          SpectralQuad(int order, double * x, double * y, double *z) ;
          SpectralQuad(int order);
          SpectralQuad();
          virtual ~SpectralQuad();
          void set_gl_points( double * x, double * y, double *z) ;
          virtual CartVect evaluate( const CartVect& params ) const;// a 2d, so 3rd component is 0, always
          virtual bool evaluate_reverse(const CartVect& x, CartVect &params, double tol = 1.0e-10, const CartVect &x0 = CartVect(0.0)) const; //a 2d, so 3rd component is 0, always
          virtual Matrix3  jacobian(const CartVect& params) const;
          void evaluate_vector(const CartVect& params, const double *field_vertex_values, int num_tuples, double *eval) const;
          void integrate_vector(const double *field_vertex_values, int num_tuples, double *integral) const;
          bool is_inside(const CartVect &params, double tol) const;

            // to compute the values that need to be cached for each element of order n
          void Init(int order);
          void freedata();
            // this will take node, vertex positions and compute the gl points
          void compute_gl_positions();
          void get_gl_points(  double *& x, double *& y, double *& z, int & size) ;
      protected:
            /* values that depend only on the order of the element , cached */
            /*  the order in all 3 directions ; it is also np in HOMME lingo*/
          static int _n;
          static real *_z[2];
          static lagrange_data _ld[2];
          static opt_data_2 _data; // we should use only 2nd component
          static real * _odwork;// work area

            // flag for initialization of data
          static bool _init;
          static real * _glpoints; // it is a space we can use to store gl positions for elements
            // on the fly; we do not have a tag yet for them, as in Nek5000 application
            // also, these positions might need to be moved on the sphere, for HOMME grids
            // do we project them or how do we move them on the sphere?

          real * _xyz[3]; // these are gl points; position?


        };// class SpectralQuad


  inline const std::vector<CartVect>& Map::get_vertices() {
    return this->vertex;
  }
  //
        inline void Map::set_vertices(const CartVect *v, int num_vs) {
          if(num_vs != (int)this->vertex.size()) {
      throw ArgError();
    }
          std::copy(v, v+num_vs, this->vertex.begin());
  }// Map::set_vertices()


    }// namespace Element
 
} // namespace moab

#endif /*MOAB_ELEM_UTIL_HPP*/
