#ifndef MOAB_SPECTRAL_HEX_HPP
#define MOAB_SPECTRAL_HEX_HPP

#include "moab/Matrix3.hpp"
#include "moab/CartVect.hpp"
extern "C"{
   #include "moab/point_locater/lotte/types.h"
   #include "moab/point_locater/lotte/poly.h"
   #include "moab/point_locater/lotte/tensor.h"
   #include "moab/point_locater/lotte/findpt.h"
   #include "moab/point_locater/lotte/extrafindpt.h"
   #include "moab/point_locater/lotte/errmem.h"
}
#include <sstream>
#include <iomanip>
#include <iostream>

namespace moab { 

namespace element_utility {

namespace {} //non-exported functionality

template< typename _Matrix>
class Spectral_hex_map {
  public:
	typedef _Matrix Matrix;
  private:
	typedef Spectral_hex_map< Matrix> Self;
  public: 
    //Constructor
    Spectral_hex_map() {}
    Spectral_hex_map( const int order,const double * x, 
		      const double * y, const double * z){
	    initialize_spectral_hex(order);
	    _xyz[ 0] = x; _xyz[ 1] = y; _xyz[ 2] = z;
    }
    //Copy constructor
    Spectral_hex_map( const Self & f ) {}
  private:
    void initialize_spectral_hex( int order){
	if (_init && _n==order){ return; }
	if( _init && _n != order){ free_data();}
	_init = true;
	_n = order;
	for( int d = 0; d < 3; d++){
		lobatto_nodes(_z[ d], _n);
		lagrange_setup(&_ld[ d], _z[ d], _n);
	}
	opt_alloc_3(&_data, _ld);
	std::size_t nf = _n*_n, ne = _n, nw = 2*_n*_n + 3*_n;
	_odwork = tmalloc(real, 6*nf + 9*ne + nw);
    }

    void free_data(){
       for(int d=0; d<3; d++){
         free(_z[d]);
         lagrange_free(&_ld[d]);
       }
       opt_free_3(&_data);
       free(_odwork);
     }

 public:
    //Natural coordinates
    template< typename Moab, typename Entity_handle, 
	      typename Points, typename Point>
    std::pair< bool, Point> operator()( const Moab & moab,
					const Entity_handle & h, 
					const Points & v, 
					const Point & p, 
					const double tol=1.e-6) {
	Point result(3, 0.0);
	solve_inverse( p, result, v);
	bool point_found = solve_inverse( p, result, v, tol) && 
						is_contained( result, tol);
	return std::make_pair( point_found, result);
    }

  private:
    void set_gl_points( double * x, double * y, double * z){
	    _xyz[ 0] = x; _xyz[ 1] = y; _xyz[ 2] = z;
    }
    template< typename Point>
    bool is_contained( const Point & p, const double tol) const{
     //just look at the box+tol here
     return ( p[0]>=-1.-tol) && (p[0]<=1.+tol) &&
            ( p[1]>=-1.-tol) && (p[1]<=1.+tol) &&
            ( p[2]>=-1.-tol) && (p[2]<=1.+tol);
    }

    template< typename Point, typename Points>
    bool solve_inverse( const Point & x, 
			Point & xi,
			const Points & points, 
			const double tol=1.e-6) {
      const double error_tol_sqr = tol*tol;
      Point delta(3,0.0);
      xi = delta;
      evaluate( xi, points, delta);
      vec_subtract( delta, x);
      std::size_t num_iterations=0;
      #ifdef SPECTRAL_HEX_DEBUG
 	std::stringstream ss;
	ss << "Point: "; 
       ss << x[ 0 ] << ", " << x[ 1] 
          << ", " << x [ 2] << std::endl;
	ss << "Hex: ";
	for(int i = 0; i < 8; ++i){
 	      	ss << points[ i][ 0] << ", " << points[ i][ 1] << ", "
		   << points[ i][ 2] << std::endl;
	}
	ss << std::endl;
      #endif
      while ( normsq( delta) > error_tol_sqr) {
	#ifdef SPECTRAL_HEX_DEBUG
	ss << "Iter #: "  << num_iterations 
	   << " Err: " << sqrt( normsq( delta)) << " Iterate: ";
	ss << xi[ 0 ] << ", " << xi[ 1] 
		<< ", " << xi[ 2] << std::endl;
	#endif
	if( ++num_iterations >= 5){ return false; }
        Matrix J;
	jacobian( xi, points, J);
        double det = moab::Matrix::determinant3( J);
        if (fabs(det) < 1.e-10){
		#ifdef SPECTRAL_HEX_DEBUG
			std::cerr << ss.str();
		#endif
		#ifndef SPECTRAL_HEX_DEBUG
		std::cerr << x[ 0 ] << ", " << x[ 1] 
			  << ", " << x [ 2] << std::endl;
		#endif
		std::cerr << "inverse solve failure: det: " << det << std::endl;
		exit( -1);
	}
        vec_subtract( xi, moab::Matrix::inverse(J, 1.0/det) * delta);
        evaluate( xi, points, delta);
	vec_subtract( delta, x);
      }
       return true;
    }

    template< typename Point, typename Points>
    Point& evaluate( const Point & p, const Points & points, Point & f) {
    	for(int d = 0; d < 3; ++d){ lagrange_0(&_ld[ d], p[ 0]); }
	for( int d = 0; d < 3; ++d){
		f[ d] = tensor_i3( _ld[ 0].J, _ld[ 0].n,
			_ld[1].J, _ld[1].n,
			_ld[2].J, _ld[2].n,
			_xyz[ d],
			_odwork);
	}
    	return f;
}

    template< typename Point, typename Points>
    Matrix& jacobian( const Point & p, const Points & points, Matrix & J) {
    	real x[ 3];
	for(int i = 0; i < 3; ++i){ _data.elx[ i] = _xyz[ i]; }
	opt_vol_set_intp_3(& _data,x);
	for(int i = 0; i < 9; ++i){ J(i%3, i/3) = _data.jac[ i]; }
	return J;
    }
    
  private:
	bool _init;
	int _n;
	real * _z[ 3];
	lagrange_data _ld[ 3];
	opt_data_3 _data;
	real * _odwork;
	real * _xyz[ 3];
}; //Class Spectral_hex_map

}// namespace element_utility

}// namespace moab
#endif //MOAB_SPECTRAL_HEX_nPP
