#ifndef MOAB_LINEAR_HEX_HPP
#define MOAB_LINEAR_HEX_HPP

#include "moab/Matrix3.hpp"
#include "moab/CartVect.hpp"

namespace moab { 

namespace element_utility {

namespace {

template< typename Vector>
double normsq( const Vector & v) { 
	const double n= (v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]); 
	std::cout << n << std::endl;
	return n;
}

template< typename Vector>
void vec_subtract( Vector & v, const Vector & u){ 
	for(int i = 0; i < 3; ++i){ v[ i] -= u[ i]; }
}

} //non-exported functionality

template< typename _Matrix>
class Linear_hex_map {
  public:
	typedef _Matrix Matrix;
  private:
	typedef Linear_hex_map< Matrix> Self;
  public: 
    //Constructor
    Linear_hex_map() {}
    //Copy constructor
    Linear_hex_map( const Self & f ) {}

 public:
    //Natural coordinates
    template< typename Entity_handle, typename Points, typename Point>
    std::pair< bool, Point> operator()( const Entity_handle & h, 
					const Points & v, 
					const Point & p, 
					const double tol) const{
	Point result(3, 0.0);
	solve_inverse( p, result, v);
	return std::make_pair( is_contained( result, tol), result);	
    }

  private:

/*    const double reference_points[8][3] = { { -1, -1, -1 },
                                      	    {  1, -1, -1 },
                                      	    {  1,  1, -1 },
                                      	    { -1,  1, -1 },
                                      	    { -1, -1,  1 },
                                      	    {  1, -1,  1 },
                                      	    {  1,  1,  1 },
                                      	    { -1,  1,  1 } }; */


    //This is a hack to avoid a .cpp file and C++11
    //reference_points(i,j) will be a 1 or -1;
    //This should unroll..
    inline const double reference_points( const std::size_t& i, 
          				  const std::size_t& j) const{
    const double rpts[8][3] = { { -1, -1, -1 },
                                {  1, -1, -1 },
                                {  1,  1, -1 },
                                { -1,  1, -1 },
                                { -1, -1,  1 },
                                {  1, -1,  1 },
                                {  1,  1,  1 },
                                { -1,  1,  1 } };
	  return rpts[ i][ j];
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
			const double tol=1.e-6) const {
      const double error_tol_sqr = tol*tol;
      //const double eps = std::numeric_limits<double>::epsilon();
      Point delta(3,0.0);
      xi = delta;
      evaluate( xi, points, delta);
      vec_subtract( delta, x);
      std::size_t num_iterations=0;
      while ( normsq( delta) > error_tol_sqr) {
	std::cout << "Xi: " << std::endl;
	std::cout << xi[ 0] << " " << xi[ 1] << " " << xi[ 2] << std::endl;
	if( ++num_iterations >= 10){ return false; }
        Matrix J;
	jacobian( xi, points, J);
        double det = moab::Matrix::determinant3( J);
        if (det < 1e-10){ 
		std::cout <<"Pt:" << std::endl;
	std::cout << x[ 0] << " " 
		  << x[ 1] << " "
		  << x[ 2] << " " << std::endl;
		std::cout <<"Hex:" << std::endl;
		for(int i = 0; i < 8; ++i){
			std::cout << points[ i][ 0] << " " 
				  << points[ i][ 1] << " "
				  << points[ i][ 2] << " " << std::endl;
		}
	}
		std::cout << "A= [";
		for (int i = 0; i < 3; ++i){
			std::cout << "[";
			for( int j = 0; j < 3; ++j){
				std::cout << J( i, j);
				if (j != 2){ std::cout << ", ";}
			}
			std::cout << "]";
			if( i != 2){
				std::cout << ",";
				std::cout << std::endl;
			}
		}
		std::cout << "]" << std::endl;
	if( det< 1e-10){
		std::cout << "inverse solve failure: det: " << det << std::endl
			  << "number of iterations: " << num_iterations 
			  << std::endl << "error sq: " 
			  << normsq( delta) << std::endl;
		exit( -1);
	}
        vec_subtract( xi, moab::Matrix::inverse(J, 1.0/det) * delta);
        evaluate( xi, points, delta);
	vec_subtract( delta, x);
      }
      std::cout << "---" << std::endl;
       return true;
    }

    template< typename Point, typename Points>
    Point& evaluate( const Point & p, const Points & points, Point & f) const{ 
	typedef typename Points::value_type Vector;
	Vector result;
	for(int i = 0; i < 3; ++i){ result[ i] = 0; }
	for (unsigned i = 0; i < 8; ++i) {
	  const double N_i = (1 + p[0]*reference_points(i,0))
	                   * (1 + p[1]*reference_points(i,1))
	                   * (1 + p[2]*reference_points(i,2));
	    result += N_i * points[ i];
	}
	result *= 0.125;
	for (int i = 0; i < 3; ++i){ f[ i] = result[ i]; }
	return f;
    }

    template< typename Point, typename Points>
    Matrix& jacobian( const Point & p, const Points & points, Matrix & J) const{
	  J = Matrix(0.0);
	  for (std::size_t i = 0; i < 8; ++i) {
	    const double   xi_p = 1 + p[0]*reference_points(i,0);
	    const double  eta_p = 1 + p[1]*reference_points(i,1);
	    const double zeta_p = 1 + p[2]*reference_points(i,2);
	    const double dNi_dxi   = reference_points(i, 0) * eta_p * zeta_p;
	    const double dNi_deta  = reference_points(i, 1) *  xi_p * zeta_p;
	    const double dNi_dzeta = reference_points(i, 2) *  xi_p *  eta_p;
	    J(0,0) += dNi_dxi   * points[i][0];
	    J(1,0) += dNi_dxi   * points[i][1];
	    J(2,0) += dNi_dxi   * points[i][2];
	    J(0,1) += dNi_deta  * points[i][0];
	    J(1,1) += dNi_deta  * points[i][1];
	    J(2,1) += dNi_deta  * points[i][2];
	    J(0,2) += dNi_dzeta * points[i][0];
	    J(1,2) += dNi_dzeta * points[i][1];
	    J(2,2) += dNi_dzeta * points[i][2];
	  }
	  return J *= 0.125;
   }
  private:
    /*
    const double reference_points[8][3] = { { -1, -1, -1 },
                                      	    {  1, -1, -1 },
                                      	    {  1,  1, -1 },
                                      	    { -1,  1, -1 },
                                      	    { -1, -1,  1 },
                                      	    {  1, -1,  1 },
                                      	    {  1,  1,  1 },
                                      	    { -1,  1,  1 } }; */
}; //Class Linear_hex_map

}// namespace element_utility

}// namespace moab
#endif //MOAB_LINEAR_HEX_nPP
