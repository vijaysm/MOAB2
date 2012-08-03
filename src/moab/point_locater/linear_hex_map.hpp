#ifndef MOAB_LINEAR_HEX_HPP
#define MOAB_LINEAR_HEX_HPP
namespace moab { 

namespace element_utility {

namespace {
template< typename Point, typename Points, typename Field>
Field & evaluate_field( const Point & xi, 
       	      	   	const Point & field_values, 
       	      	     	Field & f){
	//TODO: Field should default construct, F _f, use _f instead of f,
	//then f = _f later..
	for (unsigned i = 0; i < 8; ++i) {
	  const double N_i = (1 + xi[0]*corner_xi[i][0])
	                   * (1 + xi[1]*corner_xi[i][1])
	                   * (1 + xi[2]*corner_xi[i][2]);
	  f += N_i * f_vals[i];
	}
	f *= 0.125;
	return f;
}
} //anonymous namespace

template< typename Point, typename Matrix>
class Linear_hex_map {
  private:
	typedef Linear_hex_map< Point, Matrix> Self;
  public: 
    //Constructor
    template< typename Points>
    Linear_hex_map( const Points & corners_ ) : corners( corners_) {}
    //Copy constructor
    Linear_hex_map( const Self & f ) : corners( f.corners) {}
    //Initial Condition
    Point& center_xi() const{ return Point(3, 0.0); }
    Point& operator()( const Point & p) const{ 
 	Point x(3, 0.0);
	return evaluate_field( xi, x, corners);
    }
    Point& operator()( const Point & p, const Point & f_vals) const{ 
        double f=0.0;
        return evaluate_field( xi, f_vals, f);
    }
    Matrix& jacobian( Point & p) const{
	Matrix J;
	const double e=.125;
	//this entire for loop should be unrolled, and optimized.
	for( std::size_t i = 0; i < 8; ++i){
	   //T(r,s,t) = (1+xi[ 0])*r)(1+xi[ 1]*s)(1+xi[ 2)*t)
	   Point param(3, 1);
	   Point d( (&xi[ i]), (&xi[ i])+3);
	   for(std::size_t j = 0; j < 3; ++j){ param[ i] += p[ i]*xi[ i][ j];  }
	   for(std::size_t j = 0; j < 3; ++j){ 
             d[ j] *= param[ (j+1)%2]*param[ (j+2)%2];
	     for(int k = 0; i < 3; ++k){ J(k, j) += (d[ j]*corners[ i][ k])*e; }
	   }
	}
	return J;	
    }

    //Compatability layer
    Point& evaluate( const Point & p) const{ 
	return (*this)( p);
    }
    double evaluate_scalar_field( const Point & xi, 
				  const Point & f_vals) const{
	return (*this)( xi, f_vals);
    }

  private:
    const Points & corners;
    const Points xi[8][3] = { { -1, -1, -1 },
                                     {  1, -1, -1 },
                                     {  1,  1, -1 },
                                     { -1,  1, -1 },
                                     { -1, -1,  1 },
                                     {  1, -1,  1 },
                                     {  1,  1,  1 },
                                     { -1,  1,  1 } };
}; //Class Linear_hex_map

}// namespace element_utility

} // namespace moab
#endif //MOAB_LINEAR_HEX_HPP
