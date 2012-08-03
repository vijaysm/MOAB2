#ifndef MOAB_ELEMENT_UTILITY_HPP
#define MOAB_ELEMENT_UTILITY_HPP
namespace moab { 

namespace element_utility {
//non-exported functionality
namespace { 
   
   template< typename Point, typename Map>
   bool solve_inverse( const Point & x, Point & xi, 
   		     const Map & map, double tol ) {
     typedef typename Map::Matrix Matrix;
     const double error_tol_sqr = tol*tol;
     //TODO: rename to map.initial( xi);
     xi = map.center_xi();
     Point delta = map( xi) - x;
     Matrix J;
     while (delta % delta > error_tol_sqr) {
       J = map.jacobian( xi);
       double det = J.determinant();
       if (det < std::numeric_limits<double>::epsilon())
         return false;
       //TODO: should be inverse( J)
       xi -= J.inverse(1.0/det) * delta;
       delta = map( xi) - x;
     }
     return true;
   }

   template< typename Points, typename Point, typename Map>
   bool natural_coordinates( const Points& corner_coords,
                             const Point & x, 
   			     const Map & map,
   			     Point & xi, 
   			     double tol ) {
     Point _xi;
     if( solve_inverse( x, _xi, map, tol )){
   	 xi = _xi;
         return fabs(xi[0])-1 < etol && 
   	     	fabs(xi[1])-1 < etol && 
   	     	fabs(xi[2])-1 < etol;
     }
     _xi = xi;
     return false;
   }
} // non-exported functionality

//exported functionality
//TODO: Rewrite all of this create a Parametrizer
template< typename Points, typename Point>
bool natural_coordinates_trilinear_hex( const Points& corner_coords,
                               	        const Point & x, 
					Point & xi, 
					double tol ) {
  typedef typename moab::Matrix3 Matrix;
  moab::element_utility::Linear_hex_map< Point, Matrix> map( corner_coords);
  return natural_coordinates( corner_coords, x, map, tol);
}

template< typename Points, typename Point>
bool natural_coordinates_trilinear_tet( const Points & corner_coords,
                               	        const Point & x, 
					Point & xi, 
					double tol ) {
  typedef typename moab::Matrix3 Matrix;
  moab::element_utility::Linear_tet_map< Point, Matrix> map( corner_coords);
  return natural_coordinates( corner_coords, x, map, tol);
}

}// namespace element_utility
} // namespace moab
#endif //MOAB_ELEMENT_UTILITY_HPP
