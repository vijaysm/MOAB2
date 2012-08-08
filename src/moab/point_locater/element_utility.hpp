#ifndef MOAB_ELEMENT_UTILITY_HPP
#define MOAB_ELEMENT_UTILITY_HPP
namespace moab { 

namespace element_utility {
//non-exported functionality

template< std::size_t N, typename Points, typename Point>
bool point_inside_polygon( const Points & vertices, const Point & p){
   double anglesum=0;
   for (std::size_t i=0; i < N; ++i) {
      double costheta=0;
      double x[ 3], y[ 3];
      double mx = 0.0, my = 0.0;
      for(std::size_t j = 0; j < N; ++j){
       	x[ j] = vertices[i][ j] - p [ j];
	mx += x[ j]*x[ j];
        y[ j] = vertices[(i+1)%n][ j] - p[ j];
	my += y[ j]*y[ j];
	costheta += x[ i]*y[ i];
      }
      /* We are on a node, consider this inside */
      if (mx*my <= 1e-12){ return true; }
      costheta /= (mx*my);
      anglesum += acos(costheta);
   }
   return anglesum < PI;
}






























}

   
}// namespace element_utility
} // namespace moab
#endif //MOAB_ELEMENT_UTILITY_HPP
