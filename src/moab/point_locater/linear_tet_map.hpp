#ifndef MOAB_LINEAR_TET_HPP
#define MOAB_LINEAR_TET_HPP

namespace moab { 
namespace element_utility {

template< typename Points, typename Entity_handle, typename Matrix>
class Linear_tet_map {
  private:
	typedef Linear_tet_map< Points, 
				Entity_handle, 
				Matrix> Self;
  public: 
 	typedef typename Points::value_type Point;
    //Constructor
    Linear_tet_map( const Entity_handle  _eh, 
		    const Points & corners_) : 
		    corners( corners_), eh( _eh){}
    //Copy constructor
    Linear_tet_map( const Self & f ) : T( f.T), Tinv( f.Tinv),
				      det_T( f.det_T), det_Tinv( f.Tinv) {}
    //Initial Condition
    Point& operator()( const Point & p) const{ 
	return Tinv*(p-vertex[ 0]);
    }
    template< typename Points>
    void set_tet( const Entity_handle _eh, const Points & v){
		if (eh != _eh){
		   eh = _eh;
		   corners = v;
		   T = Matrix( v[1][0]-v[0][0], v[2][0]-v[0][0], 
			       v[3][0]-v[0][0],
                    	       v[1][1]-v[0][1], v[2][1]-v[0][1], 
			       v[3][1]-v[0][1],
                    	       v[1][2]-v[0][2], v[2][2]-v[0][2], 
			       v[3][2]-v[0][2]);
		   T_inverse = T.inverse();
		   det_T = T.determinant();
		   det_T_inverse = (det_T == 0: 
					std::numeric_limits< double>.max() : 
					1.0/det_T);
		}
    }
  private:
	const Points & corners;
	const Entity_handle eh;
	const double corner[ 4][ 3] = { {0,0,0},
				        {1,0,0},
				        {0,1,0},
				        {0,0,1} };
	Matrix T, Tinv;
	double det_T, det_Tinv;
	
}; //Class Linear_tet_map

}// namespace element_utility

} // namespace moab
#endif //MOAB_LINEAR_TET_HPP
