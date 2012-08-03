/** 
 * \class moab::point_locater
 * \author Tim Tautges
 *
 * \brief This class couples data between meshes.
 *
 * The coupler interpolates solution data at a set of points.  Data
 * being interpolated resides on a source mesh, in a tag.
 * Applications calling this coupler send in entities, usually points
 * or vertices, and receive back the tag value interpolated at those
 * points.  Entities in the source mesh containing those points 
 * do not have to reside on the same processor.
 *
 * To use, an application should:
 * - instantiate this coupler by calling the constructor collectively
 *   on all processors in the communicator
 * - call locate_points, which locates the points to be interpolated and
 *   (optionally) caches the results in this object
 * - call interpolate, which does the interpolation
 *
 * Multiple interpolations can be done after locating the points.
 *
 */
#ifndef POINT_LOCATER_HPP
#define POINT_LOCATER_HPP

namespace moab {

template< typename Elements, typename Tree>
class Point_search {

//public types
public:
	typedef typename Tree tree;
	typedef typename Elements::value_type Element;
//private types
private: 
	typedef typename Point_search< Elements, Tree>; 
//public methods
public:

//Constructor
Point_search( Elements & _source_elements,
	      Communicator & _comm,
	      const std::size_t _iterations=3): 
	      tree( _source_elements),
	      source_elements( _source_elements),
	      comm( _comm), boxes(), 
	      iterations( _iterations){}

template< typename Points> 
Error locate_points( Points & query_points, Elements & result,
                     const double rel_eps, const double abs_eps){
	return 0;
} 

template< typename Points> 
Error locate_points( Points & query_points,
                     const double rel_eps = 0.0, 
                     const double abs_eps = 0.0){
     Elements result;
     return locate_points( query_points, result, rel_eps, abs_eps);
} 

//public accessor methods
public:
	Tree &		tree() 		  const   { return tree;	    }
	Boxes &		boxes() 	  const   { return boxes;	    }
	Communicator &	communicator() 	  const   { return comm;	    }
	Elements &	source_elements() const   { return source_elements; }

//private data members  
private:
	Tree & tree;
	Communicator & comm;
	Elements & source_elements;
	Boxes boxes;
	const std::size_t iterations;
};

} // namespace moab

#endif //POINT_LOCATER_HPP
