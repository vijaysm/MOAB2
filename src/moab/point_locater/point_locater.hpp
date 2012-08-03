/** 
 * point_locater.hpp
 * Ryan H. Lewis 
 * Copyright 2012
 */
#include <vector>
#ifndef POINT_LOCATER_HPP
#define POINT_LOCATER_HPP

namespace moab {

template< typename Elements, 
	  typename _Tree, 
	  typename Communicator>
class Point_search {

//public types
public:
	typedef  _Tree Tree;
	typedef typename Elements::value_type Element;
	//arbitrary type for now.. A box should be a small struct, boxes a
	//Vector of them, or a set?
	typedef	std::vector< std::vector< std::size_t> > Boxes;
	//temporary error code
	typedef typename std::size_t Error;
//private types
private: 
	typedef Point_search< Elements, 
			      Tree, 
			      Communicator> Self; 
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
//private functionality
private:
template< typename Points, typename Bitmask>
void local_point_search( const Points & query_points, Bitmask & located_points, 
			 Elements & result,
			 const double rel_eps, const double abs_eps){}	

//public functionality
public:
template< typename Points> 
Error locate_points( Points & query_points, Elements & result,
                     const double rel_eps, const double abs_eps){
	//bitmask for located points
	std::vector< bool> located_points( query_points.size(), 0);
	local_point_search( query_points, located_points, 
			    result, rel_eps, abs_eps);
	resolve_boxes( query_points, located_points, rel_eps, abs_eps);
	//scatter-gather
	//TODO: probably too good to be true..
	local_point_search( query_points, located_points, 
			    result, rel_eps, abs_eps);
	//send back to target
	return 0;
} 

//public accessor methods
public:
Tree &		kdtree() 	const { return tree;	    	}
Communicator &	communicator() 	const { return comm;	    	}
Elements &	elements() 	const { return source_elements; }

//private data members  
private:
Tree & tree;
Communicator & comm;
Elements & source_elements;
Boxes boxes;
const std::size_t iterations;
}; //class Point_search

} // namespace moab

#endif //POINT_LOCATER_HPP
