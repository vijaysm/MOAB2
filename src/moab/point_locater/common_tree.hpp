/** 
 * common_tree.hpp
 * Ryan H. Lewis
 * (C) 2012
 * Functionality common to all trees.
 */
#include <vector>
#include <set>
#include <iostream>
#include <map>
#include <algorithm>
#include <bitset>
#include <numeric>
#include <cmath>
#include <tr1/unordered_map>
#include <limits>

#ifndef COMMON_TREE_HPP
#define COMMON_TREE_HPP

namespace moab {
namespace common_tree { 
#ifdef TREE_DEBUG
template< typename T>
void print_vector( const T & v){
	typedef typename T::const_iterator Iterator;
	std::cout << "[ ";
	for(Iterator i = v.begin(); i != v.end(); ++i){
		std::cout << *i;
		if( i+1 != v.end()){
			std::cout << ", ";
		}
	}
	std::cout << " ]" << std::endl;
}

template< typename T>
void print_vector( const T & begin, const T & end){
	std::cout << "[ ";
	for(T i = begin; i != end; ++i){
		std::cout << (*i)->second.second.to_ulong();
		if( i+1 != end){
			std::cout << ", ";
		}
	}
	std::cout << " ]" << std::endl;
}
#endif

template< typename _Box, typename _Point>
bool box_contains_point(  const _Box & box, const _Point & p){
	for( std::size_t i = 0; i < box.min.size(); ++i){
	     if( p[ i] < box.min[ i] || 
	         p[ i] > box.max[ i]){
	     	return false;
	     }
	 }
	return true;
}

namespace {
	template< typename T> 
	struct Compute_center: public std::binary_function< T, T, T> {
		inline T operator()( const T a, const T b) const{
			return (a+b)/2.0;
		}
	}; //Compute_center
} //non-exported center computation.

class Box{
	public:
	typedef std::vector< double> Vector;
	typedef std::pair< Vector::const_iterator,
			   Vector::const_iterator> Pair;
	Box(): max(3,0.0), min(3,0.0) {}
	Box( const Box & from): max( from.max), min( from.min){}
	template< typename Iterator>
	Box( const Iterator begin):max( begin,begin+3), min(begin,begin+3){}
	Box& operator=( const Box & from){
		max = from.max;
		min = from.min;
		return *this;
	}
	Vector max;
	Vector min;
}; //Box

class Box_with_center{
	public:
	typedef std::vector< double> Vector;
	typedef std::pair< Vector::const_iterator,
			   Vector::const_iterator> Pair;
	Box_with_center(): max(3,0.0), min(3,0.0), center(3,0.0){}
	Box_with_center( const Box_with_center & from): max( from.max), min( from.min), 
	center( from.min){
		std::transform( max.begin(), max.end(), center.begin(),
				center.end(), Compute_center< double>());
	}
	template< typename Iterator>
	Box_with_center( const Iterator begin):max( begin,begin+3), min(begin,begin+3){}
	Box_with_center& operator=( const Box_with_center & from){
		max = from.max;
		min = from.min;
		center = from.center;
		return *this;
	}
	Vector max;
	Vector min;
	Vector center;
}; //Box

template<typename Coordinate, typename Coordinate_iterator>
void update_bounding_max( Coordinate & max, Coordinate_iterator j){
	typedef typename Coordinate::iterator Iterator;
	for( Iterator i = max.begin(); i != max.end(); ++i, ++j){
		*i = std::max( *i, *j);
	}
}

template<typename Coordinate, typename Coordinate_iterator>
void update_bounding_min( Coordinate & min, Coordinate_iterator j){
	typedef typename Coordinate::iterator Iterator;
	for( Iterator i = min.begin(); i != min.end(); ++i, ++j){
		*i = std::min( *i, *j);
	}
}

template< typename Entity_map, typename Ordering>
void construct_ordering( Entity_map & entity_map, Ordering & entity_ordering){
	entity_ordering.reserve( entity_map.size());
	typedef typename Entity_map::iterator Map_iterator;
	for(Map_iterator i = entity_map.begin(); 
			 i != entity_map.end(); 
			 ++i){
		entity_ordering.push_back( i); 
	}
}

//Input: A bunch of entity handles
//Output: A map from handle -> Data 
//Requirements: Data contains at least a bounding box.
//And a non-default constructor which takes only a Box&
template< typename Entity_handles, 
	  typename Element_map, 
	  typename Bounding_box,
	  typename Moab>
void construct_element_map( const Entity_handles & elements, 
			    Element_map & map, 
			    Bounding_box & bounding_box,
			    Moab & moab){
	typedef typename Element_map::mapped_type Box_data;
	typedef typename Entity_handles::value_type Entity_handle;
	typedef typename Entity_handles::iterator Entity_handles_iterator;
	typedef typename std::vector< double> Coordinates;
	typedef typename Coordinates::iterator Coordinate_iterator;
	for( Entity_handles_iterator i = elements.begin(); 
				     i != elements.end(); ++i){
		//TODO: not generic enough. Why dim != 3
		const int DIM = 3;	
		int num_vertices=0;
		//Commence un-necessary deep copying.
		const Entity_handle* vertex_handle;
		moab.get_connectivity( *i, vertex_handle, num_vertices);
		Coordinates coordinate(DIM*num_vertices, 0.0);
		moab.get_coords( vertex_handle, num_vertices, &coordinate[ 0]);
		Box box( coordinate.begin());
		bounding_box = box;
		for( Coordinate_iterator j = coordinate.begin()+DIM; 
				         j != coordinate.end(); j+=DIM){
			update_bounding_max( box.max, j);
			update_bounding_min( box.min, j);
		}
		update_bounding_max( bounding_box.max, box.max.begin());
		update_bounding_min( bounding_box.min, box.min.begin());
		map.insert( std::make_pair( *i, Box_data( box)));
	}
}

} //namspace common_tree

} // namespace moab

#endif //COMMON_TREE_HPP
