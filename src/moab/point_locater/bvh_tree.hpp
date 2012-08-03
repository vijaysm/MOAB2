/** 
 * bvh_tree.hpp
 * Ryan H. Lewis
 * (C) 2012
 *
 * An element tree partitions a mesh composed of elements.
 * We subdivide the bounding box of a mesh, by putting boxes
 * on the left if there center is on the left of a split line
 * and vice versa.
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
#include "common_tree.hpp"

#ifndef BVH_TREE_HPP
#define BVH_TREE_HPP
namespace moab {

//forward declarations
template< typename _Entity_handles, 
	  typename _Box, 
	  typename _Moab,
	  typename _Parametrizer> class Bvh_tree;

//non-exported functionality
namespace {
	namespace _bvh {
	template< typename Box, typename Entity_handle>
	struct _Node{
		unsigned int dim;
		unsigned int child;
		float Lmax, Rmin;
		std::vector< std::pair< Box, Entity_handle> > entities;
	}; // _Node

	struct _Split_data {
	
	}; //_Split_data


	} // namespace _bvh
} //private namespace

template< typename _Entity_handles,
	  typename _Box,
	  typename _Moab,
	  typename _Parametrizer>
class Bvh_tree {
//public types
public:
	typedef  _Entity_handles Entity_handles;
	typedef  _Box Box;
	typedef  _Moab Moab;
	typedef  _Parametrizer Parametrizer;
	typedef typename Entity_handles::value_type Entity_handle;
	
//private types
private: 
	typedef Bvh_tree< _Entity_handles, 
			      _Box, 
			      _Moab,
			      _Parametrizer> Self;
	typedef typename std::pair< Box, Entity_handle> Leaf_element;
	typedef _bvh::_Node< Box, Entity_handle> Node;
	typedef typename std::vector< Node> Nodes;
//public methods
public:
//Constructor
Bvh_tree( Entity_handles & _entities, 
	  Moab & _moab, 
	  Box & _bounding_box, 
	  Parametrizer & _entity_contains): entity_handles_( _entities), 
				tree_(), moab( _moab), 
				bounding_box( _bounding_box),
				entity_contains( entity_contains){
	typedef typename Entity_handles::iterator Entity_handle_iterator;
	typedef typename std::map< Entity_handle, Box> Entity_map;
	typedef typename Entity_map::iterator Entity_map_iterator;
	typedef std::vector< Entity_map_iterator> Vector;
	//a fully balanced tree will have 2*_entities.size()
	//which is one doubling away..
	tree_.reserve( entity_handles_.size());
	Entity_map entity_map;
	common_tree::construct_element_map( entity_handles_, 
					    entity_map, 
					    bounding_box, 
					    moab);
 	_bounding_box = bounding_box;
	
	Vector entity_ordering;
	construct_ordering( entity_map, entity_ordering); 
	//We only build nonempty trees
	if( entity_ordering.size()){ 
	 //initially all bits are set
	 tree_.push_back( Node());
	 _bvh::_Split_data data;
	 const int depth = build_tree( entity_ordering.begin(), 
				       entity_ordering.end(), 0, data);
	 std::cout << "max tree depth: " << depth << std::endl; 
	}
}

//Copy constructor
Bvh_tree( Self & s): entity_handles_( s.entity_handles_), 
			 tree_( s.tree_), moab( s.moab), 
			 bounding_box( s.bounding_box),
			 entity_contains( s.entity_contains){}

//private functionality
private:
template< typename Iterator, typename Split_data>
int build_tree( Iterator begin, Iterator end, 
		const int index, Split_data & data){
	//find_best_split( begin, end, data);
	//std::sort( begin, end, Comparator());
	return 0;	
}

template< typename Vector, typename Node_index>
Entity_handle _find_point( const Vector & point, 
			   const Node_index & index) const{
	typedef typename Node::Entities::const_iterator Entity_iterator;
	const Node & node = tree_[ index];
	if( node.is_leaf()){
		//check each node
		for( Entity_iterator i = node.entities.begin(); 
				     i != node.entities.end(); ++i){
			if( common_tree::box_contains_point( i->first, point) &&
				entity_contains( i->second, point)){
				return i->second;
			}
		}
		return 0;
	}
	if( point[ node.dim] < node.left_line){
		return _find_point( point, node.left);
	}else if( point[ node.dim] > node.right_line){
		return _find_point( point, node.left+1);
	} else {
		const Entity_handle result =  _find_point( point, node.left);
		if( result != 0){ return result; }
		return _find_point( point, node.right+1);
	}
}

//public functionality
public:
template< typename Vector>
Entity_handle find( const Vector & point) const{
	typedef typename Vector::const_iterator Point_iterator;
	typedef typename Box::Pair Pair; 
	typedef typename Pair::first_type Box_iterator;
	return  _find_point( point, 0);
}

//public accessor methods
public:

//private data members  
private:
	const Entity_handles & entity_handles_;
	Nodes tree_;
	Moab & moab;
	Box bounding_box;
	const Parametrizer & entity_contains;

}; //class Bvh_tree

} // namespace moab

#endif //BVH_TREE_HPP
