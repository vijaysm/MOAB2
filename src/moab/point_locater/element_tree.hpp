/** 
 * element_tree.hpp
 * Ryan H. Lewis
 * (C) 2012
 *
 * An element tree partitions a mesh composed of elements.
 * We subdivide the bounding box of a mesh, and each element is
 * either entirely on the left, entirely on the right, or crossing
 * the diving line. We build a tree on the mesh with this property.
 */
#include <vector>
#include <iostream>
#include <map>
#include <algorithm>

#ifndef ELEMENT_TREE_HPP
#define ELEMENT_TREE_HPP

namespace moab {
//forward declarations

template< typename _Entity_handles, 
	  typename _Boxes, 
	  typename _Moab> class Element_tree;

//non-exported functionality
namespace _element_tree {

template< typename _Entity_handles>
class Node{
	//public types:
	public:
	typedef _Entity_handles Entity_handles;

	//private types:
	private:
	typedef Node< _Entity_handles> Self;

	//Constructors
	public:
	//Default constructor
	Node(): left_( -1), middle_( -1), right_( -1),
	        dim( -1), split( 0), left_line( 0), right_line( 0), 
		entity_handles_( NULL) {}

	//Destructor
	~Node(){ delete entity_handles_; }

	//Constructor (leaf)
	Node( Entity_handles & _entities): 
		left_( -1), middle_( -1), right_( -1),
	        dim( -1), split( 0), left_line( 0), right_line( 0), 
		entity_handles_( &_entities) {}

	//Copy constructor
	Node( const Self & from): 
	left_( from.left_), middle_( from.middle_), right_( from.right_),
	dim( from.dim), split( from.split), left_line( from.left_line), 
	right_line( from.right_line), entity_handles_( from.entity_handles_) {}

	// Functionality
	public: 
	bool leaf() const { return left_ == 0 && middle_ == 0 && right_ == 0; }
	
	//private data members:
	private:
	int  left_;
	int  middle_;
	int  right_;
	int dim;
	double split;
	double left_line;
	double right_line;
	Entity_handles * entity_handles_;

	//Element_tree can touch my privates.
	template< Entity_handles, typename B> friend class moab::Element_tree;
}; //class Node

struct Partition_data{
	std::size_t left;
	std::size_t middle;
	std::size_t right;
}; // Partition_data

template< typename Iterator>
struct Element_compare: public std::binary_function< Iterator, Iterator, bool>{
	bool operator()( const Iterator a, const Iterator b){
		return (*a).first < (*b).first;
	}
}; // Element_compae



} //namespace _element_tree

template< typename _Entity_handles,
	  typename _Boxes,
	  typename _Moab>
class Element_tree {

//public types
public:
	typedef  _Entity_handles Entity_handles;
	typedef  _Boxes Boxes;
	typedef  _Moab Moab;
	typedef typename Entity_handles::value_type Element;
	typedef typename std::pair< bool, Element> mapped_type;
	
//private types
private: 
	typedef Element_tree< _Entity_handles, 
			      _Boxes, 
			      _Moab> Self; 
	typedef typename _element_tree::Node< Entity_handles> Node;
	typedef typename Boxes::value_type Box;
	typedef typename std::pair< Box, std::size_t> Data;
	typedef typename Entity_handles::value_type Entity_handle;
	//TODO: we really want an unordered map here..
	typedef typename std::map< Entity_handle, Data> Entity_map;
	typedef typename Entity_map::const_iterator Entity_map_iterator;
	typedef typename std::vector< Entity_map_iterator> Element_list;
//public methods
public:
//Constructor
Element_tree( Entity_handles & _entities, Moab & _moab): 
	entity_handles_( _entities), tree_(), moab( _moab) {
	tree_.reserve( _entities.size());
	Entity_map map;
	construct_element_map( entity_handles_, map);
	Element_list elements( map.size());
	std::size_t index = 0;
	for(typename Entity_map::const_iterator i = map.begin(); 
					 i != map.end(); ++i){
		elements[ ++index] = i;
	}
	//We only build nonempty trees
	if( entity_handles_.size()){ 
		build_tree( elements.begin(), elements.end(), 0); 
	}
}

//Copy constructor
Element_tree( Self & s): entity_handles_( s.entity_handles_), 
			 tree_( s.tree_), moab( s.moab){}

//private functionality
private:

template< typename Entity_handles, typename Entity_map>
void construct_element_map( const Entity_handles & elements, Entity_map & map){
	typedef typename Entity_handles::const_iterator Entity_handles_iterator;
	for( Entity_handles_iterator i = elements.begin(); 
				     i != elements.end(); ++i){
		//Box box;
		//compute_bounding_box( i, box);
		//Data data( box, 0);
		//map.insert( *i, data);
	}	
}

template< typename Iterator, typename Node_index, typename Data>
bool decide_split( Iterator & begin, Iterator & end, 
		   Node_index & node, Data & data){ 
	return true;
}

template< typename Iterator, typename Node_index>
void build_tree(Iterator begin, Iterator end, Node_index node){
	_element_tree::Partition_data data;
	if ( decide_split(begin, end, node, data)){
		typedef typename Iterator::value_type Pair;
		_element_tree::Element_compare< Pair> less;
		std::sort(begin, end, less);
		
		Iterator middle_begin( begin+data.left);
		Iterator middle_end( middle_begin+data.middle);

		build_tree(begin, middle_begin, tree_[ node].left_);
		build_tree(middle_begin, middle_end, tree_[ node].middle_);
		build_tree(middle_end, end, tree_[ node].right_);
	}
}

//public functionality
public:
	

//public accessor methods
public:

//private data members  
private:
	const Entity_handles & entity_handles_;
	std::vector< Node > tree_;
	Moab & moab;

}; //class Element_tree

} // namespace moab

#endif //ELEMENT_TREE_HPP
