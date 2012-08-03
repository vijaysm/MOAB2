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


struct Box{
	typedef typename std::vector< double> Vector;
	Box(): max(3,0.0), min(3,0.0){}
	Box( const Vector & max_, const Vector & min_): max(max_), min(min_){}
	Box( const Box & from): max( from.max), min( from.min){}
	Vector max;
	Vector min;
}; //Box

template< typename _Entity_handles, 
	  typename _Boxes, 
	  typename _Moab> class Element_tree;

//non-exported functionality
namespace _element_tree {

template< typename Data>
struct Split_comparator {
  const std::size_t split_objective( const Data & a) const {
  	return std::abs( a.left - a.right) + a.middle*a.middle;
  }
  
  bool operator()( const Data & a, const Data & b) const {
  	return split_objective( a) < split_objective( b);
  }
};

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
	        dim( -1), split( 0),
		 left_line( 0), right_line( 0),
		entity_handles_( NULL) {}

	//Destructor
	~Node(){ delete entity_handles_; }

	//Constructor (leaf)
	Node( Entity_handles & _entities, Box & box): 
		left_( -1), middle_( -1), right_( -1),
	        dim( -1), split( 0), 
		left_line( 0), right_line( 0), 
		entity_handles_( &_entities) {}

	//Copy constructor
	Node( const Self & from): 
	left_( from.left_), middle_( from.middle_), right_( from.right_),
	dim( from.dim), split( from.split),
	left_line( from.left_line), right_line( from.right_line), 
	entity_handles_( from.entity_handles_) {}

	// Functionality
	public: 
	bool leaf() const { return left_ == -1 && 
			    middle_ == -1 && 
			    right_ == -1; }
	
	//private data members:
	private:
	//indices of children
	int  left_;
	int  middle_;
	int  right_;
	int dim; //split dimension
	double split; //split position
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
	Box bounding_box;
}; // Partition_data

template< typename Iterator>
struct Element_compare: public std::binary_function< Iterator, Iterator, bool>{
	bool operator()( const Iterator a, const Iterator b){
		return (*a).first < (*b).first;
	}
}; // Element_compare



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
	typedef typename std::pair< typename Boxes::value_type, 
				    std::size_t> Element_data;
	typedef typename std::deque< Node> Nodes;
	typedef typename Entity_handles::value_type Entity_handle;
	//TODO: we really want an unordered map here..
	typedef typename std::map< Entity_handle, Element_data> Entity_map;
	typedef typename Entity_map::const_iterator Entity_map_iterator;
	typedef typename std::vector< Entity_map_iterator> Element_list;
//public methods
public:
//Constructor
Element_tree( Entity_handles & _entities, Moab & _moab): 
	entity_handles_( _entities), tree_(), moab( _moab) {
	tree_.reserve( _entities.size());
	Entity_map map;
	Box bounding_box;
	construct_element_map( entity_handles_, map, box);
	Element_list elements( map.size());
	std::size_t index = 0;
	for(typename Entity_map::const_iterator i = map.begin(); 
					 i != map.end(); ++i){
		elements[ ++index] = i;
	}
	//We only build nonempty trees
	if( elements.size()){ 
		//initially all bits are set
		std::bitset< 3> directions( 7);
		build_tree( elements.begin(), elements.end(), 0, directions); 
	}
}

//Copy constructor
Element_tree( Self & s): entity_handles_( s.entity_handles_), 
			 tree_( s.tree_), moab( s.moab){}

//private functionality
private:

template<typename Coordinate>
void update_bounding_max( Coordinate & max, const Coordinate & coordinate){
	typedef typename Coordinate::iterator Iterator;
	typedef typename Coordinate::const_iterator Const_iterator;
	Const_iterator j = coordinate.begin();
	for( Iterator i = max.begin(); i != max.end(); ++i, ++j){
		*i = std::max( *i, *j);
	}
}

template<typename Coordinate>
void update_bounding_min( Coordinate & min, const Coordinate & coordinate){
	typedef typename Coordinate::iterator Iterator;
	typedef typename Coordinate::const_iterator Const_iterator;
	Const_iterator j = coordinate.begin();
	for( Iterator i = min.begin(); i != min.end(); ++i, ++j){
		*i = std::min( *i, *j);
	}
}

template< typename Entity_handles, typename Entity_map, typename Bounding_box>
void construct_element_map( const Entity_handles & elements, 
			    Entity_map & map, Bounding_box & box){
	typedef typename Entity_map::mapped_type Box_data;
	typedef typename Entity_handles::value_type Entity_handle;
	typedef typename Box_data::first_type Box;
	typedef typename Entity_handles::const_iterator Entity_handles_iterator;
	typedef typename std::vector< double> Coordinate;
	Coordinate coordinate(3, 0.0);
	Entity_handle vertex_handle;
	std::size_t num_vertices=0;
	for( Entity_handles_iterator i = elements.begin(); 
				     i != elements.end(); ++i){	
		//moab.get_connectivity( &*i, vertex_handle, num_vertices, true);
		Box box;
		for( std::size_t j = 0; j < num_vertices; ++j){
			//moab.get_coords( vertex_handle+j, 1, &coordinate[ 0]);
			update_bounding_max( box.max, coordinate);
			update_bounding_min( box.min, coordinate);
		}
		update_bounding_max( bounding_box.max, box.max);
		update_bounding_min( bounding_box.min, box.min);
		map.insert( make_pair( *i, Box_data( box, 0)));
	}
}

template< typename Iterator, typename Partition_data>
void find_optimal_split( Iterator & begin, Iterator & end, 
			 Partition_data & data, const std::size_t dim){
	//get middle line in bounding box of this dimension
	//for each elt determine if left/middle/right
	//store in element split data
	//keep track of partition data,
	//if balanced exit.
	//else move in exactly one direction.
}

template< typename Iterator, 
	  typename Node_index, 
	  typename Partition_data, 
	  typename Directions>
void determine_split( Iterator & begin, 
		      Iterator & end, 
		      const Node_index node, 
		      Partition_data & data, 
		      const Directions & directions){ 
	typedef typename Iterator::value_type Pair;
	typedef typename std::vector< Partition_data> Vector;
	typedef typename _element_tree::Split_comparator< Partition_data> 
								Comparator;
	typedef typename Vector::iterator Split;
	Vector splits( directions.count());
	Split s( splits.begin());
	for (std::size_t dir = 0; dir < directions.size(); ++dir){
		if( directions.test( dir)){
			find_optimal_split( begin, end, *s, dir);
			++s;
		}
	}
	data = *std::min_element( splits.begin(), splits.end(), Comparator());
}

//define here for now.
#define ELEMENTS_PER_LEAF 1000
template< typename Iterator, typename Node_index, 
	  typename Directions, typename Partition_data>
void build_tree( Iterator begin, Iterator end, 
		 const Node_index node, 
		 const Directions & directions, 
		 const Partition_data & _data,
		 const int depth = 0, 
		 const bool is_middle = false){
	std::size_t number_elements = std::distance(begin, end);
	if ( number_elements > ELEMENTS_PER_LEAF && 
		(!is_middle || directions.any())){
		_element_tree::Partition_data data( _data);
		determine_split( begin, end, node, data, directions); 
		typedef typename Iterator::value_type Pair;
		_element_tree::Element_compare< Pair> less;

		std::sort(begin, end, less);
		
		Iterator middle_begin( begin+data.left);
		Iterator middle_end( middle_begin+data.middle);

		//append after the current node these elements	
		tree_.insert( tree_.begin()+node+1, 3, tree_[ node]);

		//left subtree
		build_tree( begin, middle_begin, tree_[ node].left_, 
			    directions, data, depth+1,  is_middle);
		//right subtree
		build_tree( middle_end, end, tree_[ node].right_,
			    directions, data, depth+1, is_middle);

		//force the middle subtree to split
		//in a different direction from this one
		//middle subtree
		Directions new_direction( direction);
		new_direction.flip( tree_[node].dim);
		build_tree( middle_begin, middle_end, tree_[ node].middle_,
			    new_direction, data, depth+1, true);
	}else{
		//we are a leaf node
	}
}

//public functionality
public:
	

//public accessor methods
public:

//private data members  
private:
	const Entity_handles & entity_handles_;
	Nodes tree_;
	Moab & moab;

}; //class Element_tree

} // namespace moab

#endif //ELEMENT_TREE_HPP
