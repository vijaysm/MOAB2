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
#include <deque>
#include <bitset>

#ifndef ELEMENT_TREE_HPP
#define ELEMENT_TREE_HPP
namespace moab {
//forward declarations

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

struct Box{
	typedef typename std::vector< double> Vector;
	Box(): max(3,0.0), min(3,0.0){}
	Box( const Box & from): max( from.max), min( from.min){}
	Box& operator=( const Box & from){
		max = from.max;
		min = from.min;
		return *this;
	}
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
  //max - min for correctness. it *might* be slightly slower
  //but that is unlikely. compiler should optimize this.
  const std::size_t split_objective( const Data & a) const {
	const int max = 2*(a.sizes[ 2]>a.sizes[ 0]);
	return (a.sizes[ max] -a.sizes[ 2*(1-(max==2))]) + 
					a.sizes[ 1]*a.sizes[ 1];
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
	Node& operator=( const Node & from){
		left_=from.left_;
		middle_=from.middle_;
		right_=from.right_;
		dim=from.dim;
		split=from.split;
		left_line=from.left_line;
		right_line=from.right_line;
		entity_handles_=from.entity_handles_;
		return *this;
	}	
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
	//default constructor
	Partition_data(){}
	Partition_data( const Partition_data & f){
		*this=f;
	}
	Partition_data& operator=( const Partition_data & f){
		for (int i = 0; i < 3; ++i){sizes[ i] = f.sizes[ i];}
		bounding_box = f.bounding_box;
		split = f.split;
		left_line = f.left_line;
		right_line = f.right_line;
		return *this;
	}
	std::size_t sizes[ 3];
	Box bounding_box;
	double split;
	double left_line;
	double right_line;
	std::size_t& left()   { return sizes[ 0]; }
	std::size_t& middle() { return sizes[ 1]; }
	std::size_t& right()  { return sizes[ 2]; }
}; // Partition_data

template< typename Iterator>
struct Element_compare: public std::binary_function< Iterator, Iterator, bool>{
	Element_compare( const int best_split) loc( best_split){}
	//TODO: finish off the mask
	bool operator()( const Iterator a, const Iterator b) const{
		return a->second.second.to_ulong() < 
			b->second.second.to_ulong();
	}
	const int loc;
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
	//int is because we only need to store 	
	#define MAX_ITERATIONS 2
	typedef typename std::pair< Box, std::bitset<3*2*MAX_ITERATIONS> > 
								Element_data;
	typedef typename std::vector< Node> Nodes;
	typedef typename Entity_handles::value_type Entity_handle;
	//TODO: we really want an unordered map here..
	typedef typename std::map< Entity_handle, Element_data> Element_map;
	typedef typename std::vector< typename Element_map::iterator> 
								Element_list;
//public methods
public:
//Constructor
Element_tree( Entity_handles & _entities, Moab & _moab): 
	entity_handles_( _entities), tree_(), moab( _moab) {
	tree_.reserve( _entities.size());
	Element_map element_map;
	_element_tree::Partition_data _data;
	construct_element_map( entity_handles_, element_map, 
					_data.bounding_box);
	Element_list element_ordering( element_map.size());
	std::size_t index = 0;
	for(typename Element_map::iterator i = element_map.begin(); 
					  i != element_map.end(); ++i, ++index){
		element_ordering[ index] = i;
	}
	//We only build nonempty trees
	if( element_ordering.size()){ 
		//initially all bits are set
		std::bitset< 3> directions( 7);
		tree_.push_back( Node());
		build_tree( element_ordering.begin(), element_ordering.end(), 0,
			    directions, _data); 
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

template< typename Entity_handles, typename Element_map, typename Bounding_box>
void construct_element_map( const Entity_handles & elements, 
			    Element_map & map, 
			    Bounding_box & bounding_box){
	typedef typename Element_map::mapped_type Box_data;
	typedef typename Entity_handles::value_type Entity_handle;
	typedef typename Box_data::first_type Box;
	typedef typename Box_data::second_type Bitset;
	typedef typename Entity_handles::iterator Entity_handles_iterator;
	typedef typename std::vector< double> Coordinate;
	for( Entity_handles_iterator i = elements.begin(); 
				     i != elements.end(); ++i){	
		Coordinate coordinate(3, 0.0);
		const Entity_handle* vertex_handle;
		int num_vertices=0;
		moab.get_connectivity( *i, vertex_handle, num_vertices, true);
		Box box;
		for( int j = 0; j < num_vertices; ++j){
			moab.get_coords( vertex_handle+j, 1, &coordinate[ 0]);
			update_bounding_max( box.max, coordinate);
			update_bounding_min( box.min, coordinate);
		}
		update_bounding_max( bounding_box.max, box.max);
		update_bounding_min( bounding_box.min, box.min);
		map.insert( make_pair( *i, Box_data( box, Bitset( 0))));
	}
}

template< typename Iterator, typename Partition_data_iterator>
void find_optimal_split( Iterator & begin, Iterator & end, 
			 Partition_data_iterator & data_iterator, 
			 const std::size_t dim){
	typedef typename Partition_data_iterator::value_type Partition_data;
	typedef typename Iterator::value_type::value_type Map_value_type;
	typedef typename Map_value_type::second_type::second_type Bitset;
	//get middle line in bounding box of this dimension
	const Box & box( data_iterator->bounding_box);
	int iteration = 0;
	do{
	Partition_data & data = *data_iterator;
	double & split = data.split =  (box.max[ dim] - box.min[ dim])/2.0;
	double & left_line = data.left_line = split;
	double & right_line = data.right_line = split;
	//for each elt determine if left/middle/right
	for(Iterator i = begin; i != end; ++i){
		const Box & _box =  (*i)->second.first;
		Bitset & bits =  (*i)->second.second;
		//will be 0 if on left, will be 1 if in the middle
		//and 2 if on the right;
		const int side = !(split < _box.min[ dim]) + 
					(split > _box.max[ dim]);
		//keep track of partition data,
		data.sizes[ side]++;
		//if we are in the middle update the middle bounding lines
		left_line = (side==1)*std::min(left_line,_box.min[ dim]);
		right_line = (side==1)*std::max(right_line,_box.max[ dim]);
		//set the corresponding bits in the bit vector
		// looks like: [x_1 = 00 | x_2 = 00 | .. | z_1 = 00 | z_2 = 00]
		// two bits, but we encode 2 as 11 _not_ 10
		const int index = 4*dim + 2*iteration;
		bits.set( index, side>0);
		bits.set( index+1, side==2);
	}
	//if balanced exit.
	const int max = 2*(data.sizes[ 2]>data.sizes[ 0]);
	const int min = 2*(1-(max==2));
	const std::size_t total = std::distance(begin,end);
	if (data.sizes [ max] - data.sizes[ min] < .05*total){
		return;	
	}
	const int sign = max-1;
	split += sign*((data.sizes[ max] - data.sizes[ min])/
					data.sizes[ max])*(split/2.0);
	} while( ++iteration < MAX_ITERATIONS);
}

template< typename Iterator, 
	  typename Partition_data, 
	  typename Directions>
int determine_split( Iterator & begin, 
		      Iterator & end, 
		      Partition_data & data, 
		      const Directions & directions){ 
	typedef typename Iterator::value_type Pair;
	typedef typename std::vector< Partition_data> Vector;
	typedef typename _element_tree::Split_comparator< Partition_data> 
								Comparator;
	typedef typename Vector::iterator Split;
	Vector splits( 2*directions.count(), data);
	Split s( splits.begin());
	for (std::size_t dir = 0; dir < directions.size(); ++dir){
		if( directions.test( dir)){
			find_optimal_split( begin, end, s, dir);
			s+=2;
		}
	}
	Split best = std::min_element( splits.begin(), 
				       splits.end(), Comparator());
	data = *best;
	return std::distance( splits.begin(), best);
}

template< typename Node_index, typename Partition_data>
void assign_data_to_node( const Node_index & node, const Partition_data & data){
	tree_[ node].split = data.split;
	tree_[ node].left_line = data.left_line;  
	tree_[ node].right_line = data.right_line;
}
template< typename Node_index>
void extend_tree( const Node_index node){
	Nodes children( 3);
	std::size_t _index = tree_.size();
	tree_.insert( tree_.end(), children.begin(), children.end());
	tree_[ node].left_ = _index++;
	tree_[ node].middle_ = _index++;
	tree_[ node].right_ = _index;
}

//define here for now.
#define ELEMENTS_PER_LEAF 1000
template< typename Iterator, typename Node_index, 
	  typename Directions, typename Partition_data>
void build_tree( Iterator begin, Iterator end, 
		 const Node_index node, 
		 const Directions & directions, 
		 Partition_data & _data,
		 const int depth = 0, 
		 const bool is_middle = false){
	typedef typename Iterator::value_type Pair;
	std::size_t number_elements = std::distance(begin, end);
	if ( number_elements > ELEMENTS_PER_LEAF && 
		(!is_middle || directions.any())){
		int best_split = determine_split( begin, end, 
						 _data, directions); 
		assign_data_to_node( node, _data);
		_element_tree::Element_compare< Pair> less( best_split);
		std::sort( begin, end, less);
		Iterator middle_begin( begin+_data.left());
		Iterator middle_end( middle_begin+_data.middle());

		//append elements after current node
		std::cout << tree_.size() << std::endl;
		extend_tree( node);
		std::cout << "tree extended" << std::endl;
		std::cout << tree_.size() << std::endl;
		
		_element_tree::Partition_data data( _data);
		//left subtree
		build_tree( begin, middle_begin, tree_[ node].left_, 
			    directions, data, depth+1,  is_middle);
		//right subtree
		build_tree( middle_end, end, tree_[ node].right_,
			    directions, data, depth+1, is_middle);

		//force the middle subtree to split
		//in a different direction from this one
		//middle subtree
		Directions new_direction( directions);
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
