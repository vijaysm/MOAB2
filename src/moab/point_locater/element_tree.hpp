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
#include <numeric>
#include <cmath>

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
	template< typename Iterator>
	Box( Iterator begin, Iterator end):max( begin,end), min(begin,end){}
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
  //we minimizes ||left| - |right|| + |middle|^2 
  const std::size_t split_objective( const Data & a) const {
	const int max = 2*(a.second.sizes[ 2]>a.second.sizes[ 0]);
	return (a.second.sizes[ max] -a.second.sizes[ 2*(1-(max==2))]) + 
					a.second.sizes[ 1]*a.second.sizes[ 1];
  }
  bool operator()( const Data & a, const Data & b) const {
  	return split_objective( a) < split_objective( b);
  }
};

template< typename _Entity_handles, typename _Entities>
class Node{
	//public types:
	public:
	typedef _Entity_handles Entity_handles;
	typedef _Entities Entities;

	//private types:
	private:
	typedef Node< _Entity_handles, _Entities> Self;

	//Constructors
	public:
	//Default constructor
	Node(): left_( -1), middle_( -1), right_( -1),
	        dim( -1), split( 0),
		 left_line( 0), right_line( 0),
		entities( NULL) {}

	//Destructor
	~Node(){ delete entities; }

	//Copy constructor
	Node( const Self & from): 
	left_( from.left_), middle_( from.middle_), right_( from.right_),
	dim( from.dim), split( from.split),
	left_line( from.left_line), right_line( from.right_line), 
	entities( from.entities) {}

	public:
	template< typename Iterator>
	void assign_entities(const Iterator & begin, const Iterator & end){
		std::cout << "b4 segfault ? " << std::endl;
		entities = new Entities();
		std::cout << "segfault ? " << std::endl;
		Entities & _entities = *entities;
		_entities.reserve( std::distance( begin, end)); 
		std::cout << "segfault ? " << std::endl;
		for( Iterator i = begin; i != end; ++i){
			_entities.push_back( std::make_pair((*i)->second.first, 
							    (*i)->first));
		}
		std::cout << "segfault ? " << std::endl;
	}

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
		entities=from.entities;
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
	Entities * entities;

	//Element_tree can touch my privates.
	template< Entity_handles, typename B> friend class moab::Element_tree;
}; //class Node

struct Partition_data{
	//default constructor
	Partition_data():sizes(3,0){}
	Partition_data( const Partition_data & f){
		*this=f;
	}
	Partition_data( const Box & _box, int _dim): sizes(3,0),
	bounding_box( _box), split((_box.max[ _dim] + _box.min[ _dim])/2.0), 
	left_line( split), right_line( split), dim( _dim){}
	Partition_data& operator=( const Partition_data & f){
		sizes = f.sizes;
		bounding_box = f.bounding_box;
		split = f.split;
		left_line = f.left_line;
		right_line = f.right_line;
		return *this;
	}
	std::vector< std::size_t> sizes;
	Box bounding_box;
	double split;
	double left_line;
	double right_line;
	int dim;
	std::size_t& left()   { return sizes[ 0]; }
	std::size_t& middle() { return sizes[ 1]; }
	std::size_t& right()  { return sizes[ 2]; }
}; // Partition_data

template< typename Iterator>
struct Element_compare: public std::binary_function< Iterator, Iterator, bool>{
	typedef typename Iterator::value_type Pair;
	typedef typename Pair::second_type::second_type Bitset;
	Element_compare( const Bitset & _mask): mask( _mask){}
	bool operator()( const Iterator a, const Iterator b) const{
		return ((a->second.second)&mask).to_ulong() < 
			((b->second.second)&mask).to_ulong();
	}
	const Bitset mask;
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
	typedef typename std::pair< Box, typename Entity_handles::value_type> 
							Leaf_element;
	typedef typename _element_tree::Node< Entity_handles,
					      std::vector< Leaf_element> > Node;
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
		int depth = 0;
		build_tree( element_ordering.begin(), 
			    element_ordering.end(),
			    0, directions, _data, depth); 
	}
}
/*
439 void build_tree( Iterator begin, Iterator end,
440                  const Node_index node,
441                  const Directions & directions,
442                  Partition_data & _data,
443                  int & depth,
444                  const bool is_middle = false){
*/

//Copy constructor
Element_tree( Self & s): entity_handles_( s.entity_handles_), 
			 tree_( s.tree_), moab( s.moab){}

//private functionality
private:

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
	typedef typename Coordinate::iterator Coordinate_iterator;
	for( Entity_handles_iterator i = elements.begin(); 
				     i != elements.end(); ++i){	
		const Entity_handle* vertex_handle;
		int num_vertices=0;
		moab.get_connectivity( *i, vertex_handle, num_vertices);
		Coordinate coordinate(3*num_vertices, 0.0);
		moab.get_coords( vertex_handle, num_vertices, &coordinate[ 0]);
		Box box( coordinate.begin(), coordinate.begin()+3);
		bounding_box = box;
		for( Coordinate_iterator j = coordinate.begin()+3; 
				         j != coordinate.end(); j+=3){
			update_bounding_max( box.max, j);
			update_bounding_min( box.min, j);
		}
		update_bounding_max( bounding_box.max, box.max.begin());
		update_bounding_min( bounding_box.min, box.min.begin());
		map.insert( make_pair( *i, Box_data( box, Bitset( 0))));
	}
}

template< typename Iterator, typename Split_data>
void compute_split( Iterator & begin, Iterator & end, 
			 Split_data & split_data, bool iteration=false){
	typedef typename Iterator::value_type::value_type Map_value_type;
	typedef typename Map_value_type::second_type::second_type Bitset;
	//we will update the left/right line
	double & left_line = split_data.left_line;
	double & right_line = split_data.right_line;
	double & split = split_data.split;
	const int & dim = split_data.dim;
	#ifdef ELEMENT_TREE_DEBUG
 	std::cout << std::endl; 
	std::cout << "-------------------" << std::endl; 
	std::cout << "dim: " << dim << " split: " << split << std::endl;
	std::cout << "bounding_box min: "; 
	print_vector( split_data.bounding_box.min); 
	std::cout << "bounding_box max: "; 
	print_vector( split_data.bounding_box.max);
	#endif
	//for each elt determine if left/middle/right
	for(Iterator i = begin; i != end; ++i){
		const Box & box =  (*i)->second.first;
		Bitset & bits =  (*i)->second.second;
		//will be 0 if on left, will be 1 if in the middle
		//and 2 if on the right;
		const bool on_left = (box.max[ dim] < split);
		const bool on_right = (box.min[ dim] > split);
		const bool in_middle = !on_left && !on_right;
		//TODO: remove branch
		int side;
		if(in_middle){
			side=1;
			left_line  = std::min( left_line,  box.min[ dim]);
			right_line = std::max( right_line, box.max[ dim]);
		}else if( on_left){
			side=0;
		}else if( on_right){
			side=2;
		}
		//keep track of partition data,
		split_data.sizes[ side]++;
		//set the corresponding bits in the bit vector
		// looks like: [x_1 = 00 | x_2 = 00 | .. | z_1 = 00 | z_2 = 00]
		// two bits, but we encode 2 as 11 _not_ 10
		const int index = 4*dim + 2*iteration;
		bits.set( index, side>0);
		bits.set( index+1, side==2);
	}	
	#ifdef ELEMENT_TREE_DEBUG
	std::cout <<  " left_line: " << left_line;
	std::cout <<  " right_line: " << right_line << std::endl;
	std::cout << "computed partition size: ";
	print_vector( split_data.sizes);
	std::cout << "-------------------" << std::endl; 
	#endif
}

template< typename Split_data>
bool update_split_line( Split_data & data) const{
	std::size_t total = std::accumulate( data.sizes.begin(), 
					     data.sizes.end(), 0);
	const int max = 2*(data.sizes[ 2]>data.sizes[ 0]);
	const int min = 2*(1-(max==2));
	bool all_in_middle = data.sizes[ max]==0 && data.sizes[ min]==0;
	double balance_ratio = data.sizes[ max] - data.sizes[ min];
	if ( !all_in_middle && balance_ratio < .05*total){ return false; }
	if( !all_in_middle){
		//if we have some imbalance on left/right 
		//try to fix the situation 
		balance_ratio /= data.sizes[ max];
		data.split += (max-1)*balance_ratio*(data.split/2.0);
	}else{
		//if everything is in the middle, wiggle a bit in
		//the larger direction
		double left_distance = std::abs(data.left_line-data.split);
		double right_distance = std::abs(data.right_line-data.split);
		data.split += (left_distance>right_distance)?
				(-left_distance/2.0):(right_distance/2.0);
	}
	data.left_line = data.right_line = data.split;
	data.sizes.assign( data.sizes.size(), 0);
	return true;
}

template< typename Iterator, 
	  typename Split_data, 
	  typename Directions>
std::size_t determine_split( Iterator & begin, 
		      Iterator & end, 
		      Split_data & data, 
		      const Directions & directions){ 
	typedef typename Iterator::value_type Pair;
	typedef typename std::map< std::size_t, Split_data> Splits;
	typedef typename Splits::value_type Split;	
	typedef typename _element_tree::Split_comparator< Split> Comparator;
	Splits splits;
	for (std::size_t dir = 0; dir < directions.size(); ++dir){
		if( directions.test( dir)){
			Split_data split_data( data.bounding_box, dir);
			compute_split( begin, end, split_data);
			splits.insert( std::make_pair(2*dir, split_data));
			if( update_split_line( split_data)){
				compute_split( begin, end, split_data, true);
				splits.insert( std::make_pair( 2*dir+1,
							       split_data) );
			}
		}
	}
	Split best = *std::min_element( splits.begin(), splits.end(), 
					Comparator());
	data = best.second;
	return best.first;
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
#define ELEMENTS_PER_LEAF 10
template< typename Iterator, typename Node_index, 
	  typename Directions, typename Partition_data>
void build_tree( Iterator begin, Iterator end, 
		 const Node_index node, 
		 const Directions & directions, 
		 Partition_data & _data,
		 int & depth, 
		 const bool is_middle = false){
	typedef typename Iterator::value_type Map_iterator;
	typedef typename Map_iterator::value_type::second_type Element_data;
	typedef typename Element_data::second_type Bitset;
	std::size_t number_elements = std::distance(begin, end);
	if ( number_elements > ELEMENTS_PER_LEAF && 
		(!is_middle || directions.any())){
		std::size_t best_split = determine_split( begin, end, 
						 _data, directions); 
		Bitset mask( 0);
		mask.flip( best_split).flip(best_split+1);
		_element_tree::Element_compare< Map_iterator> less( mask);
		std::sort( begin, end, less);
		//update the tree
		assign_data_to_node( node, _data);
		extend_tree( node);
		Iterator middle_begin( begin+_data.left());
		Iterator middle_end( middle_begin+_data.middle());
		_element_tree::Partition_data data( _data);
		if( is_middle){
			std::cout << "middle depth: " << depth << std::endl;
		}	
		//left subtree
		std::vector< int> depths(3, depth+1);
		build_tree( begin, middle_begin, tree_[ node].left_, 
			    directions, data, depths[ 0],  is_middle);
		//right subtree
		build_tree( middle_end, end, tree_[ node].right_,
			    directions, data, depths[ 1], is_middle);
		//force the middle subtree to split
		//in a different direction from this one
		//middle subtree
		Directions new_direction( directions);
		new_direction.flip( tree_[node].dim);
		build_tree( middle_begin, middle_end, tree_[ node].middle_,
			    new_direction, data, depths[ 2], true);
		depth = *std::max_element(depths.begin(), depths.end());
	}else{
		std::cout << "here?" << std::endl;
		tree_[ node].assign_entities( begin, end);
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
