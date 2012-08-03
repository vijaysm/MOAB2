/** 
 * bvh_tree.hpp
 * Ryan H. Lewis
 * (C) 2012
 *
 * An element tree partitions a mesh composed of elements.
 * We subdivide the bounding box of a me, by putting boxes
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
#define BVH_TREE_DEBUG
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
		typedef typename  std::vector< std::pair< Box, Entity_handle> > 
								Entities;
		std::size_t dim;
		std::size_t child;
		double Lmax, Rmin;
		Entities entities;
		_Node& operator=( const _Node& f){
			dim = f.dim;
			child=f.child;
			Lmax=f.Lmax;
			Rmin=f.Rmin;
			entities=f.entities;
			return *this;
		}
	}; // _Node


	template< typename Split>
	class Split_comparator : 
			public std::binary_function< Split, Split, bool> {
		inline double objective( const Split & a) const{
			return a.Lmax*a.nl + a.Rmin*a.nr;
		}
		public:
		bool operator()( const Split & a, const Split & b) const{
			return  objective( a) < objective( b);
		}
	}; //Split_comparator

	template< typename Iterator>
	class Iterator_comparator : 
			public std::binary_function< Iterator, Iterator, bool> {
		public:
		bool operator()( const Iterator & a, const Iterator & b) const{
			return a->second.second < b->second.second ||
				( !(b->second.second < a->second.second) 
					&& a->first < b->first);
		}
	}; //Split_comparator


	class _Split_data {
		public:
		_Split_data(): dim( 0), nl( 0), nr( 0), split( 0.0), 
				Lmax( 0.0), Rmin( 0.0),bounding_box(), 
				left_box(), right_box(){}
       		_Split_data( const _Split_data & f): 
			dim( f.dim), nl( f.nl), nr( f.nr), 
			split( f.split), Lmax( f.Lmax), Rmin( f.Rmin),
			bounding_box( f.bounding_box),
			left_box( f.left_box), right_box( f.right_box){}
		std::size_t dim;
		std::size_t nl;
		std::size_t nr;
		double split;
		double Lmax, Rmin;
		common_tree::Box< double> bounding_box;
		common_tree::Box< double> left_box;
		common_tree::Box< double> right_box;
		_Split_data& operator=( const _Split_data & f){
			dim  	     = f.dim;
			nl   	     = f.nl; 
			nr   	     = f.nr;
        		split	     = f.split;
			Lmax 	     = f.Lmax;
			Rmin 	     = f.Rmin;
        		bounding_box = f.bounding_box;
			left_box     = f.left_box;
			right_box    = f.right_box;
			return *this;
		}
	}; //_Split_data

	class _Bucket {
		public:
		_Bucket(): size( 0), bounding_box(){}
		_Bucket( const _Bucket & f): 
		size( f.size), bounding_box(f.bounding_box){}
		_Bucket( const std::size_t size_): 
		size( size_), bounding_box(){}
		std::size_t size;
		common_tree::Box< double> bounding_box;
		_Bucket& operator=( const _Bucket & f){
			bounding_box = f.bounding_box;
			size = f.size;
			return *this;
		}
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
	typedef  common_tree::_Element_data< const _Box, double > Element_data;
	typedef typename std::tr1::unordered_map< Entity_handle, 
						  Element_data> Entity_map;
	typedef typename Entity_map::iterator Entity_map_iterator;
	typedef std::vector< Entity_map_iterator> Vector;
	//a fully balanced tree will have 2*_entities.size()
	//which is one doubling away..
	tree_.reserve( entity_handles_.size());
	Entity_map entity_map( entity_handles_.size());
	common_tree::construct_element_map( entity_handles_, 
					    entity_map, 
					    bounding_box, 
					    moab);

	for(Entity_map_iterator i = entity_map.begin(); 
				i != entity_map.end(); ++i){
		if( !box_contains_box( bounding_box, i->second.first)){
			std::cerr << "BB:" << bounding_box << "EB:" <<
			i->second.first << std::endl;
			std::exit( -1);
		}
	}
 	//_bounding_box = bounding_box;
	Vector entity_ordering;
	construct_ordering( entity_map, entity_ordering); 
	//We only build nonempty trees
	if( entity_ordering.size()){ 
	 //initially all bits are set
	 tree_.push_back( Node());
	 const int depth = build_tree( entity_ordering.begin(), 
				       entity_ordering.end(), 0, bounding_box);
	 #ifdef BVH_TREE_DEBUG
		 typedef typename Nodes::iterator Node_iterator;
		 typedef typename Node::Entities::iterator Entity_iterator;
		 std::size_t num_entities=0;
		 std::set< Entity_handle> entity_handles;
		 for(Node_iterator i = tree_.begin(); i != tree_.end(); ++i){
				num_entities += i->entities.size();
			for(Entity_iterator j = i->entities.begin(); 
					    j != i->entities.end(); ++j){
				entity_handles.insert( j->second);
			}
				
		 }
		if( num_entities != entity_handles_.size()){
			std::cout << "Entity Handle Size Mismatch!" << std::endl;
		}
		for( typename Entity_handles::iterator 
		i = entity_handles_.begin(); i != entity_handles_.end(); ++i){
			if ( entity_handles.find( *i) == entity_handles.end()){
				std::cout << "Tree is missing an entity! " << std::endl;
			}
		} 
					       
	 #endif
	 std::cout << "max tree depth: " << depth << std::endl; 
	}
}

//Copy constructor
Bvh_tree( Self & s): entity_handles_( s.entity_handles_), 
			 tree_( s.tree_), moab( s.moab), 
			 bounding_box( s.bounding_box),
			 entity_contains( s.entity_contains){}

//see FastMemoryEfficientCellLocationinUnstructuredGridsForVisualization.pdf 
//around page 9
#define NUM_SPLITS 4
#define NUM_BUCKETS 5 //NUM_SPLITS+1
#define SMAX 4
//Paper arithmetic is over-optimized.. this is safer.
template < typename Box>
std::size_t bucket_index( const Box & box, const Box & interval, 
			  const std::size_t dim) const{
	const double min = interval.min[ dim];
	const double length = (interval.max[ dim]-min)/NUM_BUCKETS;
	const double center = ((box.max[ dim] + box.min[ dim])/2.0)-min;
	#ifdef BVH_TREE_DEBUG
	#ifdef BVH_SHOW_INDEX
	std::cout << "[ " << min << " , " 
		  << interval.max[ dim] << " ]" <<std::endl;
	std::cout << "[ " 
		<< box.min[ dim] << " , " << box.max[ dim] << " ]" <<std::endl;
	std::cout << "Length of bucket" << length << std::endl;
	std::cout << "Center: " 
			<< (box.max[ dim] + box.min[ dim])/2.0 << std::endl;
	std::cout << "Distance of center from min:  " << center << std::endl;
	std::cout << "ratio: " << center/length << std::endl;
	std::cout << "index: " << std::ceil(center/length)-1 << std::endl;
	#endif
	#endif
	return std::ceil(center/length)-1;
}

template< typename Iterator, typename Bounding_box, typename Buckets>
void establish_buckets( const Iterator begin, const Iterator end, 
			const Bounding_box & interval, 
			Buckets & buckets) const{
	//put each element into its bucket
	for(Iterator i = begin; i != end; ++i){
		const Bounding_box & box = (*i)->second.first;
		for (std::size_t dim = 0; dim < NUM_DIM; ++dim){
			const std::size_t index = (dim*NUM_BUCKETS) + 
				bucket_index( box, interval, dim);
			_bvh::_Bucket & bucket = buckets[ index];
			if(bucket.size > 0){
			common_tree::update_bounding_box( bucket.bounding_box,
							  box);
			}else{ bucket.bounding_box = box; }
			bucket.size++;
		}
	}
	for( int i = 0; i < NUM_DIM*NUM_BUCKETS; ++i){
		const int N = NUM_BUCKETS;
		if( buckets[ i].size == 0){
			buckets[ i].bounding_box = interval;
			const int d = i/N;
			const int j = i%N;
			const double min = interval.min[ d];
			const double max = interval.max[ d]; 
			double left = min + (j/N)*(max-min);
			double right = min + ((j+1)/N)*(max-min);
			if ( (j+1) == NUM_BUCKETS){ right = interval.max[ d]; }
			else if ( j == 0){ left = interval.min[ d]; }

			buckets[ i].bounding_box.min[ d] = left;	
			buckets[ i].bounding_box.max[ d] = right;	
		}
	}
	#ifdef BVH_TREE_DEBUG
	Bounding_box elt_union = (*begin)->second.first;
	for(Iterator i = begin; i != end; ++i){
		const Bounding_box & box = (*i)->second.first;
		common_tree::update_bounding_box( elt_union, box);
		for (std::size_t dim = 0; dim < NUM_DIM; ++dim){
			const std::size_t index = (dim*NUM_BUCKETS) + 
				bucket_index( box, interval, dim);
			_bvh::_Bucket & bucket = buckets[ index];
			if(!box_contains_box( bucket.bounding_box, box)){
				std::cerr << "Buckets not covering elements!"
					  << std::endl;
			}
		}
	}
	if( !box_contains_box( elt_union, interval) ){
		std::cout << "element union: " << std::endl << elt_union; 
		std::cout << "intervals: " << std::endl << interval;
		std::cout << "union of elts does not contain original box!" 
			  << std::endl;
	}
	if ( !box_contains_box( interval, elt_union) ){
		std::cout << "original box does not contain union of elts" 
			  << std::endl;
		std::cout << interval << std::endl;
		std::cout << elt_union << std::endl;
	}
	Bounding_box test_box = buckets[ 0].bounding_box;
	for(int d = 0; d < NUM_DIM; ++d){
		for( int i = d*NUM_BUCKETS; i < (d+1)*NUM_BUCKETS; ++i){
			common_tree::update_bounding_box( test_box, 
					buckets[ i].bounding_box);
		}
	}
	if( !box_contains_box( test_box, interval) ){
		std::cout << "union of buckets does"
			  << " not contain original box!" 
			  << std::endl;
	}
	if ( !box_contains_box( interval, test_box) ){
		std::cout << "original box does "
			  << "not contain union of buckets" 
			  << std::endl;
		std::cout << interval << std::endl;
		std::cout << test_box << std::endl;
	}
	#endif
}

template< typename Splits, typename Buckets, typename Split_data>
void initialize_splits( Splits & splits, 
			const Buckets & buckets, 
			const Split_data & data, 
			const std::size_t total) const{
	typedef typename Buckets::value_type Bucket;
	typedef typename Buckets::const_iterator Bucket_iterator;
	typedef typename Splits::iterator Split_iterator;
	std::vector< double> length( data.bounding_box.min);
	std::transform( data.bounding_box.max.begin(),
			data.bounding_box.max.end(), 
			length.begin(), length.begin(), 
			std::minus< double>());
	for(std::size_t d = 0; d < NUM_DIM; ++d){
		Split_iterator s_begin = splits.begin()+d*NUM_SPLITS;
		Bucket_iterator b = buckets.begin()+d*NUM_BUCKETS;
		Bucket_iterator b_end = buckets.begin()+(d+1)*NUM_BUCKETS;
		for(Split_iterator s = s_begin ; s != s_begin+NUM_SPLITS; ++s,
									  ++b){
			s->left_box = b->bounding_box;
			s->right_box = (b+1)->bounding_box;
			s->nl = b->size;
			if( s != s_begin) {
				Split_iterator last = s-1; 
				s->nl += last->nl;
				common_tree::update_bounding_box( s->left_box, 
								last->left_box);
			}
			s->Lmax = s->left_box.max[ data.dim];
			s->nr = total - s->nl;
			s->split = std::distance(s_begin, s);
			s->dim = d;
		}
		for(Split_iterator s = s_begin+NUM_SPLITS-2; s >= s_begin; --s){
			common_tree::update_bounding_box( s->right_box, 
							  (s+1)->right_box);
			s->Rmin = s->right_box.min[ data.dim];
		}
	}
}

template< typename Iterator, typename Split_data>
void order_elements( const Iterator & begin, const Iterator & end, 
		     const Split_data & data) const{
	typedef typename Iterator::value_type Map_iterator;
	for(Iterator i = begin; i != end; ++i){
		const int index = bucket_index( (*i)->second.first,
						data.bounding_box, data.dim);
		(*i)->second.second = (index<=data.split)?0:1;
	}
	std::sort( begin, end, _bvh::Iterator_comparator< Map_iterator>());
}

template< typename Iterator, typename Split_data>
void median_order( const Iterator & begin, const Iterator & end, 
		      Split_data & data) const{
	typedef typename Iterator::value_type Map_iterator;
	for(Iterator i = begin; i != end; ++i){
		const double center = 
		       compute_box_center((*i)->second.first, data.dim);
		(*i)->second.second = center; 
	}
	std::sort( begin, end, _bvh::Iterator_comparator< Map_iterator>());
	const std::size_t total = std::distance( begin, end);
	Iterator middle = begin+(total/2);
	double middle_center = (*middle)->second.second;
	       middle_center += (*(++middle))->second.second;
	       middle_center /=2.0;
	data.split = middle_center;
	data.nl = std::distance( begin, middle)+1;
	data.nr = total-data.nl;
	middle++;
	data.left_box  = (*begin)->second.first;
	data.right_box = (*middle)->second.first;
	for(Iterator i = begin; i != middle; ++i){
		(*i)->second.second = 0;
		update_bounding_box( data.left_box, (*i)->second.first);
	}
	for(Iterator i = middle; i != end; ++i){
		(*i)->second.second = 1;
		update_bounding_box( data.right_box, 
				     (*i)->second.first);
	}
	data.Rmin = data.right_box.min[ data.dim];
	data.Lmax = data.left_box.max[ data.dim];
}

template< typename Iterator, typename Split_data>
void find_split(const Iterator & begin, 
		const Iterator & end, Split_data & data) const{
	typedef typename Iterator::value_type Map_iterator;
	typedef typename Map_iterator::value_type::second_type Box_data;
	typedef typename Box_data::first_type Bounding_box;
	typedef typename std::vector< Split_data> Splits;
	typedef typename Splits::iterator Split_iterator;
	typedef typename std::vector< _bvh::_Bucket> Buckets;
	Buckets buckets( NUM_BUCKETS*NUM_DIM);
	Splits splits( NUM_SPLITS*NUM_DIM, data);
	
	const Bounding_box interval = data.bounding_box;
	establish_buckets( begin, end, interval, buckets);
	const std::size_t total = std::distance( begin, end);
	initialize_splits( splits, buckets, data, total);
	Split_iterator best = std::min_element( splits.begin(), splits.end(),
			   _bvh::Split_comparator< Split_data>());
	data = *best;
	const bool use_median = (0 == data.nl) || (data.nr == 0);
	if (!use_median){ order_elements( begin, end, data); } 
	else{ median_order( begin, end, data); }
	#ifdef BVH_TREE_DEBUG
	bool seen_one=false,issue=false;
	std::size_t count_left=0, count_right=0;
	for( Iterator i = begin; i != end; ++i){
		double order = (*i)->second.second;
		if( order != 0 && order != 1){
			std::cerr << "Invalid order element !";
			std::cerr << order << std::endl;
			std::exit( -1);
		}
		if(order == 1){
			seen_one=1;
			count_right++;
			if(!box_contains_box( data.right_box, 
					     (*i)->second.first)){
				if(!issue){
				std::cerr << "Bounding right box issue!" 
					  << std::endl;
				}
				issue=true;
			}
		}
		if(order==0){
			count_left++;
			if(!box_contains_box( data.left_box, 
					     (*i)->second.first)){
				if(!issue){
				std::cerr << "Bounding left box issue!" 
					 << std::endl;
				}
				issue=true;
			}
			if(seen_one){
				std::cerr << "Invalid ordering!" << std::endl;
				std::cout << (*(i-1))->second.second 
					  << order << std::endl;
				exit( -1);
			}
		}
	}
	if( count_left != data.nl || count_right != data.nr) {
		std::cerr << "counts are off!" << std::endl;
		std::cerr << "total: " 
			  << std::distance( begin, end) << std::endl;
		std::cerr << "Dim: " << data.dim << std::endl;
		std::cerr << data.Lmax << " , " << data.Rmin << std::endl;
		std::cerr << "Right box: " << std::endl << data.right_box 
			  << "Left box: " << std::endl << data.left_box ;
		std::cerr << "supposed to be: " << 
					data.nl << " " << data.nr << std::endl;
		std::cerr << "accountant says: " << 
			count_left << " " << count_right << std::endl;
		std::exit( -1);
	}
	#endif
}

//private functionality
private:
template< typename Iterator>
int build_tree( const Iterator begin, const Iterator end, 
		const int index, const Box & box, 
		const int depth=0){
	#ifdef BVH_TREE_DEBUG
	for(Iterator i = begin; 
		     i != end; ++i){
		if( !box_contains_box( box, (*i)->second.first)){
			std::cerr << "depth: " << depth << std::endl;
			std::cerr << "BB:" << box << "EB:" <<
			(*i)->second.first << std::endl;
			std::exit( -1);
		}
	}
	#endif

	const std::size_t total_num_elements = std::distance( begin, end);
	Node & node = tree_[ index];
	//logic for splitting conditions
	if( total_num_elements > SMAX){
		_bvh::_Split_data data;
		data.bounding_box = box;
		find_split( begin, end, data);
		//assign data to node
		node.Lmax = data.Lmax; node.Rmin = data.Rmin;
		node.dim = data.dim; node.child = tree_.size();
		//insert left, right children;
		tree_.push_back( Node()); tree_.push_back( Node());
		const int left_depth=
		build_tree( begin, begin+data.nl, node.child, 
			    data.left_box, depth+1);
		const int right_depth=
		build_tree( begin+data.nl, end, node.child+1, 
			    data.right_box, depth+1);
		return std::max( left_depth, right_depth);
	}
	node.dim = 3;
	common_tree::assign_entities( node.entities, begin, end);
	return depth;
}

template< typename Vector, typename Node_index>
Entity_handle _find_point( const Vector & point, 
			   const Node_index & index) const{
	typedef typename Node::Entities::const_iterator Entity_iterator;
	const Node & node = tree_[ index];
	if( node.dim == 3){
		//check each node
		if( node.entities.size() == 0){
			std::cout << "Node with no entities?" << std::endl;
		}
		for( Entity_iterator i = node.entities.begin(); 
				     i != node.entities.end(); ++i){
			if( common_tree::box_contains_point( i->first, point) &&
				entity_contains( i->second, point)){
				return i->second;
			}
		}
		return 0;
	}
	if( point[ node.dim] <= node.Lmax){
		return _find_point( point, node.child);
	}else if( point[ node.dim] >= node.Rmin){
		return _find_point( point, node.child+1);
	}
	/* TODO:
	 * However, instead of always traversing either subtree
	 * first (e.g. left always before right), we first traverse the subtree whose
	 * bounding plane has the larger distance to the sought point. This results
	 * in less overall traversal, and the correct cell is identified more quickly.
	 */
	const Entity_handle result =  _find_point( point, node.child);
	if( result != 0){ return result; }
	return _find_point( point, node.child+1);
}

//public functionality
public:
template< typename Vector>
Entity_handle find( const Vector & point) const{
	typedef typename Vector::const_iterator Point_iterator;
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
	Parametrizer & entity_contains;

}; //class Bvh_tree

} // namespace moab

#endif //BVH_TREE_HPP
