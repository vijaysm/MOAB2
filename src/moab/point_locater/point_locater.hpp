/** 
 * point_locater.hpp
 * Ryan H. Lewis 
 * Copyright 2012
 */
#include <vector>
#ifndef POINT_LOCATER_HPP
#define POINT_LOCATER_HPP

namespace moab {

template< typename _Tree,
	  typename _Boxes>
class Point_search {

//public types
public:
	typedef  _Tree Tree;
	typedef  _Boxes Boxes;
//	typedef typename Tree::Elements::value_type Element;
	typedef int Error;
	
//private types
private: 
	typedef Point_search< _Tree, 
			      _Boxes> Self; 
//public methods
public:

//Constructor
Point_search( Tree & _tree,
	      Boxes & _boxes): 
	      tree_( _tree),
	      boxes( _boxes){}

//Copy constructor
Point_search( Self & s): tree_( s.tree_), 
			 boxes( s.boxes){}

//private functionality
private:

template< typename Point_map, typename List>
void resolve_boxes( const Point_map & query_points,  List & list){
       /*
	typedef typename std::vector< bool> Bitmask;
	typedef typename Points::const_iterator Point;
	typedef typename Tree::const_element_iterator Element;
	typedef typename std::vector< std::size_t> Processor_list;
	typedef typename List::value_type Tuple;
	const Element & end = tree_.element_end();
	Bitmask located_points( query_points.size(), 0);
	std::size_t index=0;
	for( Point i = query_points.begin(); i != query_points.end(); ++i,++index){
		const Element  & element = tree_.find( *i);
		if(element != end){
			located_points[ index] = 1;
		}
	}
	for(int i = 0; i < located_points.size(); ++i){
		if(!located_points[ i]){
			Processor_list processors;
			const Point & point = query_point.begin()+i;
			resolve_box_for_point( point, processors);
			for( std::size_t p = processors.begin(); 
					 p != processors.end(); ++p){
				list.push_back( Tuple( *point, *p) );
			}
		}
	}
	*/
}


//public functionality
public:
template< typename Point_map, typename Entities, typename Communicator> 
Error locate_points( Point_map & query_points, Entities & entities, Communicator & comm){
	/*
	//temporary types
	typedef typename Point_map::key_type Tuple;
	typedef typename std::vector< Tuple> List;
	List & scatter_points; 
	resolve_boxes( query_points, scatter_points);
	 //Commented out for now
	//scatter-gather
	//transfer( scatter_points);
	//find_local_points( scatter_points);
	//send back to target
	*/
	return 0;
} 

template< typename Points, typename Entities> 
Error locate_points( Points & query_points, Entities & entities) const{
	Entities result;
	result.reserve( query_points.size());	
	typedef typename Points::const_iterator Point_iterator;
	for(Point_iterator i = query_points.begin(); 
			   i != query_points.end(); ++i){
			result.push_back( tree_.find( *i));
	}
	return 0;
} 

//public accessor methods
public:
Tree &		tree() 		const { return tree_; }

//private data members  
private:
const Tree & tree_;
const Boxes & boxes;
}; //class Point_search

} // namespace moab

#endif //POINT_LOCATER_HPP
