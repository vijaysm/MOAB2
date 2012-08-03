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

#include "iBase.h"
#include "moab/Range.hpp"
#include "moab/Interface.hpp"
#include "moab/CartVect.hpp"

#include "moab/ParallelComm.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/GeomUtil.hpp"
#include "ElemUtil.hpp"
#include "moab/CN.hpp"
#include "moab/Types.hpp"
//#include "iMesh_extensions.h"

#include <stdio.h>
#include <iostream>

#include "moab/gs.hpp"
#include "moab/TupleList.hpp"

#include "assert.h"

namespace moab {

namespace _point_search {
template< typename Tuple_list, typename Vector> 
void write_data( Tuple_list & target_points, 
		 const Tuple_list & source_points, 
		 Vector & located_points){
	// count non-negatives
	// TODO: replace with int num_pts = std::count()
	int num_pts = 0;
	for (unsigned int i = 0; i < source_points.get_n(); i++){
	  if (-1 != source_points.vi_rd[3*i+2]) num_pts++;  
	}
	target_points.initialize(3, 0, 0, 1, num_pts);
	target_points.enableWriteAccess();
	for (unsigned int i = 0; i < source_points.get_n(); i++) {
	  //Don't send message for discarded points
	  if (source_points.vi_rd[3*i+2] != -1) { 
	    const int locIndex = source_points.vi_rd[3*i+1];
	    //asked 2+ procs if they have point p, they both said 
	    //yes, we'll keep the one with lowest rank
	    if(located_points[locIndex]){  
	      //TODO: check that the cases where both processors 
	      //say yes is justified 
	      continue;
	    }
	    located_points[locIndex] = 1;
	    target_points.vi_wr[3*target_points.get_n()]     = source_points.vi_rd[3*i];
	    target_points.vi_wr[3*target_points.get_n() + 1] = source_points.vi_rd[3*i+1];
	    target_points.vi_wr[3*target_points.get_n() + 2] = source_points.vi_rd[3*i+2];
	    target_points.inc_n();
	  }
	}
	target_points.disableWriteAccess();
}


//hide the ugliness
template< typename Boxes, typename Point>
bool in_box( const Boxes & boxes, const Point & point, 
	     const std::size_t i, const std::size_t j ){
	//TODO: genericify
	const double x = *point; 
	const double y = *(++point); 
	const double z = *(++point);
	//TODO: again, make not ugly
	return   box[j]   <= x && x <= boxes[6*j+3] && 
		 boxes[j+1] <= y && y <= boxes[6*j+4] && 
		 boxes[j+2] <= z && z <= boxes[6*j+5];
}

template< typename Tuple_list, typename Points>
void set_target_point( Tuple_list & target_points, 
		       const Point & query_points, 
		       const std::size_t i, const std::size_t j){
	  // if in this proc's box, will send to proc to test further
	  // check size, grow if we're at max
	  if (target_points.get_n() == target_points.get_max()){ 
	    target_points.resize( (3*target_points.get_max() + 1)/2);
	  }
	  const std::size_t size = target_points.get_n(); 
	  target_points.vi_wr[2*size] = j;
	  target_points.vi_wr[2*size+1] = i/3;
	  target_points.vr_wr[3*size] = query_points[i];
	  target_points.vr_wr[3*size+1] = query_points[i+1];
	  target_points.vr_wr[3*size+2] = query_points[i+2];
	  target_points.inc_n();
}

template< typename Points, typename Boxes, typename Tuple_list>
void determine_box( const Points & query_points,
	          const Boxes & boxes,
		  Tuple_list & target_points,
		  const std::size_t i, 
		  const std::size_t begin,
		  const std::size_t end){
	for (std::size_t j = begin; j < end; j++) {
	    // test if point is in proc's box
	  if( in_box( boxes, query_points, i,j) ) {
	      set_target_point( target_points, query_points, i, j);
	  }
	}
}
 
template< typename Points, typename Bitmask, 
	  typename Boxes, typename Tuple_list>
ErrorCode distribute_data( const Points & query_points,
			   Tuple_list & result,
			   Bitmask & located_points, 
			   const Boxes & boxes,
			   const std::size_t rank, 
			   const std::size_t size, 
			   const std::size_t num_points,
			   const double rel_eps, 
			   const double abs_eps){ 
  	
	typedef typename Points::const_iterator Point_iterator; 
	typedef typename Bitmask::const_iterator Bit_iterator; 
	bool point_located = false;
	for (Point_iterator p = query_points.begin(); 
			    p != query_points.end(); 
			    p+=3) {
	    //TODO: use stride iterators to get rid of this stupidity.
	    const std::size_t j = std::distance(query_points.begin(), p)/3;
	    // determine if local box has TODO
	    // TODO: Check that test_local_box does not assume point_located
	    // is initialized
	    test_local_box( p, result, rank, j, j, 
              		    point_located, rel_eps, abs_eps);
	     located_points[j]=point_located; 
	     if(!point_located){
	      determine_box( target_points, boxes, query_points, 
	      		     j,  0, rank);
	      determine_box( target_points, boxes, query_points, 
	         	     j, rank+1, size);
	     }
	  }
}

} //private functionality


class Point_search {

//public types
public:
	//all ugly and probably not necessary
	typedef typename std::vector< bool> Bitmask;
	typedef typename std::vector< unsigned int> Vector;

public:
  	//Are these necessary?
  	enum Method {LINEAR_FE, PLAIN_FE};
  	enum IntegType {VOLUME};

    /* constructor
     * Constructor, which also optionally initializes the coupler
     * \param pc ParallelComm object to be used with this coupler
     * \param local_elems Local elements in the source mesh
     * \param coupler_id Id of this coupler, should be the same over all procs
     * \param init_tree If true, initializes kdtree inside the constructor
     */
  
  Point_search( Interface & _impl,
            	ParallelComm & _pc,
            	Range & _local_elems,
            	int _id,
            	bool _initialize_tree = true,
		std::size_t _num_iterations=3): 
		impl( _impl), tree( &_impl), pc( _pc), 
		id( _id), local_elements( _local_elems),
		num_iterations( _num_iterations){
			//TODO: move into Tree constructor
	  		if (_initialize_tree) { 
				initialize_tree();
			}
		}

   // for each point, find box(es) containing the point,
   // appending results to target_points;
   // keeping local points separately, in local_points, which has pairs
   // of <local_index, mapped_index>, where mapped_index is the index
   // of <local_index, mapped_index>, where mapped_index is the index
   // into the mapped_points tuple list 
   template< typename Points, typename Tuple_list> 
   ErrorCode locate_points( Points & query_points, 
                            Tuple_list & target_points,
                            double rel_eps, 
                            double abs_eps){
   	target_points.enableWriteAccess();
 
   	// initialize source_points and local_pts
   	mapped_points(0, 0, 1, 3, target_points.get_max());
   	mapped_points.enableWriteAccess();
 
   	TupleList source_points(3, 0, 0, 0, target_points.get_max()); 
   	source_points.enableWriteAccess();
 
 
   	// keep track of which points have been located
   	Bitmask located_points(num_points, 0);
 
   	std::size_t rank = pc.proc_config().proc_rank();
 	std::size_t size = pc.proc_config().proc_size();
   	
 	moab::_point_search::distribute_data( query_points, target_points, 
 					      boxes, located_points,
 					      rank, size, 
 					      num_points, rel_eps, abs_eps);
   	
 	// perform scatter/gather, to gather points to source mesh procs
 	//TODO: Make less ugly.
   	(pc.proc_config().crystal_router())->gs_transfer(1, target_points, 0);
 
   	// test target points against my elements
   	for (unsigned i = 0; i < target_points.get_n(); i++){
   	  test_local_box( target_points.vr_wr+3*i, 
 			  source_points,
   	                  target_points.vi_rd[2*i], 
 			  target_points.vi_rd[2*i+1], 
 			  i, located_points, 
 			  rel_eps, abs_eps);
 	}
   	// no longer need target_points
 	//TODO: determine if this is really necessary.
   	target_points.reset();
 
   	 // send target points back to target procs
 	//TODO: Make less ugly.
   	(pc.proc_config().crystal_router())->gs_transfer(1, source_points, 0);
   	                                                                               
 	moab::_point_search::write_data( target_points, source_points, located_points);
 	#ifdef DEBUG_LOCATER
   		int mapped_points  = target_points->get_n() + local_mapped_points.size()/2;
   		int missing_points = num_points-mapped_points;
   		std::cerr << "point location: wanted " 
 			  << num_points << " points" 
 			  << "got " << local_mapped_points.size()/2 << "locally " 
 			  << target_points.get_n() "remote, "
 			  << "missing  " << missing_points << std::endl;
 		//TODO: WTF
 		//will litely break on curved geometries
   		assert(0==missing_points); 
 	#endif
   	return MB_SUCCESS;
   } 

  template< typename Points> 
  ErrorCode locate_points( Points & query_points,
                           double rel_eps = 0.0, 
                           double abs_eps = 0.0){
	TupleList result;
	return locate_points( query_points, result,
			      rel_eps, abs_eps);
  }

 
  template< typename Range, typename Tuple_list> 
  ErrorCode locate_points( Range & range,
                           Tuple_list & target_points,
                           double rel_eps = 0.0, 
                           double abs_eps = 0.0,
			   bool store_local=true){
  	// get locations
  	std::vector<double> query_points( 3*range.size());
  	ErrorCode rval = impl.get_coords( range, &query_points[0]);
  	if (store_local) { 
		target_vertices	= range;
	}
  	
  	return locate_points( query_points, 
  	      		      target_points, 
  	      		      rel_eps, abs_eps, 
  	      		      store_local);
  }


  Interface* 		mb_impl() 	const 	 { return &impl;	   }
  Tree*			my_tree() 	const 	 { return &tree;	  }
  EntityHandle 		local_root() 	const 	 { return local_root;	   }
  const Vector_double &	all_boxes() 	const  	 { return boxes;	   }
  ParallelComm*		my_pc() 	const 	 { return &pc;		   }
  const Range&		target_verts() 	const 	 { return target_vertices; }
  int 			my_id() 	const 	 { return id;		   }
  const Range&		my_local_elements() const 	 { return local_elements; }
  TupleList*		mapped_pts() 	const 	 { return &mapped_points;  }
  const Vector& 	local_mapped_points() const { return local_mapped_points;  }
  std::size_t 		num_its() 	const 	 { return num_iterations;  }
        
private:

//temporary hack for rhl sanity
#include< nat_param.hpp>

template< typename Query_points, 
	  typename Tuple_list>
ErrorCode test_local_box( Query_points & query_points,
			  Tuple_list & target_points, 
                          int from_proc, int remote_index, 
			  int index, 
                          bool & point_located,
                          double rel_eps=0.0, 
			  double abs_eps=0.0) {
  typedef typename std::vector< EntityHandle> Entities;
  typedef typename Entities::iterator Entity_iterator;
  typedef typename std::vector< CartVect> Coordinates;
  typedef typename Coordinates::iterator Coordinate_iterator;
  Entities entities;
  Coordinates nat_coords;
  
  if(!target_points.get_writeEnabled()){ target_points.enableWriteAccess(); }

  //TODO: incorrect conditional
  //TODO: does this need to be done on _every_ test_local_box call?   
  // relative epsilon given, translate to absolute epsilon using box dimensions
  if (rel_eps && !abs_eps) {
    CartVect minmax[2];
    tree.get_tree_box(local_root, minmax[0].array(), minmax[1].array());
    abs_eps = rel_eps * (minmax[1] - minmax[0]).length();
  }
  
  ErrorCode result = nat_param(query_points, entities, nat_coords, abs_eps);
  if (MB_SUCCESS != result) { return result; }

    // if we didn't find any ents and we're looking locally, nothing more to do
  if (entities.empty()){
      if (target_points->get_n() == target_points->get_max())P
	target_points->resize(target_points->get_max() + (1+target_points->get_max())/2);
      }
      target_points->vi_wr[3*target_points->get_n()] = from_proc;
      target_points->vi_wr[3*target_points->get_n()+1] = remote_index;
      target_points->vi_wr[3*target_points->get_n()+2] = -1;
      target_points->inc_n();
    point_located = false;
    return MB_SUCCESS;
  }

    // grow if we know we'll exceed size
  if (mapped_points->get_n()+entities.size() >= mapped_points->get_max())
    mapped_points->resize(mapped_points->get_max() + (1+mapped_points->get_max())/2);

  mapped_points.enableWriteAccess();
  Coordinate_iterator ncit = nat_coords.begin();
  for (Entity_iterator eit = entities.begin(); 
      	       eit != entities.end(); eit++, ncit++) {
      // store in tuple mapped_points
    mapped_points.vr_wr[3*mapped_points->get_n()] = (*ncit)[0];
    mapped_points.vr_wr[3*mapped_points->get_n()+1] = (*ncit)[1];
    mapped_points.vr_wr[3*mapped_points->get_n()+2] = (*ncit)[2];
    mapped_points.vul_wr[mapped_points->get_n()] = *eit;
    mapped_points.inc_n();

      // also store local point, mapped point indices
    {
      if (target_points->get_n() == target_points->get_max()){ 
	target_points->resize(target_points->get_max() + (1+target_points->get_max())/2);
	}
        // store in tuple source_points
      target_points->vi_wr[3*target_points->get_n()] = from_proc;
      target_points->vi_wr[3*target_points->get_n()+1] = remote_index;
      target_points->vi_wr[3*target_points->get_n()+2] = mapped_points->get_n()-1;
      target_points->inc_n();
    }
    else {
	local_mapped_points.push_back(index);
	local_mapped_points.push_back(mapped_points.get_n()-1);
    }
  }
  point_located = true;
  if(target_points.get_writeEnabled()){ target_points.disableWriteAccess(); }
  //if(target_points && !canWrite) target_points->disableWriteAccess();
  return MB_SUCCESS;
}
	//temporary hack for rhl sanity
	#include<initalize_tree.hpp>

	Interface & impl;
	Tree & tree;
	EntityHandle local_root;
	std::vector<double> boxes;
	ParallelComm & pc;
	int id;
	Range local_elements;
	Range target_vertices;
	TupleList & mapped_points;
	TupleList & local_target_points;
	std::vector<unsigned int> local_mapped_points;
	std::size_t num_iterations;

 // a cached spectral element for source and target , separate
 // assume that their numberof GL points (order+1) does not change
 // if it does change, we need to reinitialize it
 void * _spectralSource;
 void * _spectralTarget;
 moab::Tag _xm1Tag, _ym1Tag, _zm1Tag;

};

} // namespace moab

#endif //POINT_LOCATER_HPP
