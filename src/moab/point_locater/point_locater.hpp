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

#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}
#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}
#define ERRORMPI(a,b) {if (MPI_SUCCESS != err) {std::cerr << a << std::endl; return b;}}

#define MASTER_PROC 0

namespace moab {

class Point_search
{
public:

  enum Method {LINEAR_FE, PLAIN_FE};

  enum IntegType {VOLUME};

    /* constructor
     * Constructor, which also optionally initializes the coupler
     * \param pc ParallelComm object to be used with this coupler
     * \param local_elems Local elements in the source mesh
     * \param coupler_id Id of this coupler, should be the same over all procs
     * \param init_tree If true, initializes kdtree inside the constructor
     */
  Point_search( Interface *impl,
            	ParallelComm *pc,
            	Range &local_elems,
            	int coupler_id,
            	bool init_tree = true): 
		mbImpl(impl), myPc(pc), 
		myId(coupler_id), numIts(3) {
	  		assert(NULL != impl && (pc || !local_elems.empty()));
	    		// keep track of the local points, at least for now
	  		myRange = local_elems;
	    		// now initialize the tree
	  		if (init_tree) initialize_tree();
		}
  ErrorCode locate_points( double * xyz, 
		  	   int num_points,
                           double rel_eps = 0.0, 
                           double abs_eps = 0.0,
                           TupleList * tl = NULL,
                           bool store_local = true){
  	assert(tl || store_local);

  	// allocate tuple_list to hold point data: (p, i, , xyz), i = point index
  	TupleList target_pts;
  	target_pts.initialize(2, 0, 0, 3, num_points);
  	target_pts.enableWriteAccess();

  	// initialize source_pts and local_pts
  	TupleList source_pts;
  	mappedPts = new TupleList(0, 0, 1, 3, target_pts.get_max());
  	mappedPts->enableWriteAccess();

  	source_pts.initialize(3, 0, 0, 0, target_pts.get_max()); 
  	source_pts.enableWriteAccess();

  	mappedPts->set_n(0);
  	source_pts.set_n(0);
  	ErrorCode result;

  	  // keep track of which points have been located
  	std::vector<unsigned char> located_pts(num_points, 0);

  	  // for each point, find box(es) containing the point,
  	  // appending results to tuple_list;
  	  // keep local points separately, in local_pts, which has pairs
  	  // of <local_index, mapped_index>, where mapped_index is the index
  	  // of <local_index, mapped_index>, where mapped_index is the index
  	  // into the mappedPts tuple list

  	unsigned int my_rank = (myPc ? myPc->proc_config().proc_rank() : 0);
  	bool point_located;
  	
  	for (int i = 0; i < 3*num_points; i+=3) 
  	{
  	    // test point locally first
  	  result = test_local_box(xyz+i, my_rank, i/3, i/3, point_located, rel_eps, abs_eps);
  	  if (MB_SUCCESS != result) {
  	    return result;
  	  }
  	  if (point_located) {
  	    located_pts[i/3] = 0x1;
  	    continue;
  	  }

  	  
  	    // if not located locally, test other procs' boxes
  	  for (unsigned int j = 0; j < (myPc ? myPc->proc_config().proc_size() : 0); j++)
  	  {
  	    if (j == my_rank) continue;
  	    
  	      // test if point is in proc's box
  	    if (allBoxes[6*j] <= xyz[i] && xyz[i] <= allBoxes[6*j+3] && 
  	        allBoxes[6*j+1] <= xyz[i+1] && xyz[i+1] <= allBoxes[6*j+4] && 
  	        allBoxes[6*j+2] <= xyz[i+2] && xyz[i+2] <= allBoxes[6*j+5])
  	    {
  	        // if in this proc's box, will send to proc to test further
  	        // check size, grow if we're at max
  	      if (target_pts.get_n() == target_pts.get_max()) 
  	        target_pts.resize(target_pts.get_max() + (1+target_pts.get_max())/2);
  	
  	      target_pts.vi_wr[2*target_pts.get_n()] = j;
  	      target_pts.vi_wr[2*target_pts.get_n()+1] = i/3;

  	      target_pts.vr_wr[3*target_pts.get_n()] = xyz[i];
  	      target_pts.vr_wr[3*target_pts.get_n()+1] = xyz[i+1];
  	      target_pts.vr_wr[3*target_pts.get_n()+2] = xyz[i+2];
  	      target_pts.inc_n();
  	    }
  	  }
  	}

  	  // perform scatter/gather, to gather points to source mesh procs
  	if (myPc) {
  	  (myPc->proc_config().crystal_router())->gs_transfer(1, target_pts, 0);

  	    // after scatter/gather:
  	    // target_pts.set_n( # points local proc has to map );
  	    // target_pts.vi_wr[2*i] = proc sending point i
  	    // target_pts.vi_wr[2*i+1] = index of point i on sending proc
  	    // target_pts.vr_wr[3*i..3*i+2] = xyz of point i
  	    //
  	    // Mapping builds the tuple list:
  	    // source_pts.set_n (target_pts.get_n() )
  	    // source_pts.vi_wr[3*i] = target_pts.vi_wr[2*i] = sending proc
  	    // source_pts.vi_wr[3*i+1] = index of point i on sending proc
  	    // source_pts.vi_wr[3*i+2] = index of mapped point (-1 if not mapped)
  	    //
  	    // Also, mapping builds local tuple_list mappedPts:
  	    // mappedPts->set_n( # mapped points );
  	    // mappedPts->vul_wr[i] = local handle of mapped entity
  	    // mappedPts->vr_wr[3*i..3*i+2] = natural coordinates in mapped entity

  	    // test target points against my elements
  	  for (unsigned i = 0; i < target_pts.get_n(); i++) 
  	  {
  	    result = test_local_box(target_pts.vr_wr+3*i, 
  	                            target_pts.vi_rd[2*i], target_pts.vi_rd[2*i+1], i, 
  	                            point_located, rel_eps, abs_eps, &source_pts);
  	    if (MB_SUCCESS != result) return result;
  	  }

  	    // no longer need target_pts
  	  target_pts.reset();

  	    // send target points back to target procs
  	  (myPc->proc_config().crystal_router())->gs_transfer(1, source_pts, 0);
  	}
  	
  	// store proc/index tuples in targetPts, and/or pass back to application;
  	// the tuple this gets stored to looks like:
  	// tl.set_n( # mapped points );
  	// tl.vi_wr[3*i] = remote proc mapping point
  	// tl.vi_wr[3*i+1] = local index of mapped point
  	// tl.vi_wr[3*i+2] = remote index of mapped point
  	//
  	// Local index is mapped into either myRange, holding the handles of
  	// local mapped entities, or myXyz, holding locations of mapped pts

  	// count non-negatives
  	int num_pts = 0;
  	for (unsigned int i = 0; i < source_pts.get_n(); i++)
  	  if (-1 != source_pts.vi_rd[3*i+2]) num_pts++;  

  	  // store information about located points
  	TupleList *tl_tmp;
  	if (!store_local) 
  	  tl_tmp = tl;
  	else {
  	  targetPts = new TupleList();
  	  tl_tmp = targetPts;
  	}

  	tl_tmp->initialize(3, 0, 0, 1, num_pts);
  	tl_tmp->enableWriteAccess();

  	for (unsigned int i = 0; i < source_pts.get_n(); i++) {
  	  if (-1 != source_pts.vi_rd[3*i+2]) { //why bother sending message saying "i don't have the point" if it gets discarded?

  	    int locIndex = source_pts.vi_rd[3*i+1];
  	    if(located_pts[locIndex]){  
  	      //asked 2+ procs if they have point p, they both said yes, we'll keep the one with lowest rank
  	      //todo: check that the cases where both say yes are justified (seemed to happen too often in tests)
  	      continue;
  	    }

  	    located_pts[locIndex] = 1;

  	    tl_tmp->vi_wr[3*tl_tmp->get_n()]     = source_pts.vi_rd[3*i];
  	    tl_tmp->vi_wr[3*tl_tmp->get_n() + 1] = source_pts.vi_rd[3*i+1];
  	    tl_tmp->vi_wr[3*tl_tmp->get_n() + 2] = source_pts.vi_rd[3*i+2];
  	    tl_tmp->inc_n();
  	  }
  	}

  	int mappedPoints  = tl_tmp->get_n() + localMappedPts.size()/2;
  	int missingPoints = num_points-mappedPoints;
  	printf("point location: wanted %d got %u locally, %d remote, missing %d\n", 
  	       num_points, (uint)localMappedPts.size()/2,  tl_tmp->get_n(), missingPoints);
  	assert(0==missingPoints); //will litely break on curved geometries
  	
  	  // no longer need source_pts
  	source_pts.reset();

  	  // copy into tl if passed in and storing locally
  	if (tl && store_local) {
  	  tl = new TupleList(3, 0, 0, 1, num_pts);
  	  tl->enableWriteAccess();
  	  memcpy(tl->vi_wr, tl_tmp->vi_rd, 3*tl_tmp->get_n()*sizeof(int));
  	  tl->set_n( tl_tmp->get_n() );
  	  tl->disableWriteAccess();
  	}

  	tl_tmp->disableWriteAccess();

  	  // done
  	return MB_SUCCESS;
  }
  
  ErrorCode locate_points( Range &targ_verts,
                           double rel_eps = 0.0, 
                           double abs_eps = 0.0,
                           TupleList *tl = NULL,
                           bool store_local = true){
  	// get locations
  	std::vector<double> locs( 3*targ_verts.size());
  	ErrorCode rval = mbImpl->get_coords( targ_verts, &locs[0]);
  	if (MB_SUCCESS != rval) {
  	        return rval;
  	}

  	if (store_local) {
  	        targetVerts = targ_verts;
  	}
  	
  	return locate_points( &locs[0], 
  	      	  	      targ_verts.size(), 
  	      		      rel_eps, 
  	      		      abs_eps, 
  	      		      tl, 
  	      		      store_local);
  }

  
    /* Get functions */
  // Any good compiler should be inlining all of this for you.
  // If not, you should file a bug against the compiler
  // these are const methods, and they are small. 
  inline Interface *mb_impl() const {return mbImpl;}
  inline AdaptiveKDTree *my_tree() const {return myTree;};
  inline EntityHandle local_root() const {return localRoot;}
  inline const std::vector<double> &all_boxes() const {return allBoxes;}
  inline ParallelComm *my_pc() const {return myPc;}
  inline const Range &target_verts() const {return targetVerts;}
  inline int my_id() const {return myId;}
  inline const Range &my_range() const {return myRange;}
  inline TupleList *mapped_pts() const {return mappedPts;}
  inline const std::vector<unsigned int> &local_mapped_pts() const { 
	  return localMappedPts;
  }
  inline int num_its() const {return numIts;}
        
private:
ErrorCode nat_param( double xyz[3], 
                     std::vector<EntityHandle> &entities, 
                     std::vector<CartVect> &nat_coords,
                     double epsilon = 0.0){
  AdaptiveKDTreeIter treeiter;
  ErrorCode result = myTree->get_tree_iterator(localRoot, treeiter); 
  if (MB_SUCCESS != result) {
    std::cout << "Problems getting iterator" << std::endl;
    return result;
  }

  EntityHandle closest_leaf;
  if (epsilon) {
    std::vector<double> dists;
    std::vector<EntityHandle> leaves;
    result = myTree->leaves_within_distance(localRoot, xyz, epsilon, leaves, &dists);
    if (leaves.empty()) 
      // not found returns success here, with empty list, just like case with no epsilon
      return MB_SUCCESS;
      // get closest leaf
    double min_dist = *dists.begin();
    closest_leaf = *leaves.begin();
    std::vector<EntityHandle>::iterator vit = leaves.begin()+1;
    std::vector<double>::iterator dit = dists.begin()+1;
    for (; vit != leaves.end() && min_dist; vit++, dit++) {
      if (*dit < min_dist) {
        min_dist = *dit;
        closest_leaf = *vit;
      }
    }
  }
  else {
    result = myTree->leaf_containing_point(localRoot, xyz, treeiter);
    if(MB_ENTITY_NOT_FOUND==result) //point is outside of myTree's bounding box
      return MB_SUCCESS; 
    else if (MB_SUCCESS != result) {
      std::cout << "Problems getting leaf \n";
      return result;
    }
    closest_leaf = treeiter.handle();
  }

    // find natural coordinates of point in element(s) in that leaf
  CartVect tmp_nat_coords; 
  Range range_leaf;
  result = mbImpl->get_entities_by_dimension(closest_leaf, 3, range_leaf, false);
  if(result != MB_SUCCESS) std::cout << "Problem getting leaf in a range" << std::endl;

  // loop over the range_leaf
  for(Range::iterator iter = range_leaf.begin(); iter != range_leaf.end(); iter++)
  {
    //test to find out in which entity the point is
    //get the EntityType and create the appropriate Element::Map subtype
    // if spectral, do not need coordinates, just the GL points
    EntityType etype = mbImpl->type_from_handle(*iter);
    if (NULL!= this->_spectralSource && etype==MBHEX)
    {
      EntityHandle eh = *iter;
      const double * xval;
      const double * yval;
      const double * zval;
      ErrorCode rval = mbImpl-> tag_get_by_ptr(_xm1Tag, &eh, 1,(const void **) &xval );
      if (moab::MB_SUCCESS != rval)
      {
        std::cout << "can't get xm1 values \n";
        return MB_FAILURE;
      }
      rval = mbImpl-> tag_get_by_ptr(_ym1Tag, &eh, 1, (const void **)&yval );
      if (moab::MB_SUCCESS != rval)
      {
        std::cout << "can't get ym1 values \n";
        return MB_FAILURE;
      }
      rval = mbImpl-> tag_get_by_ptr(_zm1Tag, &eh, 1, (const void **)&zval );
      if (moab::MB_SUCCESS != rval)
      {
        std::cout << "can't get zm1 values \n";
        return MB_FAILURE;
      }
      Element::SpectralHex * spcHex = ( Element::SpectralHex * ) _spectralSource;

      spcHex->set_gl_points((double*)xval, (double*)yval, (double*)zval);
      try{
        tmp_nat_coords =spcHex->ievaluate(CartVect(xyz));
      }
      catch (Element::Map::EvaluationError) {
        std::cout << "point "<< xyz[0] << " " << xyz[1] << " " << xyz[2] <<
            " is not converging inside hex " << mbImpl->id_from_handle(eh) << "\n";
        continue; // it is possible that the point is outside, so it will not converge
      }

    }
    else
    {
      const EntityHandle *connect;
      int num_connect;

        //get connectivity
      result = mbImpl->get_connectivity(*iter, connect, num_connect, true);

        //get coordinates of the vertices
      std::vector<CartVect> coords_vert(num_connect);
      result = mbImpl->get_coords(connect, num_connect, &(coords_vert[0][0]));
      if (MB_SUCCESS != result) {
        std::cout << "Problems getting coordinates of vertices\n";
        return result;
      }

      if (etype == MBHEX) {
        Element::LinearHex hexmap(coords_vert);
        try {
          tmp_nat_coords = hexmap.ievaluate(CartVect(xyz), epsilon);
        }
        catch (Element::Map::EvaluationError) {
          continue;
        }
      }
      else if (etype == MBTET){
        Element::LinearTet tetmap(coords_vert);
        try {
          tmp_nat_coords = tetmap.ievaluate(CartVect(xyz));
        }
        catch (Element::Map::EvaluationError) {
          continue;
        }
      }
      else {
        std::cout << "Entity not Hex or Tet" << std::endl;
        continue;
      }
    }
      //if we get here then we've found the coordinates.
      //save them and the entity and return success.
    entities.push_back(*iter);
    nat_coords.push_back(tmp_nat_coords);
    return MB_SUCCESS;
  }

  //didn't find any elements containing the point
  return MB_SUCCESS;
}

ErrorCode test_local_box( double *xyz, 
                                   int from_proc, int remote_index, int index, 
                                   bool &point_located,
                                   double rel_eps=0.0, 
				   double abs_eps=0.0,
                                   TupleList *tl=NULL) {
  std::vector<EntityHandle> entities;
  std::vector<CartVect> nat_coords;
  bool canWrite;
  if (tl) {
    canWrite = tl->get_writeEnabled();
    if(!canWrite) tl->enableWriteAccess();
  }

  if (rel_eps && !abs_eps) {
      // relative epsilon given, translate to absolute epsilon using box dimensions
    CartVect minmax[2];
    myTree->get_tree_box(localRoot, minmax[0].array(), minmax[1].array());
    abs_eps = rel_eps * (minmax[1] - minmax[0]).length();
  }
  
  ErrorCode result = nat_param(xyz, entities, nat_coords, abs_eps);
  if (MB_SUCCESS != result) return result;

    // if we didn't find any ents and we're looking locally, nothing more to do
  if (entities.empty()){
    if(tl){

      if (tl->get_n() == tl->get_max())
	tl->resize(tl->get_max() + (1+tl->get_max())/2);

      tl->vi_wr[3*tl->get_n()] = from_proc;
      tl->vi_wr[3*tl->get_n()+1] = remote_index;
      tl->vi_wr[3*tl->get_n()+2] = -1;
      tl->inc_n();

    }
    point_located = false;
    return MB_SUCCESS;
  }

    // grow if we know we'll exceed size
  if (mappedPts->get_n()+entities.size() >= mappedPts->get_max())
    mappedPts->resize(mappedPts->get_max() + (1+mappedPts->get_max())/2);;


  std::vector<EntityHandle>::iterator eit = entities.begin();
  std::vector<CartVect>::iterator ncit = nat_coords.begin();

  mappedPts->enableWriteAccess();
  for (; eit != entities.end(); eit++, ncit++) {
      // store in tuple mappedPts
    mappedPts->vr_wr[3*mappedPts->get_n()] = (*ncit)[0];
    mappedPts->vr_wr[3*mappedPts->get_n()+1] = (*ncit)[1];
    mappedPts->vr_wr[3*mappedPts->get_n()+2] = (*ncit)[2];
    mappedPts->vul_wr[mappedPts->get_n()] = *eit;
    mappedPts->inc_n();

      // also store local point, mapped point indices
    if (tl) 
    {
      if (tl->get_n() == tl->get_max()) 
	tl->resize(tl->get_max() + (1+tl->get_max())/2);

        // store in tuple source_pts
      tl->vi_wr[3*tl->get_n()] = from_proc;
      tl->vi_wr[3*tl->get_n()+1] = remote_index;
      tl->vi_wr[3*tl->get_n()+2] = mappedPts->get_n()-1;
      tl->inc_n();
    }
    else {
      localMappedPts.push_back(index);
      localMappedPts.push_back(mappedPts->get_n()-1);
    }
  }

  point_located = true;
  
  if(tl && !canWrite) tl->disableWriteAccess();

  return MB_SUCCESS;
}

  ErrorCode initialize_tree(){
	  Range local_ents;
	  AdaptiveKDTree::Settings settings;
	  settings.candidatePlaneSet = AdaptiveKDTree::SUBDIVISION;
	
	    //get entities on the local part
	  ErrorCode result = MB_SUCCESS;
	  if (myPc) result = myPc->get_part_entities( local_ents, 
			  				3);
	  else local_ents = myRange;

	  if (MB_SUCCESS != result || local_ents.empty()) {
	    std::cout << "Problems getting source entities" 
		      << std::endl;
	    return result;
	  }
	
	    // build the tree for local processor
	  for (int i = 0; i < numIts; i++) {
	    myTree = new AdaptiveKDTree( mbImpl);
	    result = myTree->build_tree( local_ents, 
			    		 localRoot, 
					 &settings);
	    if (MB_SUCCESS != result) {
	      std::cout << "Problems building tree";
	      if (numIts != i) {
	        delete myTree;
	        settings.maxEntPerLeaf *= 2;
	        std::cout << "; increasing elements/leaf to " 
	                  << settings.maxEntPerLeaf << std::endl;;
	      }
	      else {
	        std::cout << "; exiting" << std::endl;
	        return result;
	      }
	    }
	    else
	      break; // get out of tree building
	  }
	
	    // get the bounding box for local tree
	  if (myPc){
	    allBoxes.resize(6*myPc->proc_config().proc_size());
	  }
	  else{
		 allBoxes.resize(6);
	  }
	  unsigned int my_rank = (myPc ? myPc->proc_config().proc_rank() : 0);
	  result = myTree->get_tree_box( localRoot, 
			  		 &allBoxes[6*my_rank], 
					 &allBoxes[6*my_rank+3]);
	  if (MB_SUCCESS != result) return result;
	  
	    // now communicate to get all boxes
	    // use "in place" option
	  if (myPc) {
	    int mpi_err = MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
	                                &allBoxes[0], 6, MPI_DOUBLE, 
	                                myPc->proc_config().proc_comm());
	    if (MPI_SUCCESS != mpi_err) return MB_FAILURE;
	  }
	
	  /*  std::ostringstream blah;
	  for(int i=0; i<allBoxes.size(); i++)
	  blah << allBoxes[i] << " ";
	  std::cout<<blah.str()<<"\n";*/
	
	
	#ifndef NDEBUG
	  double min[3] = {0,0,0}, max[3] = {0,0,0};
	  unsigned int dep;
	  myTree->get_info(localRoot, min, max, dep);
	  std::cout << "Proc " << my_rank << ": box min/max, tree depth = ("
	            << min[0] << "," << min[1] << "," << min[2] << "), ("
	            << max[0] << "," << max[1] << "," << max[2] << "), "
	            << dep << std::endl;
	#endif  
	
	  return result;
  }
  Interface *mbImpl;
  AdaptiveKDTree * myTree;
  EntityHandle localRoot;
  std::vector<double> allBoxes;
  ParallelComm *myPc;
  int myId;
  Range myRange;
  Range targetVerts;
  TupleList *mappedPts;
  TupleList *targetPts;
  std::vector<unsigned int> localMappedPts;
  int numIts;

 // a cached spectral element for source and target , separate
 // assume that their numberof GL points (order+1) does not change
 // if it does change, we need to reinitialize it
 void * _spectralSource;
 void * _spectralTarget;
 moab::Tag _xm1Tag, _ym1Tag, _zm1Tag;

};

} // namespace moab

#endif //POINT_LOCATER_HPP
