#include "moab/ParallelMergeMesh.hpp"
#include "moab/Core.hpp"
#include "moab/CartVect.hpp"
#include "moab/Range.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/Skinner.hpp"
#include "moab/MergeMesh.hpp"
#include "moab/ParallelComm.hpp"

extern "C" 
{
  //#include "errmem.h"
  //#include "types.h"
  //#include "sort.h"
  //#include "tuple_list.h"
#include "crystal.h"
#include "transfer.h"
}

#include <algorithm>
#include "moab_mpi.h"


namespace moab{

  //Constructor
  /*1) Get Merge Data and tolerance*/
  ParallelMergeMesh::ParallelMergeMesh(ParallelComm *pc, 
				       const double epsilon) :
    myPcomm(pc), myEps(epsilon)
  {
  }

  //Perform Merge
  ErrorCode ParallelMergeMesh::merge() 
  {
    Interface *mb = myPcomm->get_moab();
    int dim;
    ErrorCode rval = mb->get_dimension(dim);
    if(rval != MB_SUCCESS){
      return rval;
    }

    /*2)Merge Mesh Locally*/
    //Get all dim dimensional entities
    Range ents;
    rval = mb->get_entities_by_dimension(0,dim,ents);
    if(rval != MB_SUCCESS){
      return rval;
    }

    //Merge Mesh Locally
    MergeMesh merger(mb);
    merger.merge_entities(ents,myEps);
    if(rval != MB_SUCCESS){
      return rval;
    }

    //Rebuild the ents range
    ents.clear();
    rval = mb->get_entities_by_dimension(0,dim,ents);
    if(rval != MB_SUCCESS){
      return rval;
    }

    /*3)Get Skin
      -Get Range of 0 dimensional entities*/
    Range skinents;
    Skinner skinner(mb);
    rval = skinner.find_skin(ents,0,skinents);
    if(rval != MB_SUCCESS){
      return rval;
    }
 
    /*4)Get Bounding Box*/
    AdaptiveKDTree kd(mb);
    double box[6], gbox[6];
    rval = kd.bounding_box(skinents,box, box+dim);
    if(rval != MB_SUCCESS){
      return rval;
    }

    //Invert the max
    for(int i=dim; i<dim*2; i++){
      box[i] *= -1;
    }

    /* 5)Communicate to all processors*/
    MPI_Allreduce(box, gbox, dim*2, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    /* 6)Assemble Global Bounding Box*/
    //Flip the max back
    for(int i=dim; i<dim*2; i++){
      gbox[i] *= -1;
    }

    /* 7)Partition Global Bounding Box 
       -Simple 1 Dimensional Partition
       -An element (x,y,z) will be sent to floor((x-min(x))/length) */
    double length = gbox[dim]-gbox[0];
    length /= static_cast<double>(myPcomm->size());

    /* 8)Get Skin Coordinates, Vertices */
    double *x = new double[skinents.size()]; 
    double *y = new double[skinents.size()]; 
    double *z = new double[skinents.size()];
    rval = mb->get_coords(skinents,x,y,z);
    if(rval != MB_SUCCESS){
      //Prevent Memory Leak
      delete []x; delete []y; delete []z;
      return rval;
    }

    /* 9)Assemble Tuples */
    double eps2 = myEps*myEps;
    tuple_list tup;
    tuple_list_init_max(&tup,1,0,1,3,skinents.size());
    int maxProc = myPcomm->size()-1;
    unsigned long j=0, tup_i=0,tup_ul=0,tup_r=0;
    for(Range::iterator it = skinents.begin();it != skinents.end();it++){
      //Calculate the toprocessor
      int toproc1 = static_cast<int>(floor((x[j]-gbox[0])/length));    
      int toproc2 = static_cast<int>(floor((x[j]+eps2-gbox[0])/length));
      
      //Make sure no entities go to an invalid processor
      toproc1 = toproc1<=maxProc?toproc1:maxProc;
      toproc2 = toproc2<=maxProc?toproc2:maxProc;

      //Make sure we have room for at least 2 more tuples
      while(tup.n+1 >= tup.max){
	tuple_list_grow(&tup);
      }

      //Assemble a tuple
      //Tuple is of the form (toProc, handle, x,y, z)
      tup.vi[tup_i++]=toproc1; //toProc
      tup.vul[tup_ul++]=*it; //Handle
      tup.vr[tup_r++]=x[j]; //x
      tup.vr[tup_r++]=y[j]; //y
      tup.vr[tup_r++]=z[j]; //z
      tup.n++;
      //Assemble a second tuple if necessary
      if(toproc1 != toproc2){
	tup.vi[tup_i++]=toproc2;// 2nd toProc
	tup.vul[tup_ul++]=*it;// handle
	tup.vr[tup_r++]=x[j];// x 
	tup.vr[tup_r++]=y[j];// y
	tup.vr[tup_r++]=z[j];// z
	tup.n++;
      }
      j++;
    }

    //Delete the coord arrays
    delete []x; delete []y; delete []z;

    /* 10)Gather-Scatter Tuple
       -tup comes out as (remProc,handle,x,y,z) */
    crystal_data cd; 
    moab_crystal_init(&cd,myPcomm->comm());
    //1 represents dynamic tuple, 0 represents index of the processor to send to
    moab_gs_transfer(1, &tup, 0, &cd);
  
    /* 11)Sort By X,Y,Z
       -Currently Utilizing a custom sort*/
    buffer buf;
    unsigned long long max_size = 
      skinents.size()*(MAX_SHARING_PROCS+1)*sizeof(double);
    buffer_init(&buf, max_size);
    //Sort by x,y,z
    temp_tuple_sort_real(&tup,eps2);

    /* 12)Match: new tuple list*/
    tuple_list matches;
    tuple_list_init_max(&matches,2,0,2,0,skinents.size());
  
    //Counters for accessing tuples more efficiently
    unsigned long i = 0, mat_i=0, mat_ul=0;
    j=0; tup_r=0;
    if(myPcomm->rank()==0 && false){
      std::cout<<"Ver6"<<std::endl;
    }
    while((i+1)<tup.n){
      //Proximity Comparison
      double xi = tup.vr[tup_r], 
	yi = tup.vr[tup_r+1],
	zi = tup.vr[tup_r+2];

      bool done = false;
      while(!done && (j+1)<tup.n){
	j++; tup_r+= tup.mr;
	CartVect cv(tup.vr[tup_r]-xi,
		    tup.vr[tup_r+1]-yi,
		    tup.vr[tup_r+2]-zi);
	if(cv.length_squared() > eps2){
	  done = true;
	}
      }

      //Allocate the tuple list before adding matches
      while(matches.n+(j-i)*(j-i-1) >= matches.max){
	tuple_list_grow(&matches);
      }
 
      //We now know that tuples [i to j) exclusice match.  
      //If n tuples match, n*(n-1) match tuples will be made
      //tuples are of the form (proc1,proc2,handle1,handle2)
      int kproc = i*tup.mi;
      unsigned khand = i*tup.mul;
      for(unsigned k = i; k<j; k++){
	int lproc = kproc+tup.mi;
	unsigned lhand = khand+tup.mul;
	for(unsigned l=k+1; l<j; l++){
	  matches.vi[mat_i++]=tup.vi[kproc];//proc1
	  matches.vi[mat_i++]=tup.vi[lproc];//proc2
	  matches.vul[mat_ul++]=tup.vul[khand];//handle1
	  matches.vul[mat_ul++]=tup.vul[lhand];//handle2
	  matches.n++;

	  matches.vi[mat_i++]=tup.vi[lproc];//proc1
	  matches.vi[mat_i++]=tup.vi[kproc];//proc2
	  matches.vul[mat_ul++]=tup.vul[lhand];//handle1
	  matches.vul[mat_ul++]=tup.vul[khand];//handle2
	  matches.n++;
	  lproc += tup.mi;
	  lhand += tup.mul;
	}
	kproc += tup.mi;
	khand += tup.mul;
      }//End for(int k...
      i = j;
    }//End while(i<tup.n)
    //Cleanup

    tuple_list_free(&tup);
    /* 13)Gather-Scatter Again */
  
    //Transfer
    //1 represents dynamic list, 0 represents proc index to send tuple to
    moab_gs_transfer(1,&matches,0,&cd);
    moab_crystal_free(&cd);

    //Sorts are necessary to check for doubles
    //Sort by remote handle
    moab_tuple_list_sort(&matches,3,&buf);
    //Sort by matching proc
    moab_tuple_list_sort(&matches,1,&buf);
    //Sort by local handle
    moab_tuple_list_sort(&matches,2,&buf);
    buffer_free(&buf);

    /* 14)get vertex/tag it
       -Using update_remote_data
       -Still only using update_remote_data one match at a time
       -I have a different version written, but there may be a complication*/
    unsigned long n = 0; 
    mat_i = 0; mat_ul = 0;
    while(n < matches.n){
      //Assert that we do not have a local match
      ulong localHand = matches.vul[mat_ul];
      sint *proc = &matches.vi[mat_i+1];
      ulong *remHand = &matches.vul[mat_ul+1];
      //Make sure we didn't miss any local shared ents
      assert(*proc != myPcomm->rank());

      //Use update_remote_data to apply tag data
      myPcomm->update_remote_data((EntityHandle)localHand,
				  (int *)proc,
				  (EntityHandle *)remHand,
				  1,1);

      //Look for doubles in this list, move past them
      //Shouldn't happen in most cases, but not impossible
      n++; mat_ul+= matches.mul; mat_i += matches.mi;
      while(n < matches.n && 
	    localHand == matches.vul[mat_ul] &&//[n*matches.mul] &&
	    *proc == matches.vi[mat_i+1] &&    //[n*matches.mi+1] && 
	    *remHand == matches.vul[mat_ul+1]){//[n*matches.mul+1]){
	n++; mat_ul+= matches.mul; mat_i += matches.mi;
      }
    }
    tuple_list_free(&matches);
    /* 15)Match higher-dimensional entities

       -Copying parameters from resolve_shared_ents*/
    rval = myPcomm->create_interface_sets(dim,dim-1);
    if(rval != MB_SUCCESS){
      return rval;
    }

    std::set<unsigned int> psets;
    rval = myPcomm->get_interface_procs(psets, true);
    if(rval != MB_SUCCESS){
      return rval;
    }

    //rval = myPcomm->exchange_ghost_cells(-1,-1,0,true,true);
    //if(rval != MB_SUCCESS){
    //return rval;
    //}

    rval = myPcomm->create_iface_pc_links();
    return rval;
  }

  //Swap around tuples
  void ParallelMergeMesh::temp_tuple_swap_real(tuple_list *tup, 
					       unsigned a, unsigned b)
  {
    if(a==b) return;

    for(unsigned int i=0; i<tup->mi;i++){
      sint t =tup->vi[a*tup->mi+i];
      tup->vi[a*tup->mi+i] = tup->vi[b*tup->mi+i];
      tup->vi[b*tup->mi+i] = t; 
    }
    for(unsigned int i=0; i<tup->ml;i++){
      slong t =tup->vl[a*tup->ml+i];
      tup->vl[a*tup->ml+i] = tup->vl[b*tup->ml+i];
      tup->vl[b*tup->ml+i] = t; 
    }
    for(unsigned int i=0; i<tup->mul;i++){
      ulong t =tup->vul[a*tup->mul+i];
      tup->vul[a*tup->mul+i] = tup->vul[b*tup->mul+i];
      tup->vul[b*tup->mul+i] = t; 
    }
    for(unsigned int i=0; i<tup->mr;i++){
      real t =tup->vr[a*tup->mr+i];
      tup->vr[a*tup->mr+i] = tup->vr[b*tup->mr+i];
      tup->vr[b*tup->mr+i] = t; 
    }
  }

  //Simple selection sort to test real
  void ParallelMergeMesh::temp_tuple_sort_real(tuple_list *tup, 
						double eps2)
  {
    for(unsigned int i=0;i<tup->n-1;i++){
      int min = i;
      for(unsigned int j = i+1; j<tup->n;j++){
	unsigned p = 0;
	while(p < tup->mr){
	  if(fabs(tup->vr[(min*tup->mr)+p]-tup->vr[(j*tup->mr)+p]) <= eps2){
	    p++;
	    continue;
	  }
	  if(tup->vr[(min*tup->mr)+p] > tup->vr[(j*tup->mr)+p]){
	    min = j;
	    break;
	  }
	  else{
	    break;
	  }
	}
      }//End for j
      temp_tuple_swap_real(tup,i,min);
    }//End for i
  }

}//End namespace moab

