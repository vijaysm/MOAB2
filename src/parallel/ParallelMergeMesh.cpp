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
#include "errmem.h"
#include "types.h"
#include "sort.h"
#include "tuple_list.h"
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
  Range ents;
  rval = mb->get_entities_by_dimension(0,dim,ents);
  if(rval != MB_SUCCESS){
    return rval;
  }

  //Merge Mesh Locally
  MergeMesh *merger = new MergeMesh(mb);
  merger->merge_entities(ents,myEps);
  if(rval != MB_SUCCESS){
    delete merger;
    return rval;
  }
  delete merger;

  //Rebuild the ents range
  ents.clear();
  rval = mb->get_entities_by_dimension(0,dim,ents);
  if(rval != MB_SUCCESS){
    return rval;
  }

  /*3)Get Skin
    -Get Range of 0 dimensional entities
    -Quicker to skin first?*/
  Range skinents;
  Skinner *skinner = new Skinner(mb);
  rval = skinner->find_skin(ents,0,skinents, false);
  if(rval != MB_SUCCESS){
    delete skinner;
    return rval;
  }
  delete skinner;
 
  /*4)Get Bounding Box*/
  AdaptiveKDTree kd(mb);
  double min[6], gmin[3], max[6], gmax[3];
  rval = kd.bounding_box(skinents,min,min+3);
  if(rval != MB_SUCCESS){
    return rval;
  }

  /* 5)Communicate to all processors*/
  /* 6)Assemble Global Bounding Box*/
// tjt - put all into the same vector, then do 1 allreduce
//  MPI_Allreduce(min, gmin, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
//  MPI_Allreduce(max, gmax, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(min, gmin, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  /* 7)Partition Global Bounding Box 
     -Simple 1 Dimensional Partition for now
     -An element (x,y,z) will be sent to floor((x-min(x))/length) */
  double length = gmax[0]-gmin[0];
  length /= static_cast<double>(myPcomm->size());

  /* 8)Get Skin Coordinates, Vertices */
  double *x = new double[skinents.size()]; 
  double *y = new double[skinents.size()]; 
  double *z = new double[skinents.size()];

  //Want to use vector, but doesn't appear doable for a range
  rval = mb->get_coords(skinents,x,y,z);
  if(rval != MB_SUCCESS){
    //Prevent Memory Leak
    delete []x; delete []y; delete []z;
    return rval;
  }

  /* 9)Assemble Tuples */
  tuple_list tup;
  tuple_list_init_max(&tup,1,0,1,3,skinents.size()*2);
  int maxProc = myPcomm->size()-1;
  unsigned long j=0, tup_i=0,tup_ul=0,tup_r=0;
  for(Range::iterator it = skinents.begin();it != skinents.end();it++){
    //Calculate the toprocessor
    int toproc1 = static_cast<int>(floor((x[j]-gmin[0])/length));
    int toproc2 = static_cast<int>(floor((x[j]+myEps-gmin[0])/length));

    //Make sure no entities go to an invalid processor
    toproc1 = toproc1<=maxProc?toproc1:maxProc;
    toproc2 = toproc2<=maxProc?toproc2:maxProc;

    //Make sure we have room for at least 2 more tuples
    if(tup.n+1 >= tup.max){
      tuple_list_grow(&tup);
    }

    //Assemble a tuple
    //Tuple is of the form (toProc, handle, x,y, z)
    tup.vi[tup_i++]=toproc1;//tup.vi[tup.n*tup.mi]=toproc1; 
    tup.vul[tup_ul++]=*it;//tup.vul[tup.n*tup.mul]=*it;
    tup.vr[tup_r++]=x[j];//tup.vr[tup.n*tup.mr]=x[j]; 
    tup.vr[tup_r++]=y[j];//tup.vr[tup.n*tup.mr+1]=y[j]; 
    tup.vr[tup_r++]=z[j];//tup.vr[tup.n*tup.mr+2]=z[j];
    tup.n++;

    //Assemble a second tuple if necessary
    if(toproc1 != toproc2){
      tup.vi[tup_i++]=toproc2;//tup.vi[tup.n*tup.mi]=toproc1; 
      tup.vul[tup_ul++]=*it;//tup.vul[tup.n*tup.mul]=*it;
      tup.vr[tup_r++]=x[j];//tup.vr[tup.n*tup.mr]=x[j]; 
      tup.vr[tup_r++]=y[j];//tup.vr[tup.n*tup.mr+1]=y[j]; 
      tup.vr[tup_r++]=z[j];//tup.vr[tup.n*tup.mr+2]=z[j];
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
     -Sorting is stable according to sort_imp.c*/
  buffer buf;
  int max_size = skinents.size()*(MAX_SHARING_PROCS+1);
  buffer_init(&buf, max_size);
  moab_tuple_list_sort(&tup,5,&buf);//Sort by z
  moab_tuple_list_sort(&tup,4,&buf);//y
  moab_tuple_list_sort(&tup,3,&buf);//x

  /* 12)Match: new tuple list*/
  tuple_list matches;
  tuple_list_init_max(&matches,2,0,2,0,skinents.size());
  double eps2 = myEps*myEps;
  
  //Counters for accessing tuples more efficiently
  unsigned long i = 0, mat_i=0, mat_ul=0;
  j=0; tup_r=0;

  while(i<tup.n){
    //Proximity Comparison
    double xi = tup.vr[tup_r], 
           yi = tup.vr[tup_r+1],
           zi = tup.vr[tup_r+2];

    bool done = false;
    while(j<tup.n && !done){
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
    for(unsigned k = i; k<j; k++){
      int kproc = k*tup.mi, lproc = i*tup.mi;
      unsigned khand = k*tup.mul, lhand = i*tup.mul;
      for(unsigned l=i; l<j; l++){
	if(l != k){
	  /*matches.vi[matches.n*matches.mi]=tup.vi[k*tup.mi];//proc1
	  matches.vi[matches.n*matches.mi+1]=tup.vi[l*tup.mi];//proc2
	  matches.vul[matches.n*matches.mul]=tup.vul[k*tup.mul];//handle1
	  matches.vul[matches.n*matches.mul+1]=tup.vul[l*tup.mul];//handle2*/

	  //Note that tup.mi == tup.mul == 0
	  matches.vi[mat_i++]=tup.vi[kproc];//proc1
	  matches.vi[mat_i++]=tup.vi[lproc];//proc2
	  matches.vul[mat_ul++]=tup.vul[khand];//handle1
	  matches.vul[mat_ul++]=tup.vul[lhand];//handle2
	  matches.n++;
	}
	lproc += tup.mi;
	lhand += tup.mul;
      }
    }//End for(int k...
    i=j;
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
     -Using update_remote_data*/
  //An extra bool to monitor if duplicates are created
  bool hasDuplicates = false;
  unsigned long n = 0; 
  mat_i = 0; mat_ul = 0;
  while(n < matches.n){
    //Assert that we do not have a local match
    ulong localHand = matches.vul[mat_ul];//[n*matches.mul];
    sint *proc = &matches.vi[mat_i+1];
    ulong *remHand = &matches.vul[mat_ul+1];
    //Make sure we didn't miss any local shared ents
    assert(*proc != myPcomm->get_id());

    //Use update_remote_data to apply tag data
    myPcomm->update_remote_data((EntityHandle)localHand,
				(int *)proc,
				(EntityHandle *)remHand,
				1,1);

    //Look for doubles in this list, move past them
    n++; mat_ul+= matches.mul; mat_i += matches.mi;
    while(n < matches.n && 
	  localHand == matches.vul[mat_ul] &&//[n*matches.mul] &&
	  *proc == matches.vi[mat_i+1] &&//[n*matches.mi+1] && 
	  *remHand == matches.vul[mat_ul+1]){//[n*matches.mul+1]){
      hasDuplicates = true;
      n++; mat_ul+= matches.mul; mat_i += matches.mi;
    }
  }
  if(hasDuplicates){
    std::cerr<<"Duplicates were in the matched list"<<std::endl;
  }
  tuple_list_free(&matches);
  
  /* 15)Match higher-dimensional entities
     -Copying parameters from resolve_shared_ents*/
  rval = myPcomm->exchange_ghost_cells(-1,-1,0,true,true);
  if(rval != MB_SUCCESS){
    return rval;
  }
  return MB_SUCCESS;
}

}//End namespace moab
