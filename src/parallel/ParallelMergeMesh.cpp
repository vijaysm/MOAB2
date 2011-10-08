#include "moab/ParallelMergeMesh.hpp"
#include "moab/Core.hpp"
#include "moab/CartVect.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/Skinner.hpp"
#include "moab/MergeMesh.hpp"
#include "moab/CN.hpp"
#include "float.h"
#include <algorithm>

#ifdef USE_MPI
  #include "moab_mpi.h"
#endif

namespace moab{

  //Constructor
  /*Get Merge Data and tolerance*/
  ParallelMergeMesh::ParallelMergeMesh(ParallelComm *pc, 
				       const double epsilon) :
    myPcomm(pc), myEps(epsilon)
  {
    myMB = pc->get_moab();
    mySkinEnts.resize(4);
    cdAllocated = false;
  }

  
  //Have a wrapper function on the actual merge to avoid memory leaks
  //Merges elements within a proximity of epsilon
  ErrorCode ParallelMergeMesh::merge() 
  {
    ErrorCode rval = PerformMerge();
    CleanUp();
    return rval;
  }

  //Perform the merge
  ErrorCode ParallelMergeMesh::PerformMerge()
  {
    //Get the mesh dimension
    int dim;
    ErrorCode rval = myMB->get_dimension(dim);
    if(rval != MB_SUCCESS){
      return rval;
    }
    
    //Get the local skin elements
    rval = PopulateMySkinEnts(dim);
    if(rval != MB_SUCCESS){
      return rval;
    }

    //Determine the global bounding box
    double gbox[6];
    rval = GetGlobalBox(gbox);
    if(rval != MB_SUCCESS){
      return rval;
    }

    /* Assemble The Destination Tuples */
    //Get a list of tuples which contain (toProc, handle, x,y,z)
    tuple_list_init_max(&myTup,1,0,1,3,mySkinEnts[0].size());
    rval = PopulateMyTup(gbox);
    if(rval != MB_SUCCESS){
      return rval;
    }

    /* Gather-Scatter Tuple
       -tup comes out as (remoteProc,handle,x,y,z) */
    moab_crystal_init(&myCD, myPcomm->comm());
    cdAllocated = true;

    //1 represents dynamic tuple, 0 represents index of the processor to send to
    moab_gs_transfer(1, &myTup, 0, &myCD);

    /* Sort By X,Y,Z
       -Utilizes a custom quick sort incoroporating eplison*/
    SortTuplesByReal(&myTup,myEps);

    //Initilize another tuple list for matches
    tuple_list_init_max(&myMatches,2,0,2,0,mySkinEnts[0].size());

    //ID the matching tuples
    rval = PopulateMyMatches();
    if(rval != MB_SUCCESS){
      return rval;
    }

    //We can free up the tuple myTup now
    tuple_list_free(&myTup);
    myTup.max = 0;

    /*Gather-Scatter Again*/
    //1 represents dynamic list, 0 represents proc index to send tuple to
    moab_gs_transfer(1,&myMatches,0,&myCD);
    //We can free up the crystal router now
    moab_crystal_free(&myCD);
    cdAllocated = false;

    //Sort the matches tuple list
    SortMyMatches();

    //Tag the shared elements
    rval = TagSharedElements(dim);

    //Free up the matches tuples
    tuple_list_free(&myMatches);
    myMatches.max = 0;

    return rval;
  }

  //Sets mySkinEnts with all of the skin entities on the processor
  ErrorCode ParallelMergeMesh::PopulateMySkinEnts(int dim)
  {
    /*Merge Mesh Locally*/
    //Get all dim dimensional entities
    Range ents;
    ErrorCode rval = myMB->get_entities_by_dimension(0,dim,ents);
    if(rval != MB_SUCCESS){
      return rval;
    }

    //Merge Mesh Locally
    MergeMesh merger(myMB, false);
    merger.merge_entities(ents,myEps);
    //We can return if there is only 1 proc
    if(rval != MB_SUCCESS || myPcomm->size() == 1){
      return rval;
    }

    //Rebuild the ents range
    ents.clear();
    rval = myMB->get_entities_by_dimension(0,dim,ents);
    if(rval != MB_SUCCESS){
      return rval;
    }

    /*Get Skin
      -Get Range of all dimensional entities
      -skinEnts[i] is the skin entities of dimension i*/  
    Skinner skinner(myMB);
    for(int skin_dim = dim; skin_dim >= 0; skin_dim--){
      rval = skinner.find_skin(ents,skin_dim,mySkinEnts[skin_dim]);
      if(rval != MB_SUCCESS){
	return rval;
      }
    }
    return MB_SUCCESS;
  }

 //Determine the global assembly box
  ErrorCode ParallelMergeMesh::GetGlobalBox(double *gbox)
  {
    ErrorCode rval;

    /*Get Bounding Box*/
    double box[6];
    if(mySkinEnts[0].size() != 0){
      AdaptiveKDTree kd(myMB);
      rval = kd.bounding_box(mySkinEnts[0],box, box+3);
      if(rval != MB_SUCCESS){
	return rval;
      }
    }
    //If there are no entities...
    else{
      for(int i=0;i<6;i++){
	if(i < 3){
	  box[i] = DBL_MAX;
	}
	else{
	  box[i] = -DBL_MAX;
	}
      }
    }

    //Invert the max
    for(int i=3; i<6;i++){
      box[i] *= -1;
    }

    /*Communicate to all processors*/
    MPI_Allreduce(box, gbox, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    /*Assemble Global Bounding Box*/
    //Flip the max back
    for(int i=3; i<6; i++){
      gbox[i] *= -1;
    }
    return MB_SUCCESS;
  }

  //Assemble the tuples with their processor destination
  ErrorCode ParallelMergeMesh::PopulateMyTup(double * gbox)
  {
    /*Figure out how do partition the global box*/
    double lengths[3];
    int parts[3];
    ErrorCode rval = PartitionGlobalBox(gbox, lengths, parts);
    if(rval != MB_SUCCESS){
      return rval;
    }

    /* Get Skin Coordinates, Vertices */
    double *x = new double[mySkinEnts[0].size()]; 
    double *y = new double[mySkinEnts[0].size()]; 
    double *z = new double[mySkinEnts[0].size()];
    rval = myMB->get_coords(mySkinEnts[0],x,y,z);
    if(rval != MB_SUCCESS){
      //Prevent Memory Leak
      delete []x; delete []y; delete []z;
      return rval;
    }    

    //Initialize variable to be used in the loops
    std::vector<int> toProcs;
    int xPart, yPart, zPart, xEps, yEps, zEps, baseProc;
    unsigned long long tup_i=0, tup_ul=0, tup_r=0, count=0;
    //These are boolean to determine if the vertice is on close enought to a given border
    bool xDup, yDup, zDup;
    //Go through each vertice
    for(Range::iterator it = mySkinEnts[0].begin(); it != mySkinEnts[0].end(); it++){
      xDup = false; yDup = false; zDup = false;
      //Figure out which x,y,z partition the element is in.
      xPart = static_cast<int>(floor((x[count]-gbox[0])/lengths[0]));
      xPart = (xPart<parts[0]?xPart:parts[0]-1);//Make sure it stays within the bounds

      yPart = static_cast<int>(floor((y[count]-gbox[1])/lengths[1]));
      yPart = (yPart<parts[1]?yPart:parts[1]-1);//Make sure it stays within the bounds

      zPart = static_cast<int>(floor((z[count]-gbox[2])/lengths[2]));
      zPart = (zPart<parts[2]?zPart:parts[2]-1);//Make sure it stays within the bounds

      //Figure out the partition with the addition of Epsilon
      xEps = static_cast<int>(floor((x[count]-gbox[0]+myEps)/lengths[0]));
      yEps = static_cast<int>(floor((y[count]-gbox[1]+myEps)/lengths[1]));
      zEps = static_cast<int>(floor((z[count]-gbox[2]+myEps)/lengths[2]));

      //Figure out if the vertice needs to be sent to multiple procs
      xDup = (xPart != xEps && xEps < parts[0]);
      yDup = (yPart != yEps && yEps < parts[1]);
      zDup = (zPart != zEps && zEps < parts[2]);

      //Add appropriate processors to the vector
      baseProc = xPart+ yPart * parts[0] + zPart * parts[0] * parts[1]; 
      toProcs.push_back(baseProc);
      if(xDup){
	toProcs.push_back(baseProc + 1);//Get partition to the right
      }
      if(yDup){
	toProcs.push_back(baseProc + parts[0]);//Partition up 1
      }
      if(zDup){
	toProcs.push_back(baseProc + parts[0] * parts[1]);//Partition above 1
      }
      if(xDup && yDup){
	toProcs.push_back(baseProc + parts[0] + 1);//Partition up 1 and right 1
      }
      if(xDup && zDup){
	toProcs.push_back(baseProc + parts[0] * parts[1] + 1);//Partition right 1 and above 1
      }
      if(yDup && zDup){
	toProcs.push_back(baseProc + parts[0] * parts[1] + parts[0]);//Partition up 1 and above 1
      }
      if(xDup && yDup && zDup){
	toProcs.push_back(baseProc + parts[0] * parts[1] + parts[0] + 1);//Partition right 1, up 1, and above 1
      }
      //Grow the tuple list if necessary
      while(myTup.n + toProcs.size() >= myTup.max){
	tuple_list_grow(&myTup);
      }
      //Add each proc as a tuple
      for(std::vector<int>::iterator proc = toProcs.begin();
	  proc != toProcs.end();
	  proc++){
	myTup.vi[tup_i++] = *proc;
	myTup.vul[tup_ul++] = *it;
	myTup.vr[tup_r++] = x[count];
	myTup.vr[tup_r++] = y[count];
	myTup.vr[tup_r++] = z[count];
	myTup.n++;
      }

      count++;
      toProcs.clear();
    }
    return MB_SUCCESS;
  }

  //Partition the global box by the number of procs
  ErrorCode ParallelMergeMesh::PartitionGlobalBox(double *gbox, double *lengths, int *parts)
  {
    //Determine the length of each side
    double xLen = gbox[3]-gbox[0];
    double yLen = gbox[4]-gbox[1];
    double zLen = gbox[5]-gbox[2];
    unsigned numProcs = myPcomm->size();
    
    //Partition sides from the longest to shortest lengths
    //If x is the longest side
    if(xLen >= yLen && xLen >= zLen){
      parts[0] = PartitionSide(xLen, yLen * zLen, numProcs, true);
      numProcs /= parts[0];
      //If y is second longest
      if(yLen >= zLen){
	parts[1] = PartitionSide(yLen, zLen, numProcs, false);
	parts[2] = numProcs/parts[1];
      }
      //If z is the longer
      else{
	parts[2] = PartitionSide(zLen, yLen, numProcs, false);
	parts[1] = numProcs/parts[2];
      }
    }
    //If y is the longest side
    else if (yLen >= zLen){
      parts[1] = PartitionSide(yLen, xLen * zLen, numProcs, true);
      numProcs /= parts[1];
      //If x is the second longest
      if(xLen >= zLen){
	parts[0] = PartitionSide(xLen, zLen, numProcs, false);
	parts[2] = numProcs/parts[0];
      }
      //If z is the second longest
      else{
	parts[2] = PartitionSide(zLen, xLen, numProcs, false);
	parts[0] = numProcs/parts[2];
      }
    }
    //If z is the longest side
    else{
      parts[2] = PartitionSide(zLen, xLen * yLen, numProcs, true);
      numProcs /= parts[2];
      //If x is the second longest
      if(xLen >= yLen){
	parts[0] = PartitionSide(xLen, yLen, numProcs, false);
	parts[1] = numProcs/parts[0];
      }
      //If y is the second longest
      else{
	parts[1] = PartitionSide(yLen, xLen, numProcs, false);
	parts[0] = numProcs/parts[1];
      }
    }
    
    //Divide up each side to give the lengths
    lengths[0] = xLen/(double)parts[0];
    lengths[1] = yLen/(double)parts[1];
    lengths[2] = zLen/(double)parts[2];
    return MB_SUCCESS;
  }
  //Partition a side based on the length ratios
  int ParallelMergeMesh::PartitionSide(double sideLen, double restLen, unsigned numProcs, bool altRatio)
  {
    //If theres only 1 processor, then just return 1
    if(numProcs == 1){
      return 1;
    }
    //Initialize with the ratio of 1 proc
    double ratio = -DBL_MAX;
    unsigned factor = 1;
    //We need to be able to save the last ratio and factor (for comparison)
    double oldRatio = ratio;
    double oldFactor = 1;
    
    //This is the ratio were shooting for
    double goalRatio = sideLen/restLen;

    //Calculate the divisor and numerator power
    //This avoid if statements in the loop and is useful since both calculations are similar
    double divisor, p;
    if(altRatio){
      divisor = (double)numProcs * sideLen;
      p = 3;
    }
    else{
      divisor = (double)numProcs;
      p = 2;
    }
    
    //Find each possible factor
    for (unsigned i = 2; i <= numProcs/2; i++){
      //If it is a factor...
      if (numProcs % i == 0){
	//We need to save the past factor
	oldRatio = ratio;
	oldFactor = factor;
	//There are 2 different ways to calculate the ratio:
	//Comparing 1 side to 2 sides: (i*i*i)/(numProcs*x)
	//Justification:  We have a ratio x:y:z (side Lengths) == a:b:c (procs).  So a=kx, b=ky, c=kz.
	//Also, abc=n (numProcs) => bc = n/a.  Also, a=kx => k=a/x => 1/k=x/a
	//And so x/(yz) == (kx)/(kyz) == (kx)/(kykz(1/k)) == a/(bc(x/a)) == a/((n/a)(x/a)) == a^3/(nx).
	//Comparing 1 side to 1 side: (i*i)/numprocs
	//Justification: i/(n/i) == i^2/n
	ratio = pow((double)i, p)/divisor;
	factor = i;
	//Once we have passed the goal ratio, we can break since we'll only move away from the goal ratio
	if(ratio >= goalRatio){
	  break;
	}
      }
    }
    //If we havent reached the goal ratio yet, check out factor = numProcs
    if(ratio < goalRatio){
      oldRatio = ratio;
      oldFactor = factor;
      factor = numProcs;
      ratio = pow((double)numProcs, p)/divisor;
    }
    
    //Figure out if our oldRatio is better than ratio
    if(fabs(ratio - goalRatio) > fabs(oldRatio-goalRatio)){
      factor = oldFactor;
    }
    //Return our findings
    return factor;
  }
  
  //Id the tuples that are matching
  ErrorCode ParallelMergeMesh::PopulateMyMatches()
  {
    //Counters for accessing tuples more efficiently
    unsigned long i = 0, mat_i=0, mat_ul=0, j=0, tup_r=0;
    double eps2 = myEps*myEps;
    while((i+1)<myTup.n){
      //Proximity Comparison
      double xi = myTup.vr[tup_r], 
	yi = myTup.vr[tup_r+1],
	zi = myTup.vr[tup_r+2];

      bool done = false;
      while(!done){
	j++; tup_r += myTup.mr;
	if(j >= myTup.n){
	  break;
	}
	CartVect cv(myTup.vr[tup_r]-xi,
		    myTup.vr[tup_r+1]-yi,
		    myTup.vr[tup_r+2]-zi);
	if(cv.length_squared() > eps2){
	  done = true;
	}
      }
      //Allocate the tuple list before adding matches
      while(myMatches.n+(j-i)*(j-i-1) >= myMatches.max){
	tuple_list_grow(&myMatches);
      }
 
      //We now know that tuples [i to j) exclusive match.  
      //If n tuples match, n*(n-1) match tuples will be made
      //tuples are of the form (proc1,proc2,handle1,handle2)
      if(i+1 < j){
	int kproc = i*myTup.mi;
	unsigned long khand = i*myTup.mul;
	for(unsigned long k = i; k<j; k++){
	  int lproc = kproc+myTup.mi;
	  unsigned long lhand = khand+myTup.mul;
	  for(unsigned long l=k+1; l<j; l++){
	    myMatches.vi[mat_i++]=myTup.vi[kproc];//proc1
	    myMatches.vi[mat_i++]=myTup.vi[lproc];//proc2
	    myMatches.vul[mat_ul++]=myTup.vul[khand];//handle1
	    myMatches.vul[mat_ul++]=myTup.vul[lhand];//handle2
	    myMatches.n++;
	    
	    myMatches.vi[mat_i++]=myTup.vi[lproc];//proc1
	    myMatches.vi[mat_i++]=myTup.vi[kproc];//proc2
	    myMatches.vul[mat_ul++]=myTup.vul[lhand];//handle1
	    myMatches.vul[mat_ul++]=myTup.vul[khand];//handle2
	    myMatches.n++;
	    lproc += myTup.mi;
	    lhand += myTup.mul;
	  }
	  kproc += myTup.mi;
	  khand += myTup.mul;
	}//End for(int k...
      }
      i = j;
    }//End while(i+1<tup.n)
    return MB_SUCCESS;
  }

  //Sort the matching tuples so that vertices can be tagged accurately
  ErrorCode ParallelMergeMesh::SortMyMatches()
  {
    buffer buf;
    unsigned long long max_size = mySkinEnts[0].size()*(MAX_SHARING_PROCS+1)*sizeof(double);
    buffer_init(&buf, max_size);
    //Sorts are necessary to check for doubles
    //Sort by remote handle
    moab_tuple_list_sort(&myMatches,3,&buf);
    //Sort by matching proc
    moab_tuple_list_sort(&myMatches,1,&buf);
    //Sort by local handle
    moab_tuple_list_sort(&myMatches,2,&buf);
    buffer_free(&buf);
    return MB_SUCCESS;
  }

  //Tag the shared elements using existing PComm functionality
  ErrorCode ParallelMergeMesh::TagSharedElements(int dim)
  {
    //Manipulate the matches list to tag vertices and entities
    //Set up proc ents
    Range proc_ents;
    ErrorCode rval;

    // get the entities in the partition sets
    for (Range::iterator rit = myPcomm->partitionSets.begin(); rit != myPcomm->partitionSets.end(); rit++) {
      Range tmp_ents;
      rval = myMB->get_entities_by_handle(*rit, tmp_ents, true);
      if (MB_SUCCESS != rval){
	return rval;
      }
      proc_ents.merge(tmp_ents);
    }
    if (myMB->dimension_from_handle(*proc_ents.rbegin()) !=
	myMB->dimension_from_handle(*proc_ents.begin())) {
      Range::iterator lower = proc_ents.lower_bound(CN::TypeDimensionMap[0].first),
	upper = proc_ents.upper_bound(CN::TypeDimensionMap[dim-1].second);
      proc_ents.erase(lower, upper);
    }
    

    //This vector doesnt appear to be used but its in resolve_shared_ents
    int maxp = -1;
    std::vector<int> sharing_procs(MAX_SHARING_PROCS);
    std::fill(sharing_procs.begin(), sharing_procs.end(), maxp);

    // get ents shared by 1 or n procs
    std::map<std::vector<int>, std::vector<EntityHandle> > proc_nranges;
    Range proc_verts;
    rval = myMB->get_adjacencies(proc_ents, 0, false, proc_verts,
				   Interface::UNION);
    if(rval != MB_SUCCESS){
      return rval;
    }

    rval = myPcomm->tag_shared_verts(myMatches, proc_nranges, proc_verts);
    if(rval != MB_SUCCESS){
      return rval;
    }
    
    // get entities shared by 1 or n procs
    rval = myPcomm->tag_shared_ents(dim,dim-1, &mySkinEnts[0],proc_nranges);
    if(rval != MB_SUCCESS){
      return rval;
    }
    
    // create the sets for each interface; store them as tags on
    // the interface instance
    Range iface_sets;
    rval = myPcomm->create_interface_sets(proc_nranges, dim, dim-1);
    if(rval != MB_SUCCESS){
      return rval;
    }
    // establish comm procs and buffers for them
    std::set<unsigned int> procs;
    rval = myPcomm->get_interface_procs(procs, true);
    if(rval != MB_SUCCESS){
      return rval;
    }

    // resolve shared entity remote handles; implemented in ghost cell exchange
    // code because it's so similar
    rval = myPcomm->exchange_ghost_cells(-1, -1, 0, true, true);
    if(rval != MB_SUCCESS){
      return rval;
    }
    // now build parent/child links for interface sets
    rval = myPcomm->create_iface_pc_links();
    return rval;
  }

  //Make sure to free up any allocated data
  //Need to avoid a double free
  void ParallelMergeMesh::CleanUp()
  {
    //Delete the matches tuple if necessary
    if(myMatches.max > 0){
      tuple_list_free(&myMatches);
    }
    //Delete the myTup if necessary
    if(myTup.max > 0){
      tuple_list_free(&myTup);
    }
    //Free up the crystal router
    if(cdAllocated){
      moab_crystal_free(&myCD);
    }
  }

  //Swap around tuples
  void ParallelMergeMesh::SwapTuples(tuple_list *tup, 
				     unsigned long a, unsigned long b)
  {
    if(a==b) return;
    //Swap mi
    unsigned long a_val = a*tup->mi, b_val=b*tup->mi;
    for(unsigned long i=0; i< tup->mi;i++){
      sint t =tup->vi[a_val];
      tup->vi[a_val] = tup->vi[b_val];
      tup->vi[b_val] = t; 
      a_val++;
      b_val++;
    }
    //Swap ml
    a_val = a*tup->ml;
    b_val = b*tup->ml;
    for(unsigned long i=0; i< tup->ml;i++){
      slong t =tup->vl[a_val];
      tup->vl[a_val] = tup->vl[b_val];
      tup->vl[b_val] = t;
      a_val++;
      b_val++;
    }
    //Swap mul
    a_val = a*tup->mul;
    b_val = b*tup->mul;
    for(unsigned long i=0; i< tup->mul;i++){
      ulong t =tup->vul[a_val];
      tup->vul[a_val] = tup->vul[b_val];
      tup->vul[b_val] = t; 
      a_val++;
      b_val++;
    }
    //Swap mr
    a_val = a*tup->mr;
    b_val = b*tup->mr;
    for(unsigned long i=0; i< tup->mr;i++){
      real t =tup->vr[a_val];
      tup->vr[a_val] = tup->vr[b_val];
      tup->vr[b_val] = t; 
      a_val++;
      b_val++;
    }
  }

  //Simple quick  sort to test real
  void ParallelMergeMesh::SortTuplesByReal(tuple_list *tup,
					  double eps)
  {
    //Call the recursive function
    PerformRealSort(tup, 0, tup->n,eps);
  }

  //Perform the sorting of a tuple by real
  //To sort an entire tuple_list, call (tup,0,tup.n.epsilon) 
  void ParallelMergeMesh::PerformRealSort(tuple_list *tup, 
					  unsigned long left, 
					  unsigned long right,
					  double eps)
  {  
    //If list size is only 1 or 0 return
    if(left+1 >= right){
      return;
    }
    unsigned long swap = left, tup_l = left*tup->mr, tup_t = tup_l + tup->mr;

    //Swap the median with the left position for a (hopefully) better split
    SwapTuples(tup,left,(left+right)/2);

    //Partition the data
    for(unsigned long t=left+1;t<right;t++){
      //If the left value(pivot) is greater than t_val, swap it into swap
      if(TupleGreaterThan(tup,tup_l,tup_t,eps)){
	swap++;
	SwapTuples(tup,swap,t);
      }
      tup_t+=tup->mr;
    }

    //Swap so that position swap is in the correct position
    SwapTuples(tup,left,swap);

    //Sort left and right of swap
    PerformRealSort(tup,left,swap,eps);
    PerformRealSort(tup,swap+1,right,eps);
  }

  //Note, this takes the actual tup->vr[] index (aka i*tup->mr)
  bool ParallelMergeMesh::TupleGreaterThan(tuple_list *tup, 
					   unsigned long vrI, 
					   unsigned long vrJ, 
					   double eps){
    unsigned check=0;
    while(check < tup->mr){
      //If the values are the same
      if(fabs(tup->vr[vrI+check]-tup->vr[vrJ+check]) <= eps){
	check++;
	continue;
      }
      //If the I greater than J 
      else if(tup->vr[vrI+check] > tup->vr[vrJ+check]){
	return true;
      }
      //If J greater than I
      else{
	return false;
      }
    }
    //Values are the same
    return false;
  }

}//End namespace moab
