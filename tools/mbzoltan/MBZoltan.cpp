/* $Id$
 *
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

/** Get a mesh from MOAB and write a Zoltan partition set back into MOAB and to
 *  a file
 *
 * 
 */

#include <iostream>
#include <assert.h>

#include "MBZoltan.hpp"
#include "moab/Interface.hpp"
#include "Internals.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/Range.hpp"
#include "moab/WriteUtilIface.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "moab/ParallelComm.hpp"
#include "MBTagConventions.hpp"
#include "moab/CN.hpp"

using namespace moab;

#define RR if (MB_SUCCESS != result) return result

static double *Points=NULL;
static int *GlobalIds=NULL;
static int NumPoints=0;
static int *NumEdges=NULL;
static int *NborGlobalId=NULL;
static int *NborProcs=NULL;

const bool debug = false;

MBZoltan::MBZoltan( Interface *impl , 
                    const bool use_coords,
                    int argc, 
                    char **argv ) 
                   : mbImpl(impl), 
                     myZZ(NULL), 
                     newMoab(false), 
                     useCoords(use_coords),
                     argcArg(argc), 
                     argvArg(argv)
{
  mbpc = ParallelComm::get_pcomm(mbImpl, 0);
  if (!mbpc)
    mbpc = new ParallelComm( impl, MPI_COMM_WORLD, 0 );
}

MBZoltan::~MBZoltan() 
{
  if (NULL == myZZ) delete myZZ;
}

ErrorCode MBZoltan::balance_mesh(const char *zmethod,
                                   const char *other_method,
                                   const bool write_as_sets,
                                   const bool write_as_tags) 
{
  if (!strcmp(zmethod, "RCB") && !strcmp(zmethod, "RIB") &&
      !strcmp(zmethod, "HSFC") && !strcmp(zmethod, "Hypergraph") &&
      !strcmp(zmethod, "PHG") && !strcmp(zmethod, "PARMETIS") &&
      !strcmp(zmethod, "OCTPART")) 
  {
    std::cout << "ERROR node " << mbpc->proc_config().proc_rank() << ": Method must be "
              << "RCB, RIB, HSFC, Hypergraph (PHG), PARMETIS, or OCTPART"
              << std::endl;
    return MB_FAILURE;
  }

  std::vector<double> pts; // x[0], y[0], z[0], ... from MOAB
  std::vector<int> ids; // point ids from MOAB
  std::vector<int> adjs, length;
  Range elems;

  // Get a mesh from MOAB and divide it across processors.

  ErrorCode result;
  
  if (mbpc->proc_config().proc_rank() == 0) {
    result = assemble_graph(3, pts, ids, adjs, length, elems); RR;
  }
  
  myNumPts = mbInitializePoints((int)ids.size(), &pts[0], &ids[0], &adjs[0],
                                &length[0]);

  // Initialize Zoltan.  This is a C call.  The simple C++ code 
  // that creates Zoltan objects does not keep track of whether 
  // Zoltan_Initialize has been called.

  float version;

  Zoltan_Initialize(argcArg, argvArg, &version); 

  // Create Zoltan object.  This calls Zoltan_Create.  
  if (NULL == myZZ) myZZ = new Zoltan(MPI_COMM_WORLD);

  if (NULL == zmethod || !strcmp(zmethod, "RCB"))
    SetRCB_Parameters();
  else if (!strcmp(zmethod, "RIB"))
    SetRIB_Parameters();
  else if (!strcmp(zmethod, "HSFC"))
    SetHSFC_Parameters();
  else if (!strcmp(zmethod, "Hypergraph") || !strcmp(zmethod, "PHG"))
    if (NULL == other_method)
      SetHypergraph_Parameters("auto");
    else
      SetHypergraph_Parameters(other_method);
  else if (!strcmp(zmethod, "PARMETIS")) {
    if (NULL == other_method)
      SetPARMETIS_Parameters("RepartGDiffusion");
    else
      SetPARMETIS_Parameters(other_method);
  }
  else if (!strcmp(zmethod, "OCTPART")) {
    if (NULL == other_method)
      SetOCTPART_Parameters("2");
    else
      SetOCTPART_Parameters(other_method);
  }
  
  // Call backs:

  myZZ->Set_Num_Obj_Fn(mbGetNumberOfAssignedObjects, NULL);
  myZZ->Set_Obj_List_Fn(mbGetObjectList, NULL);
  myZZ->Set_Num_Geom_Fn(mbGetObjectSize, NULL);
  myZZ->Set_Geom_Multi_Fn(mbGetObject, NULL);
  myZZ->Set_Num_Edges_Multi_Fn(mbGetNumberOfEdges, NULL);
  myZZ->Set_Edge_List_Multi_Fn(mbGetEdgeList, NULL);

  // Perform the load balancing partitioning

  int changes;
  int numGidEntries;
  int numLidEntries;
  int numImport;
  ZOLTAN_ID_PTR importGlobalIds;
  ZOLTAN_ID_PTR importLocalIds;
  int *importProcs;
  int *importToPart;
  int numExport;
  ZOLTAN_ID_PTR exportGlobalIds;
  ZOLTAN_ID_PTR exportLocalIds;
  int *exportProcs;
  int *exportToPart;

  int rc = myZZ->LB_Partition(changes, numGidEntries, numLidEntries, 
                            numImport, importGlobalIds, importLocalIds, 
                            importProcs, importToPart,
                            numExport, exportGlobalIds, exportLocalIds, 
                            exportProcs, exportToPart);

  rc = mbGlobalSuccess(rc);
  
  if (!rc){
    mbPrintGlobalResult("==============Result==============", 
                        myNumPts, numImport, numExport, changes);
  }
  else{
    return MB_FAILURE;
  }
  
  // take results & write onto MOAB partition sets

  int *assignment;

  mbFinalizePoints((int)ids.size(), numExport, exportLocalIds,
                   exportProcs, &assignment);
  
  if (mbpc->proc_config().proc_rank() == 0) {
    ErrorCode result = write_partition(mbpc->proc_config().proc_size(), elems, assignment,
                                         write_as_sets, write_as_tags);

    if (MB_SUCCESS != result) return result;

    free((int *) assignment);
  }

  // Free the memory allocated for lists returned by LB_Parition()

  myZZ->LB_Free_Part(&importGlobalIds, &importLocalIds, &importProcs, 
                     &importToPart);
  myZZ->LB_Free_Part(&exportGlobalIds, &exportLocalIds, &exportProcs, 
                     &exportToPart);

  // Implementation note:  A Zoltan object contains an MPI communicator.
  //   When the Zoltan object is destroyed, it uses it's MPI communicator.
  //   So it is important that the Zoltan object is destroyed before
  //   the MPI communicator is destroyed.  To ensure this, dynamically
  //   allocate the Zoltan object, so you can explicitly destroy it.
  //   If you create a Zoltan object on the stack, it's destructor will
  //   be invoked atexit, possibly after the communicator's
  //   destructor.

  return MB_SUCCESS;
}

ErrorCode MBZoltan::partition_mesh(const int nparts,
                                     const char *zmethod,
                                     const char *other_method,
                                     const bool write_as_sets,
                                     const bool write_as_tags,
                                     const int part_dim) 
{
    // should only be called in serial
  if (mbpc->proc_config().proc_size() != 1) {
    std::cout << "MBZoltan::partition_mesh must be called in serial." 
              << std::endl;
    return MB_FAILURE;
  }
  
  if (NULL != zmethod && strcmp(zmethod, "RCB") && strcmp(zmethod, "RIB") &&
      strcmp(zmethod, "HSFC") && strcmp(zmethod, "Hypergraph") &&
      strcmp(zmethod, "PHG") && strcmp(zmethod, "PARMETIS") &&
      strcmp(zmethod, "OCTPART")) 
  {
    std::cout << "ERROR node " << mbpc->proc_config().proc_rank() << ": Method must be "
              << "RCB, RIB, HSFC, Hypergraph (PHG), PARMETIS, or OCTPART"
              << std::endl;
    return MB_FAILURE;
  }

  std::vector<double> pts; // x[0], y[0], z[0], ... from MOAB
  std::vector<int> ids; // point ids from MOAB
  std::vector<int> adjs, length;
  Range elems;

  // Get a mesh from MOAB and divide it across processors.

  ErrorCode result;
  
  result = assemble_graph(part_dim, pts, ids, adjs, length, elems); RR;
  
  myNumPts = mbInitializePoints((int)ids.size(), &pts[0], &ids[0], &adjs[0],
                                &length[0]);

  // Initialize Zoltan.  This is a C call.  The simple C++ code 
  // that creates Zoltan objects does not keep track of whether 
  // Zoltan_Initialize has been called.

  float version;

  Zoltan_Initialize(argcArg, argvArg, &version); 

  // Create Zoltan object.  This calls Zoltan_Create.  
  if (NULL == myZZ) myZZ = new Zoltan(MPI_COMM_WORLD);

  if (NULL == zmethod || !strcmp(zmethod, "RCB"))
    SetRCB_Parameters();
  else if (!strcmp(zmethod, "RIB"))
    SetRIB_Parameters();
  else if (!strcmp(zmethod, "HSFC"))
    SetHSFC_Parameters();
  else if (!strcmp(zmethod, "Hypergraph") || !strcmp(zmethod, "PHG")) {
    if (NULL == other_method)
      SetHypergraph_Parameters("auto");
    else
      SetHypergraph_Parameters(other_method);
  }
  else if (!strcmp(zmethod, "PARMETIS")) {
    if (NULL == other_method)
      SetPARMETIS_Parameters("RepartGDiffusion");
    else
      SetPARMETIS_Parameters(other_method);
  }
  else if (!strcmp(zmethod, "OCTPART")) {
    if (NULL == other_method)
      SetOCTPART_Parameters("2");
    else
      SetOCTPART_Parameters(other_method);
  }
  
    // set # requested partitions
  char buff[10];
  sprintf(buff, "%d", nparts);
  int retval = myZZ->Set_Param("NUM_GLOBAL_PARTITIONS", buff);
  if (ZOLTAN_OK != retval) return MB_FAILURE;

    // request only partition assignments
  retval = myZZ->Set_Param("RETURN_LISTS", "PARTITION ASSIGNMENTS");
  if (ZOLTAN_OK != retval) return MB_FAILURE;
  
  // Call backs:

  myZZ->Set_Num_Obj_Fn(mbGetNumberOfAssignedObjects, NULL);
  myZZ->Set_Obj_List_Fn(mbGetObjectList, NULL);
  myZZ->Set_Num_Geom_Fn(mbGetObjectSize, NULL);
  myZZ->Set_Geom_Multi_Fn(mbGetObject, NULL);
  myZZ->Set_Num_Edges_Multi_Fn(mbGetNumberOfEdges, NULL);
  myZZ->Set_Edge_List_Multi_Fn(mbGetEdgeList, NULL);

  // Perform the load balancing partitioning

  int changes;
  int numGidEntries;
  int numLidEntries;
  int dumnum1;
  ZOLTAN_ID_PTR dum_local, dum_global;
  int *dum1, *dum2;
  int num_assign;
  ZOLTAN_ID_PTR assign_gid, assign_lid;
  int *assign_procs, *assign_parts;

  retval = myZZ->LB_Partition(changes, numGidEntries, numLidEntries, 
                              dumnum1, dum_global, dum_local, dum1, dum2,
                              num_assign, assign_gid, assign_lid,
                              assign_procs, assign_parts);
  if (ZOLTAN_OK != retval) return MB_FAILURE;
  
  // take results & write onto MOAB partition sets
  result = write_partition(nparts, elems, assign_parts,
                           write_as_sets, write_as_tags);

  if (MB_SUCCESS != result) return result;


  // Free the memory allocated for lists returned by LB_Parition()
  myZZ->LB_Free_Part(&assign_gid, &assign_lid, &assign_procs, &assign_parts);

  return MB_SUCCESS;
}

ErrorCode MBZoltan::assemble_graph(const int dimension,
                                     std::vector<double> &coords,
                                     std::vector<int> &moab_ids,
                                     std::vector<int> &adjacencies, 
                                     std::vector<int> &length,
                                     Range &elems) 
{
    // assemble a graph with vertices equal to elements of specified dimension, edges
    // signified by list of other elements to which an element is connected

    // get the elements of that dimension
  ErrorCode result = mbImpl->get_entities_by_dimension(0, dimension, elems);
  if (MB_SUCCESS != result || elems.empty()) return result;
  
    // assign global ids
  result = mbpc->assign_global_ids(0, dimension); RR;

    // now assemble the graph, calling MeshTopoUtil to get bridge adjacencies through d-1 dimensional
    // neighbors
  MeshTopoUtil mtu(mbImpl);
  Range adjs;
    // can use a fixed-size array 'cuz the number of lower-dimensional neighbors is limited
    // by MBCN
  int neighbors[5*MAX_SUB_ENTITIES];
  double avg_position[3];
  int moab_id;
  
    // get the global id tag hanlde
  Tag gid;
  result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, 4, MB_TAG_DENSE, MB_TYPE_INTEGER,
                              gid, NULL, true);
  if (MB_SUCCESS != result) return result;
  
  for (Range::iterator rit = elems.begin(); rit != elems.end(); rit++) {

      // get bridge adjacencies
    adjs.clear();
    result = mtu.get_bridge_adjacencies(*rit, (dimension > 0 ? dimension-1 : 3), 
                                        dimension, adjs); RR;
    
    
      // get the graph vertex ids of those
    if (!adjs.empty()) {
      assert(adjs.size() < 5*MAX_SUB_ENTITIES);
      result = mbImpl->tag_get_data(gid, adjs, neighbors); RR;
    }

      // copy those into adjacencies vector
    length.push_back((int)adjs.size());
    std::copy(neighbors, neighbors+adjs.size(), std::back_inserter(adjacencies));


      // get average position of vertices
    result = mtu.get_average_position(*rit, avg_position); RR;
    
      // get the graph vertex id for this element
    result = mbImpl->tag_get_data(gid, &(*rit), 1, &moab_id); RR;

      // copy those into coords vector
    moab_ids.push_back(moab_id);
    std::copy(avg_position, avg_position+3, std::back_inserter(coords));
  }

  if (debug) {
    std::cout << "Length vector: " << std::endl;
    std::copy(length.begin(), length.end(), std::ostream_iterator<int>(std::cout, ", "));
    std::cout << std::endl;
    std::cout << "Adjacencies vector: " << std::endl;
    std::copy(adjacencies.begin(), adjacencies.end(), std::ostream_iterator<int>(std::cout, ", "));
    std::cout << std::endl;
    std::cout << "Moab_ids vector: " << std::endl;
    std::copy(moab_ids.begin(), moab_ids.end(), std::ostream_iterator<int>(std::cout, ", "));
    std::cout << std::endl;
    std::cout << "Coords vector: " << std::endl;
    std::copy(coords.begin(), coords.end(), std::ostream_iterator<double>(std::cout, ", "));
    std::cout << std::endl;
  }

  return MB_SUCCESS;
}

ErrorCode MBZoltan::write_partition(const int nparts,
                                      Range &elems, 
                                      const int *assignment,
                                      const bool write_as_sets,
                                      const bool write_as_tags) 
{
  ErrorCode result;

    // get the partition set tag
  Tag part_set_tag;
  int dum_id = -1;
  result = mbImpl->tag_create("PARALLEL_PARTITION", 4, MB_TAG_SPARSE, MB_TYPE_INTEGER, part_set_tag, &dum_id,
                              true); RR;
  
    // get any sets already with this tag, and clear them
  Range tagged_sets;
  result = mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &part_set_tag, NULL, 1,
                                                tagged_sets, Interface::UNION); RR;
  if (!tagged_sets.empty()) {
    result = mbImpl->clear_meshset(tagged_sets); RR;
    if (!write_as_sets) {
      result = mbImpl->tag_delete_data(part_set_tag, tagged_sets); RR;
    }
  }

  if (write_as_sets) {
      // first, create partition sets and store in vector
    EntityHandle *part_sets = new EntityHandle[nparts];
    if (NULL == part_sets) return MB_FAILURE;
  
    if (nparts > (int) tagged_sets.size()) {
        // too few partition sets - create missing ones
      int num_new = nparts - tagged_sets.size();
      for (int i = 0; i < num_new; i++) {
        EntityHandle new_set;
        result = mbImpl->create_meshset(MESHSET_SET, new_set); RR;
        tagged_sets.insert(new_set);
      }
    }
    else if (nparts < (int) tagged_sets.size()) {
        // too many partition sets - delete extras
      int num_del = tagged_sets.size() - nparts;
      for (int i = 0; i < num_del; i++) {
        EntityHandle old_set = tagged_sets.pop_back();
        result = mbImpl->delete_entities(&old_set, 1); RR;
      }
    }
  
      // assign partition sets to vector
    Range::iterator rit = tagged_sets.begin();
    int i = 0;
    for (; i < nparts; rit++, i++) part_sets[i] = *rit;
  
      // write a tag to those sets denoting they're partition sets, with a value of the
      // proc number
    int *dum_ids = new int[nparts];
    for (int i = 0; i < nparts; i++) dum_ids[i] = i;
  
    result = mbImpl->tag_set_data(part_set_tag, part_sets, nparts, dum_ids); RR;

      // assign entities to the relevant sets
    for (rit = elems.begin(), i = 0; rit != elems.end(); rit++, i++) {
      result = mbImpl->add_entities(part_sets[assignment[i]], &(*rit), 1); RR;
    }

      // check for empty sets, warn if there are any
    std::vector<unsigned int> empty_sets;
    for (i = 0; i < nparts; i++) {
      int num_ents = 0;
      result = mbImpl->get_number_entities_by_handle(part_sets[i], num_ents);
      if (MB_SUCCESS != result || !num_ents) empty_sets.push_back(i);
    }
    if (!empty_sets.empty()) {
      std::cout << "WARNING: " << empty_sets.size() << " empty sets in partition: ";
      for (std::vector<unsigned int>::iterator vit = empty_sets.begin(); vit != empty_sets.end(); vit++)
        std::cout << *vit << " ";
      std::cout << std::endl;
    }
  }
  
  if (write_as_tags) {
      // allocate integer-size partitions
    result = mbImpl->tag_set_data(part_set_tag, elems, assignment); RR;
  }
  
  return MB_SUCCESS;
}

void MBZoltan::SetRCB_Parameters()
{
  if (mbpc->proc_config().proc_rank() == 0) std::cout << "\nRecursive Coordinate Bisection" << std::endl;
  // General parameters:

  myZZ->Set_Param("DEBUG_LEVEL", "0");     // no debug messages
  myZZ->Set_Param("LB_METHOD", "RCB");     // recursive coordinate bisection

  // RCB parameters:

  myZZ->Set_Param("RCB_OUTPUT_LEVEL", "2");
  myZZ->Set_Param("KEEP_CUTS", "1");              // save decomposition
  myZZ->Set_Param("RCB_RECTILINEAR_BLOCKS", "1"); // don't split point on boundary
}

void MBZoltan::SetRIB_Parameters()
{
  if (mbpc->proc_config().proc_rank() == 0) std::cout << "\nRecursive Inertial Bisection" << std::endl;
  // General parameters:

  myZZ->Set_Param("DEBUG_LEVEL", "0");     // no debug messages
  myZZ->Set_Param("LB_METHOD", "RIB");     // Recursive Inertial Bisection

  // RIB parameters:

  myZZ->Set_Param("KEEP_CUTS", "1");              // save decomposition
  myZZ->Set_Param("AVERAGE_CUTS", "1");
}

void MBZoltan::SetHSFC_Parameters()
{
  if (mbpc->proc_config().proc_rank() == 0) std::cout << "\nHilbert Space Filling Curve" << std::endl;
  // General parameters:

  myZZ->Set_Param("DEBUG_LEVEL", "0");     // no debug messages
  myZZ->Set_Param("LB_METHOD", "HSFC");    // perform Hilbert space filling curve

  // HSFC parameters:

  myZZ->Set_Param("KEEP_CUTS", "1");              // save decomposition
}

void MBZoltan::SetHypergraph_Parameters(const char *phg_method)
{
  if (mbpc->proc_config().proc_rank() == 0) std::cout << "\nHypergraph (or PHG): " << std::endl;
  // General parameters:

  myZZ->Set_Param("DEBUG_LEVEL", "0");     // no debug messages
  myZZ->Set_Param("LB_METHOD", "Hypergraph");     // Hypergraph (or PHG)

  // Hypergraph (or PHG) parameters:
  myZZ->Set_Param("PHG_COARSEPARTITION_METHOD",phg_method);//CoarsePartitionMethod
}

void MBZoltan::SetPARMETIS_Parameters(const char *parmetis_method)
{
  if (mbpc->proc_config().proc_rank() == 0) std::cout << "\nPARMETIS: " << parmetis_method << std::endl;
  // General parameters:

  myZZ->Set_Param("DEBUG_LEVEL", "0");     // no debug messages
  myZZ->Set_Param("LB_METHOD", "PARMETIS");     // the ParMETIS library

  // PARMETIS parameters:

  myZZ->Set_Param("PARMETIS_METHOD", parmetis_method); // method in the library
}

void MBZoltan::SetOCTPART_Parameters(const char *oct_method)
{
  if (mbpc->proc_config().proc_rank() == 0) std::cout << "\nOctree Partitioning: " << oct_method
			   << std::endl;
  // General parameters:

  myZZ->Set_Param("DEBUG_LEVEL", "0");     // no debug messages
  myZZ->Set_Param("LB_METHOD", "OCTPART");     // octree partitioning

  // OCTPART parameters:

  myZZ->Set_Param("OCT_METHOD", oct_method); // the SFC to be used
  myZZ->Set_Param("OCT_OUTPUT_LEVEL", "3");
}

int MBZoltan::mbInitializePoints(int npts, double *pts, int *ids, 
                                 int *adjs, int *length)
{
  unsigned int i;
  int j;
  int *numPts, *nborProcs;
  int sum, ptsPerProc, ptsAssigned, mySize;
  MPI_Status stat;
  double *sendPts;
  int *sendIds;
  int *sendEdges;
  int *sendNborId;
  int *sendProcs;

  if (mbpc->proc_config().proc_rank() == 0)
  {
      /* divide pts to start */

    MALLOC(numPts, sizeof(int) * mbpc->proc_config().proc_size(), int *);
    ptsPerProc = npts / mbpc->proc_config().proc_size();
    ptsAssigned = 0;

    for (i=0; i < mbpc->proc_config().proc_size()-1; i++)
    {
      numPts[i] = ptsPerProc;
      ptsAssigned += ptsPerProc;
    }

    numPts[mbpc->proc_config().proc_size()-1] = npts - ptsAssigned;

    mySize = numPts[mbpc->proc_config().proc_rank()];
    sendPts = pts + (3 * numPts[0]);
    sendIds = ids + numPts[0];
    sendEdges = length + numPts[0];
    sum = 0;

    for (j=0; j<numPts[0]; j++)
      sum += length[j];

    sendNborId = adjs + sum;

    for (j=numPts[0]; j<npts; j++)
      sum += length[j];

    MALLOC(nborProcs, sizeof(int) * sum, int *);

    for (j=0; j<sum; j++)
      if ((i = adjs[j]/ptsPerProc) < mbpc->proc_config().proc_size())
        nborProcs[j] = i;
      else
        nborProcs[j] = mbpc->proc_config().proc_size() - 1;

    sendProcs = nborProcs + (sendNborId - adjs);

    for (i=1; i<mbpc->proc_config().proc_size(); i++)
    {
      MPI_Send(&numPts[i], 1, MPI_INT, i, 0x00,MPI_COMM_WORLD);
      MPI_Send(sendPts, 3 * numPts[i], MPI_DOUBLE, i, 0x01,MPI_COMM_WORLD);
      MPI_Send(sendIds, numPts[i], MPI_INT, i, 0x03,MPI_COMM_WORLD);
      MPI_Send(sendEdges, numPts[i], MPI_INT, i, 0x06,MPI_COMM_WORLD);
      sum = 0;

      for (j=0; j<numPts[i]; j++)
        sum += sendEdges[j];

      MPI_Send(sendNborId, sum, MPI_INT, i, 0x07,MPI_COMM_WORLD);
      MPI_Send(sendProcs, sum, MPI_INT, i, 0x08,MPI_COMM_WORLD);
      sendPts += (3 * numPts[i]);
      sendIds += numPts[i];
      sendEdges += numPts[i];
      sendNborId += sum;
      sendProcs += sum;
    }

    free(numPts);
  }
  else
  {
    MPI_Recv(&mySize, 1, MPI_INT, 0, 0x00, MPI_COMM_WORLD, &stat);
    MALLOC(pts, sizeof(double) * 3 * mySize, double *);
    MALLOC(ids, sizeof(int) * mySize, int *);
    MALLOC(length, sizeof(int) * mySize, int *);
    MPI_Recv(pts, 3 * mySize, MPI_DOUBLE, 0, 0x01, MPI_COMM_WORLD, &stat);
    MPI_Recv(ids, mySize, MPI_INT, 0, 0x03, MPI_COMM_WORLD, &stat);
    MPI_Recv(length, mySize, MPI_INT, 0, 0x06, MPI_COMM_WORLD, &stat);
    sum = 0;

    for (j=0; j<mySize; j++)
      sum += length[j];

    MALLOC(adjs, sizeof(int) * sum, int *);
    MALLOC(nborProcs, sizeof(int) * sum, int *);
    MPI_Recv(adjs, sum, MPI_INT, 0, 0x07, MPI_COMM_WORLD, &stat);
    MPI_Recv(nborProcs, sum, MPI_INT, 0, 0x08, MPI_COMM_WORLD, &stat);
  }     
          
  Points = pts;
  GlobalIds = ids;  
  NumPoints = mySize;
  NumEdges = length;
  NborGlobalId = adjs;
  NborProcs = nborProcs;
          
  return mySize;
}     

void MBZoltan::mbFinalizePoints(int npts, int numExport,
                                ZOLTAN_ID_PTR exportLocalIDs, int *exportProcs,
                                int **assignment)
{
  int *MyAssignment;
  int i;
  int numPts;
  MPI_Status stat;
  int *recvA;

  /* assign pts to start */

  if (mbpc->proc_config().proc_rank() == 0) {
    MALLOC(MyAssignment, sizeof(int) * npts, int *);
  }
  else {
    MALLOC(MyAssignment, sizeof(int) * NumPoints, int *);
  }

  for (i=0; i<NumPoints; i++)
    MyAssignment[i] = mbpc->proc_config().proc_rank();

  for (i=0; i<numExport; i++)
    MyAssignment[exportLocalIDs[i]] = exportProcs[i];

  if (mbpc->proc_config().proc_rank() == 0)
    {
      /* collect pts */

      recvA = MyAssignment + NumPoints;

      for (i=1; i< (int) mbpc->proc_config().proc_size(); i++)
	{
	  MPI_Recv(&numPts, 1, MPI_INT, i, 0x04, MPI_COMM_WORLD, &stat);
	  MPI_Recv(recvA, numPts, MPI_INT, i, 0x05, MPI_COMM_WORLD, &stat);
	  recvA += numPts;
	}

      *assignment = MyAssignment;
    }
  else
    {
      MPI_Send(&NumPoints, 1, MPI_INT, 0, 0x04,MPI_COMM_WORLD);
      MPI_Send(MyAssignment, NumPoints, MPI_INT, 0, 0x05,MPI_COMM_WORLD);
      free(MyAssignment);
    }     
}

int MBZoltan::mbGlobalSuccess(int rc)
{
  int fail = 0;
  unsigned int i;
  int *vals;
  MALLOC(vals, mbpc->proc_config().proc_size() * sizeof(int), int *);

  MPI_Allgather(&rc, 1, MPI_INT, vals, 1, MPI_INT, MPI_COMM_WORLD);

  for (i=0; i<mbpc->proc_config().proc_size(); i++){
    if (vals[i] != ZOLTAN_OK){
      if (0 == mbpc->proc_config().proc_rank()){
        mbShowError(vals[i], "Result on process ");
      }
      fail = 1;
    }
  }

  free(vals);
  return fail;
}

void MBZoltan::mbPrintGlobalResult(const char *s, 
                                   int begin, int import, int exp, int change)
{
  unsigned int i;
  int *v1;
  MALLOC(v1, 4 * sizeof(int), int *);
  int *v2 = NULL;
  int *v;

  v1[0] = begin;
  v1[1] = import;
  v1[2] = exp;
  v1[3] = change;

  if (mbpc->proc_config().proc_rank() == 0){
    MALLOC(v2, 4 * mbpc->proc_config().proc_size() * sizeof(int), int*);
  }

  MPI_Gather(v1, 4, MPI_INT, v2, 4, MPI_INT, 0, MPI_COMM_WORLD);

  if (mbpc->proc_config().proc_rank() == 0){
    fprintf(stdout,"======%s======\n",s);
    for (i=0, v=v2; i<mbpc->proc_config().proc_size(); i++, v+=4){
      fprintf(stdout,"%d: originally had %d, import %d, exp %d, %s\n",
        i, v[0], v[1], v[2],
        v[3] ? "a change of partitioning" : "no change");
    }
    fprintf(stdout,"==========================================\n");
    fflush(stdout); 

    free(v2);
  }

  free(v1);
}

void MBZoltan::mbShowError(int val, const char *s)
{
  if (s)
    {
    printf("%s ",s);
    }
  switch (val)
    {
    case ZOLTAN_OK:
      printf("%d: SUCCESSFUL\n", mbpc->proc_config().proc_rank());
      break;
    case ZOLTAN_WARN:
      printf("%d: WARNING\n", mbpc->proc_config().proc_rank());
      break;
    case ZOLTAN_FATAL:
      printf("%d: FATAL ERROR\n", mbpc->proc_config().proc_rank());
      break;
    case ZOLTAN_MEMERR:
      printf("%d: MEMORY ALLOCATION ERROR\n", mbpc->proc_config().proc_rank());
      break;
    default:
      printf("%d: INVALID RETURN CODE\n", mbpc->proc_config().proc_rank());
      break;
    }
  return;
}

/**********************
** call backs
**********************/

int mbGetNumberOfAssignedObjects(void *userDefinedData, int *err)
{
  *err = 0;
  return NumPoints;
}

void mbGetObjectList(void *userDefinedData, int numGlobalIds, int numLids,
  ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int wgt_dim, float *obj_wgts,
  int *err)
{
  int i;
    
  for (i=0; i<NumPoints; i++)
    {
    gids[i] = GlobalIds[i];
    lids[i] = i;
    }
    
  *err = 0;
    
  return;
}

int mbGetObjectSize(void *userDefinedData, int *err)
{
  *err = 0; 
  return 3;
} 

void mbGetObject(void *userDefinedData, int numGlobalIds, int numLids, int numObjs,
  ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int numDim, double *pts, int *err)
{ 
  int i, id, id3;
  int next = 0;
  
  if (numDim != 3)    
    {
    *err = 1;         
    return;
    }

  for (i=0; i<numObjs; i++)
    {
    id = lids[i];
  
    if ((id < 0) || (id >= NumPoints))
      {
      *err = 1;
      return;
      }

    id3 = lids[i] * 3;

    pts[next++] = Points[id3];
    pts[next++] = Points[id3 + 1];
    pts[next++] = Points[id3 + 2];
    }
} 

void mbGetNumberOfEdges(void *userDefinedData, int numGlobalIds, int numLids,
			int numObjs, 
			ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,	int *numEdges,
			int *err)
{
  int i, id;
  int next = 0;

  for (i=0; i<numObjs; i++)
    {
      id = lids[i];

      if ((id < 0) || (id >= NumPoints))
	{
	  *err = 1;
	  return;
	}

      numEdges[next++] = NumEdges[id];
    }
}

void mbGetEdgeList(void *userDefinedData, int numGlobalIds, int numLids,
		   int numObjs,
		   ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int *numEdges,
		   ZOLTAN_ID_PTR nborGlobalIds, int *nborProcs, int wgt_dim,
		   float *edge_wgts, int *err)
{
  int i, id, idSum, j;
  int next = 0;

  for (i=0; i<numObjs; i++)
    {
      id = lids[i];

      if ((id < 0) || (id >= NumPoints))
	{
	  *err = 1;
	  return;
	}

      idSum = 0;

      for (j=0; j<id; j++)
	  idSum += NumEdges[j];

      for (j=0; j<NumEdges[id]; j++)
	{
	  nborGlobalIds[next] = NborGlobalId[idSum];
	  nborProcs[next++] = NborProcs[idSum++];
	}
    }
}

