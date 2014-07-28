/*
 * wrap_intx.cpp
 *  Will implement the intersection method that will be callable from fortran too
 *  will be added to the library mbcslam.a
 *
 *
 *  Created on: Dec 14, 2013
 *      Author: iulian
 */
#include "iMesh.h"
#include "iMeshP.h"
#include "MBiMesh.hpp"
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "Intx2MeshOnSphere.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/ParallelComm.hpp"
#include "MBTagConventions.hpp"
#include <mpi.h>

using namespace moab;
double radius = 1.;
double gtol = 1.e-9;
bool debug = true;

#define  NC  3

#ifdef __cplusplus
extern "C" {
#endif

void update_tracer(iMesh_Instance instance,
    iBase_EntitySetHandle imesh_euler_set, int * ierr) {
  Range ents;
  moab::Interface * mb = MOABI;
  *ierr = 1;

  EntityHandle euler_set = (EntityHandle) imesh_euler_set;

  Intx2MeshOnSphere worker(mb);
  worker.SetRadius(radius);

  worker.SetErrorTolerance(gtol);

  EntityHandle covering_lagr_set;

  ErrorCode rval = mb->create_meshset(MESHSET_SET, covering_lagr_set);
  ERRORV(rval, "can't create covering set ");

  // we need to update the correlation tag and remote tuples
  rval = worker.create_departure_mesh_2nd_alg(euler_set, covering_lagr_set);
  ERRORV(rval, "can't populate covering set ");

  if (debug) {
    rval = mb->write_file("lagr.h5m", 0, 0, &covering_lagr_set, 1);
    ERRORV(rval, "can't write covering set ");
  }

  //
  rval = enforce_convexity(mb, covering_lagr_set);
  ERRORV(rval, "can't write covering set ");

  EntityHandle outputSet;
  rval = mb->create_meshset(MESHSET_SET, outputSet);
  ERRORV(rval, "can't create output set ");

  rval = worker.intersect_meshes(covering_lagr_set, euler_set, outputSet);
  ERRORV(rval, "can't intersect ");

  if (debug) {
    rval = mb->write_file("output.vtk", 0, 0, &outputSet, 1);
    ERRORV(rval, "can't write covering set ");
  }

  // tagElem is the average computed at each element, from nodal values
  Tag tagElem = 0;
  std::string tag_name2("TracerAverage");
  rval = mb->tag_get_handle(tag_name2.c_str(), 1, MB_TYPE_DOUBLE, tagElem,
      MB_TAG_DENSE | MB_TAG_CREAT);
  ERRORV(rval, "can't get tracer tag ");

  // area of the euler element is fixed, store it; it is used to recompute the averages at each
  // time step
  Tag tagArea = 0;
  std::string tag_name4("Area");
  rval = mb->tag_get_handle(tag_name4.c_str(), 1, MB_TYPE_DOUBLE, tagArea,
      MB_TAG_DENSE | MB_TAG_CREAT);
  ERRORV(rval, "can't get area tag");

  rval = worker.update_tracer_data(outputSet, tagElem, tagArea);
  ERRORV(rval, "can't update tracer ");

  // everything can be deleted now from intx data; polygons, etc.

  *ierr = 0;
  return;
}

ErrorCode create_coarse_mesh(Interface * mb, ParallelComm * pcomm,
    EntityHandle coarseSet, double * coords, int * corners, int nc, int nelem) {

  int rank = pcomm->proc_config().proc_rank();

  Tag tagid; // global id at corners
  ErrorCode rval = mb->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER,
      tagid, MB_TAG_DENSE | MB_TAG_CREAT);
  ERRORR(rval, "can't create corner tag ");

  int dum_id = -1;
  Tag partitionTag;
  mb->tag_get_handle(PARALLEL_PARTITION_TAG_NAME, 1, MB_TYPE_INTEGER,
      partitionTag, MB_TAG_SPARSE | MB_TAG_CREAT, &dum_id);

  // there are nelem*4 corners, and 3*(nc2+1)*(nc2+1)*nelem coordinates
  // create first a coarse mesh,
  int size_corners = 4 * nelem;
  //int size_coords=3*(nc+1)*(nc+1)*nelem;
  // order first the corners array, and eliminate duplicates
  std::vector<int> corn1(size_corners);
  std::copy(corners, corners + size_corners, corn1.begin());
  std::sort(corn1.begin(), corn1.end());
  corn1.erase(std::unique(corn1.begin(), corn1.end()), corn1.end());

  int num_nodes_coarse = (int) corn1.size();

  std::map<int, int> index_map;
  for (size_t i = 0; i < corn1.size(); i++)
    index_map[corn1[i]] = i;

  std::vector<int> index_used(corn1.size(), 0);

  ReadUtilIface *read_iface;
  rval = mb->query_interface(read_iface);
  ERRORR(rval, "can't get query intf ");

  std::vector<double *> coordv;
  EntityHandle start_vert, start_elem, *connect;
// create verts, num is 2(nquads+1) because they're in a 1d row; will initialize coords in loop over quads later
  rval = read_iface->get_node_coords(3, num_nodes_coarse, 0, start_vert,
      coordv);
  ERRORR(rval, "can't get node coords ");
// fill it up
  int stride = (nc + 1) * (nc + 1); // 16, for nc ==3
  //int order_in_coord[] = {0, 3, 15, 12}; // for nc=3
  /*

   *  first j, then i, so this is the order of the points in coords array, now:
   *
   *   nc(nc+1),   ...      (nc+1)*(nc+1)-1
   *
   *   2(nc+1),
   *   nc+1,   nc+2,            2*(nc+1)-1
   *   0,      1 ,    2, ..........., nc
   */
  int order_in_coord[] = { 0, nc, (nc + 1) * (nc + 1) - 1, nc * (nc + 1) };
  for (int i = 0; i < size_corners; i++) {
    int id = corners[i];
    int index = index_map[id];

    int j = i % 4;
    int ind_coord = i / 4;
    ind_coord = (ind_coord * stride + order_in_coord[j]) * 3;
    if (index_used[index]) {
      if (fabs(coordv[0][index] - coords[ind_coord]) > 0.000001)
        std::cout << " id:" << corners[i] << " i:" << i << " j:" << j << " "
            << coordv[0][index] << " " << coords[ind_coord] << "\n";
      continue;
    }
    coordv[0][index] = coords[ind_coord];
    coordv[1][index] = coords[ind_coord + 1];
    coordv[2][index] = coords[ind_coord + 2];
    index_used[index] = 1;
    EntityHandle vertexh = start_vert + index;
    rval = mb->tag_set_data(tagid, &vertexh, 1, (void*) &id);
    ERRORR(rval, "can't set tag id on vertex ");
  }
// create quads; one quad for each edge
  rval = read_iface->get_element_connect(nelem, 4, MBQUAD, 0, start_elem,
      connect);
  ERRORR(rval, "can't create elements ");

  for (int i = 0; i < nelem; i++) {
    for (int j = 0; j < 4; j++) {
      int index_v = index_map[corners[i * 4 + j]];
      connect[i * 4 + j] = start_vert + index_v;
    }
  }

  Range quads(start_elem, start_elem + nelem - 1);

  mb->add_entities(coarseSet, quads);

  // create all edges adjacent to quads
  Range coarseEdges;
  rval = mb->get_adjacencies(quads, 1, true, coarseEdges, Interface::UNION);
  ERRORR(rval, "can't create edges ");

  mb->add_entities(coarseSet, coarseEdges);

  rval = pcomm->resolve_shared_ents(coarseSet, 2, 1); // resolve vertices and edges
  ERRORR(rval, "can't resolve shared vertices and edges ");

  mb->tag_set_data(partitionTag, &coarseSet, 1, &rank);

  /* // do we have any edges?
   Range edges;
   mb->get_entities_by_dimension(0, 1, edges);
   if (rank == 0)
   std::cout << " number of edges after resolve share ents:" << edges.size()
   << "\n";*/

  rval = mb->write_file("coarse.h5m", 0, "PARALLEL=WRITE_PART", &coarseSet, 1);
  ERRORR(rval, "can't write in parallel coarse mesh");
  // coarse mesh, from corners

  return rval;
}
ErrorCode fill_coord_on_edges(Interface * mb, std::vector<double*> & coordv2, double * coords,
    Range & edges, EntityHandle start_v, Range & coarseQuads, int nc, Tag & fineVertOnEdgeTag)
{
  ErrorCode rval=MB_SUCCESS;

  assert(NC==nc);

  int edges_index[4][NC-1]; // indices in the coords array for vertices, oriented positively
  /*
   *  first j, then i, so this is the order of the points in coords array, now:
   *
   *   nc(nc+1),   ...      (nc+1)*(nc+1)-1
   *
   *   2(nc+1),
   *   nc+1,   nc+2,            2*(nc+1)-1
   *   0,      1 ,    2, ..........., nc
   */
  // for first edge, oriented positive, the coords indices (*3) are 1, 2, (nc-1)
  // second edge                                                  2*(nc+1)-1, 3*(nc+1) -1, ..., nc*(nc+1)-1
  // third edge:                                                  (nc+1) * (nc+1) -2, ..., nc*(nc+1) +1
  // fourth edge:                                                 (nc-1)*(nc+1), ..., (nc+1)
  for (int j=1; j<=nc-1; j++)
  {
    edges_index[0][j-1] = j;                        // for nc = 3: 1, 2
    edges_index[1][j-1] = (j+1)*(nc+1) - 1;         //             7, 11
    edges_index[2][j-1] = (nc+1)*(nc+1) - j - 1;    //             14, 13
    edges_index[3][j-1] = (nc-j) * (nc+1) ;         //             8, 4
  }
  //int num_quads=(int)coarseQuads.size();
  int stride = (nc+1)*(nc+1);
  int indexv=0;
  for (Range::iterator eit= edges.begin(); eit!=edges.end(); eit++)
  {
    EntityHandle edge=*eit;
    std::vector<EntityHandle> faces;
    rval = mb->get_adjacencies(&edge, 1, 2, false, faces);
    ERRORR(rval, "can't get adjacent faces.");
    if (faces.size()<1)
      return MB_FAILURE;
    int sense=0, side_number=-1, offset=-1;
    EntityHandle quad=faces[0]; // just consider first quad
    rval = mb->side_number(quad, edge, side_number,  sense, offset);
    ERRORR(rval, "can't get side number");
    int indexq=coarseQuads.index(quad);

    if (indexq==-1)
      return MB_FAILURE;

    EntityHandle firstRefinedV=start_v+indexv;
    rval = mb->tag_set_data(fineVertOnEdgeTag, &edge, 1, &firstRefinedV);
    ERRORR(rval, "can't set refined vertex tag");
    // copy the right coordinates from the coords array to coordv2 array

    double * start_quad = &coords[3*stride*indexq];
    if (sense>0)
    {
      for (int k=1; k<=nc-1; k++)
      {
        int index_in_quad=edges_index[side_number][k-1]*3;
        coordv2[0][indexv]= start_quad[ index_in_quad ];
        coordv2[1][indexv]= start_quad[ index_in_quad + 1];
        coordv2[2][indexv]= start_quad[ index_in_quad + 2];
        indexv++;
      }
    }
    else
    {
      // sense < 0, so we will traverse the edge in inverse sense
      for (int k=1; k<=nc-1; k++)
      {
        int index_in_quad= edges_index[side_number][nc-1-k]*3;
        coordv2[0][indexv]= start_quad[ index_in_quad ];
        coordv2[1][indexv]= start_quad[ index_in_quad + 1];
        coordv2[2][indexv]= start_quad[ index_in_quad + 2];
        indexv++;
      }
    }
  }
  return rval ;
}
ErrorCode create_fine_mesh(Interface * mb, ParallelComm * pcomm,
    EntityHandle coarseSet, EntityHandle fine_set, double * coords,
   int nc, int nelem) {

  int rank = pcomm->proc_config().proc_rank();
  int stride = (nc + 1) * (nc + 1); // 16, for nc ==3
  // there are stride*3 coordinates for each coarse quad, representing the fine mesh

  Tag fineVertTag;
  // will store the first refined vertex on an edge in a entity handle tag
  EntityHandle def = 0;
  ErrorCode rval = mb->tag_get_handle("__firstvertex", 1, MB_TYPE_HANDLE,
      fineVertTag, MB_TAG_DENSE | MB_TAG_CREAT, &def);
  ERRORR(rval, "can't get tag on first vertex");

  Range coarseQuads;
  Range edges;
  rval = mb->get_entities_by_dimension(coarseSet, 2, coarseQuads);
  ERRORR(rval, "can't get coarse quads  ");

  rval = mb->get_entities_by_dimension(coarseSet, 1, edges);
  ERRORR(rval, "can't get coarse edges  ");

  Range verts;
  rval = mb->get_connectivity(coarseQuads, verts);
  ERRORR(rval, "can't get coarse vertices ");

  std::cout <<" local coarse mesh on rank " << rank << "  "<< coarseQuads.size() << " quads, "
     << edges.size() << " edges, " << verts.size() <<  " vertices.\n";

  int dum_id = -1;
  Tag partitionTag;
  mb->tag_get_handle(PARALLEL_PARTITION_TAG_NAME, 1, MB_TYPE_INTEGER,
      partitionTag, MB_TAG_SPARSE | MB_TAG_CREAT, &dum_id);
  // fine mesh, with all coordinates
  std::vector<double *> coordv2;

  ReadUtilIface *read_iface;
  rval = mb->query_interface(read_iface);
  ERRORR(rval, "can't get query intf ");

  EntityHandle start_vert;
  // create verts, (nc+1)*(nc+1)*nelem - verts.size()
  //
  int numVertsOnEdges = (nc-1)*(int)edges.size();
  int numExtraVerts = numVertsOnEdges + (nc-1)*(nc-1)*nelem; // internal fine vertices
  rval = read_iface->get_node_coords(3, numExtraVerts, 0,
      start_vert, coordv2);
  ERRORR(rval, "can't get coords fine mesh");

  // fill coordinates for vertices on the edges, then vertices in the interior of coarse quads
  // we know that all quads are in order, their index corresponds to index in coords array
  rval = fill_coord_on_edges(mb, coordv2, coords, edges, start_vert, coarseQuads, nc, fineVertTag);
  ERRORR(rval, "can't fill edges vertex coords on fine mesh");

  EntityHandle start_elem, *connect;
  rval = read_iface->get_element_connect(nelem * nc * nc, 4, MBQUAD, 0,
      start_elem, connect);
  ERRORR(rval, "can't create elements fine mesh");

  int iv = (nc-1)*edges.size(); // iv is in the coordv2 array indices
  start_vert=start_vert+(nc-1)*edges.size();
  // now fill coordinates on interior nodes; also mark the start for each interior vertex in a coarse quad
  for (int ie = 0; ie < nelem; ie++) {
    // just fill coordinates for an array of (nc-1)*(nc-1) vertices
    EntityHandle firstVert = start_vert+(nc-1)*(nc-1)*ie;
    EntityHandle eh = coarseQuads[ie];
    rval = mb->tag_set_data(fineVertTag, &eh, 1, &firstVert);
    ERRORR(rval, "can't set refined vertex tag");

    int index_coords = stride * ie;
    for (int j = 1; j <= (nc - 1); j++) {
      for (int i = 1; i <= (nc - 1) ; i++)
      {
        int indx2 = 3 * (index_coords + (nc+1) * j + i);
        coordv2[0][iv] = coords[indx2];
        coordv2[1][iv] = coords[indx2 + 1];
        coordv2[2][iv] = coords[indx2 + 2];
        iv++;
      }
    }
  }

  Range quads3(start_elem, start_elem + nelem * nc * nc - 1);

  int ic = 0;
  for (int ie = 0; ie < nelem; ie++) {
    // just fill coordinates for an array of (nc-1)*(nc-1) vertices
    EntityHandle arr2[NC+1][NC+1]; //
    /*
     *     (nc,0)         (nc,nc)
     *
     *
     *     (1,0)           (1,nc)
     *     (0,0) (0,1)     (0,nc)
     */
    EntityHandle coarseQ = coarseQuads[ie];
    const EntityHandle * conn4=NULL;
    int nnodes=0;
    rval = mb->get_connectivity(coarseQ, conn4, nnodes);
    ERRORR(rval, "can't get conn of coarse quad");
    if (nnodes!=4)
      return MB_FAILURE;

    arr2[ 0][ 0] = conn4[0];
    arr2[nc][ 0] = conn4[1];
    arr2[nc][nc] = conn4[2];
    arr2[ 0][nc] = conn4[3];

    // get the coarse edges
    std::vector<EntityHandle> aedges;
    rval = mb->get_adjacencies(&coarseQ, 1, 1, false, aedges);
    ERRORR(rval, "can't get adje edges of coarse quad");
    assert((int)aedges.size()==4);

    for (int k=0; k<4; k++)
    {
      EntityHandle edh = aedges[k];
      /*
          edges_index[0][j-1] = j;                        // for nc = 3: 1, 2
          edges_index[1][j-1] = (j+1)*(nc+1) - 1;         //             7, 11
          edges_index[2][j-1] = (nc+1)*(nc+1) - j - 1;    //             14, 13
          edges_index[3][j-1] = (nc-j) * (nc+1) ;         //             8, 4
       */
      int sense=0, side_number=-1, offset=-1;
      rval = mb->side_number(coarseQ, edh, side_number,  sense, offset);
      ERRORR(rval, "can't get side number");
      EntityHandle firstV; // first vertex on edge, if edge oriented positively
      rval = mb->tag_get_data(fineVertTag, &edh, 1, &firstV);
      ERRORR(rval, "can't get first vertex tag on edge");
      if (sense>0)
      {
        if (0==side_number)
        {
          for (int i=1; i<=nc-1; i++)
            arr2[i][0] = firstV+i-1;
        }
        else if (1==side_number)
        {
          for (int j=1; j<=nc-1; j++)
            arr2[nc][ j] = firstV+j-1;
        }
        else if (2==side_number)
        {
          for (int i=nc-1; i>=1; i--)
            arr2[i][nc] = firstV + nc-1-i;
        }
        else if (3==side_number)
        {
          for (int j=nc-1; j >= 1; j--)
            arr2[0][j] = firstV+nc-1-j;
        }
      }
      else // if (sense<0)
      {
        if (0==side_number)
        {
          for (int i=1; i<=nc-1; i++)
            arr2[i][0] = firstV+nc-1-i;
        }
        else if (1==side_number)
        {
          for (int j=1; j<=nc-1; j++)
            arr2[nc][j] = firstV+nc-1-j;
        }
        else if (2==side_number)
        {
          for (int i=nc-1; i>=1; i--)
            arr2[i][nc] = firstV + i-1;
        }
        else if (3==side_number)
        {
          for (int j=nc-1; j >= 1; j--)
            arr2[0][j] = firstV+ j - 1;
        }
      }
    }
    // fill the interior points matrix
    EntityHandle firstV; // first vertex on interior of coarse quad
    rval = mb->tag_get_data(fineVertTag, &coarseQ, 1, &firstV);
    ERRORR(rval, "can't get first vertex tag on coarse tag");
    int inc=0;
    for (int j=1; j<=nc-1; j++)
    {
      for (int i=1; i<=nc-1; i++)
      {
        arr2[i][j] = firstV+inc;
        inc++;
      }
    }
    // fill now the matrix of quads; vertices are from 0 to (nc+1)*(nc+1)-1
    for (int j1 = 0; j1 < nc; j1++) {
      for (int i1 = 0; i1 < nc; i1++) {
        connect[ic++] = arr2[i1  ][j1  ]; // first one
        connect[ic++] = arr2[i1+1][j1  ]; //
        connect[ic++] = arr2[i1+1][j1+1]; // opp diagonal
        connect[ic++] = arr2[i1  ][j1+1]; //
      }
    }
  }

  mb->add_entities(fine_set, quads3);

  mb->tag_set_data(partitionTag, &fine_set, 1, &rank);

  mb->write_file("fine.h5m", 0, "PARALLEL=WRITE_PART", &fine_set, 1);
  ERRORR(rval, "can't write set 3, fine ");

  return rval;
}

// this is called from Fortran 90
void create_mesh(iMesh_Instance instance,
    iBase_EntitySetHandle * imesh_euler_set, double * coords, int * corners,
    int nc, int nelem, MPI_Fint comm, int * ierr) {
  /* double * coords=(double*) icoords;
   int * corners = (int*) icorners;*/
  *ierr = 1;
  Interface * mb = MOABI;
  MPI_Comm mpicomm = MPI_Comm_f2c(comm);
  // instantiate parallel comm now or not?
  ParallelComm *pcomm = new ParallelComm(mb, mpicomm);

  EntityHandle coarseSet;
  ErrorCode rval = mb->create_meshset(MESHSET_SET, coarseSet);
  ERRORV(rval, "can't create coarse set ");

  rval = create_coarse_mesh(mb, pcomm, coarseSet, coords, corners, nc, nelem);
  ERRORV(rval, "can't create coarse set ");

  EntityHandle fine_set;
  rval = mb->create_meshset(MESHSET_SET, fine_set);
  ERRORV(rval, "can't create coarse set ");

  rval = create_fine_mesh(mb, pcomm, coarseSet, fine_set, coords, nc, nelem);
  ERRORV(rval, "can't create coarse set ");

  *imesh_euler_set = (iBase_EntitySetHandle) fine_set; // with coarse mesh

  *ierr = 0;
  return;
}
#ifdef __cplusplus
} // extern "C"
#endif

