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
struct vertex_id {
  int id;
  int position_corner;
};


void create_mesh(iMesh_Instance instance,
    iBase_EntitySetHandle * imesh_euler_set, double * coords, int * corners,
    int nc, int nelem, MPI_Fint comm, int * ierr) {
  /* double * coords=(double*) icoords;
   int * corners = (int*) icorners;*/
  *ierr = 1;
  Interface * mb = MOABI;
  MPI_Comm  mpicomm = MPI_Comm_f2c(comm);
  // instantiate parallel comm now or not?
  ParallelComm *pcomm= new ParallelComm(mb, mpicomm);

  int rank = pcomm->proc_config().proc_rank();
  EntityHandle euler_set;
  ErrorCode rval = mb->create_meshset(MESHSET_SET, euler_set);
  ERRORV(rval, "can't create euler set ");

  *imesh_euler_set = (iBase_EntitySetHandle) euler_set; // with coarse mesh

  Tag tagid; // global id at corners
  rval = mb->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, tagid,
      MB_TAG_DENSE | MB_TAG_CREAT);
  ERRORV(rval, "can't create corner tag ");

  int dum_id = -1;
  Tag partitionTag;
  mb->tag_get_handle(PARALLEL_PARTITION_TAG_NAME,
                                                  1, MB_TYPE_INTEGER, partitionTag,
                                                    MB_TAG_SPARSE|MB_TAG_CREAT, &dum_id);

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
  ERRORV(rval, "can't get query intf ");

  std::vector<double *> coordv;
  EntityHandle start_vert, start_elem, *connect;
// create verts, num is 2(nquads+1) because they're in a 1d row; will initialize coords in loop over quads later
  rval = read_iface->get_node_coords(3, num_nodes_coarse, 0, start_vert,
      coordv);
  ERRORV(rval, "can't get node coords ");
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
  int order_in_coord[] = { 0, nc , (nc + 1) * (nc + 1) - 1, nc * (nc + 1)};
  for (int i = 0; i < size_corners; i++) {
    int id = corners[i];
    int index = index_map[id];

    int j = i % 4;
    int ind_coord = i / 4;
    ind_coord = (ind_coord * stride + order_in_coord[j]) * 3;
    if (index_used[index])
    {
      if ( fabs(coordv[0][index]-coords[ind_coord])>0.000001 )
        std::cout<<" id:" << corners[i] << " i:" << i << " j:" << j << " " <<
        coordv[0][index] << " " << coords[ind_coord] << "\n";
      continue;
    }
    coordv[0][index] = coords[ind_coord];
    coordv[1][index] = coords[ind_coord + 1];
    coordv[2][index] = coords[ind_coord + 2];
    index_used[index]=1;
    EntityHandle vertexh = start_vert + index;
    rval = mb->tag_set_data(tagid, &vertexh, 1, (void*) &id);
    ERRORV(rval, "can't set tag id on vertex ");
  }
// create quads; one quad for each edge
  rval = read_iface->get_element_connect(nelem, 4, MBQUAD, 0, start_elem,
      connect);
  ERRORV(rval, "can't create elements ");

  for (int i = 0; i < nelem; i++) {
    for (int j = 0; j < 4; j++) {
      int index_v = index_map[corners[i * 4 + j]];
      connect[i * 4 + j] = start_vert + index_v;
    }
  }

  Range quads(start_elem, start_elem + nelem - 1);

  mb->add_entities(euler_set, quads);

  rval = pcomm->resolve_shared_ents(euler_set, 2, 0); // resolve just vertices
  ERRORV(rval, "can't resolve shared vertices ");

  mb->tag_set_data(partitionTag, &euler_set, 1, &rank);

  rval = mb->write_file("coarse.h5m", 0, "PARALLEL=WRITE_PART", &euler_set, 1);
  ERRORV(rval, "can't write in parallel coarse mesh");
  // coarse mesh, from corners

  // fine mesh, with all coordinates
  std::vector<double *> coordv2;

  // create verts, (nc+1)*(nc+1)*nelem
  rval = read_iface->get_node_coords(3, nelem * (nc + 1) * (nc + 1), 0,
      start_vert, coordv2);
  ERRORV(rval, "can't get coords fine mesh");

  rval = read_iface->get_element_connect(nelem * nc * nc, 4, MBQUAD, 0,
      start_elem, connect);
  ERRORV(rval, "can't create elements fine mesh");
  Range quads3(start_elem, start_elem + nelem * nc * nc - 1);
  int iv = 0;
  int ic = 0;
  for (int i = 0; i < nelem; i++) {
    // just duplicate nodes,
    int index_coords = stride * i;
    for (int j = 0; j < (nc + 1) * (nc + 1); j++) {
      int indx2 = 3 * (index_coords + j);
      coordv2[0][iv] = coords[indx2];
      coordv2[1][iv] = coords[indx2 + 1];
      coordv2[2][iv] = coords[indx2 + 2];
      iv++;
    }
    EntityHandle v00 = start_vert + stride * i;
    // fill now the matrix of quads; vertices are from 0 to (nc+1)*(nc+1)-1
    for (int j1 = 0; j1 < nc; j1++) {
      for (int i1 = 0; i1 < nc; i1++) {
        connect[ic++] = v00 + j1 * (nc + 1)       + i1; // first one
        connect[ic++] = v00 + j1 * (nc + 1)       + i1 + 1; //
        connect[ic++] = v00 + (j1 + 1) * (nc + 1) + i1 + 1; // opp diagonal
        connect[ic++] = v00 + (j1 + 1) * (nc + 1) + i1 ; //
      }
    }
  }
  EntityHandle set3;
  rval = mb->create_meshset(MESHSET_SET, set3);

  ERRORV(rval, "can't create set 3, fine mesh ");

  mb->add_entities(set3, quads3);



  mb->tag_set_data(partitionTag, &set3, 1, &rank);

  mb->write_file("fine.h5m", 0, "PARALLEL=WRITE_PART", &set3, 1);
  ERRORV(rval, "can't write set 3, fine ");

  *ierr = 0;
  return;
}
#ifdef __cplusplus
} // extern "C"
#endif

