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
#include "moab/ParallelMergeMesh.hpp"
#include <sstream>
#include <mpi.h>

using namespace moab;
double radius = 1.;
double gtol = 1.e-9;
bool debug = false;

// this mapping to coordinates will keep an index into the coords and dep_coords array
// more exactly, the fine vertices are in a Range fineVerts;
// then the vertex index i in fineVerts , fineVerts[i] will have 3 coordinates at
// index  mapping_to_coords[i] * 3
// more exactly, x = coords[ mapping_to_coords[i] * 3    ]
//               y = coords[ mapping_to_coords[i] * 3 + 1]
//               z = coords[ mapping_to_coords[i] * 3 + 2]

// in the same way, departure position for vertex depVerts[i] will be at index
//   mapping_to_coords[i] * 2 in dep_coords array (it has just latitude and longitude)
//

int * mapping_to_coords = NULL;
int numVertices = 0;
// this will be instantiated at create mesh step
// should be cleaned up at the end
Intx2MeshOnSphere * pworker = NULL;

Range allCells; /// these will be cells on current task, after ghost exchange, and they will be used for tags exchange

// more variables that should be saved between time steps, etc
Tag areaTag=0;
Tag rhoTag = 0;
Tag barycenterTag = 0;
Tag centerOfMassTag = 0;
Tag cellIntTag = 0;
Tag rhoCoefTag = 0;
Tag planeTag = 0;
Tag rhoBoundsTag = 0;
Tag weightsTag = 0;
Tag gid = 0;
ParallelComm *pcomm = NULL;
// for printing, debug style
int num_update_calls=0;
int nlev = 4;
int level_to_print =1; // + 1
bool dprint = false;
// end printing
// should get rid of this; instead of using array[NC+1][NC+1], use  row based indexing (C-style):
// parray =  new int[ (NC+1)*(NC+1) ] , and instead of array[i][j], use parray[ i*(NC+1) + j ]
#define  NC  3

/*
 *  methods defined here:
 *  void update_tracer(iMesh_Instance instance,
 iBase_EntitySetHandle imesh_euler_set, int * ierr);

 void create_mesh(iMesh_Instance instance,
 iBase_EntitySetHandle * imesh_euler_set, double * coords, int * corners,
 int nc, int nelem, MPI_Fint comm, int * ierr) ;

 void intersection_at_level(iMesh_Instance instance,
 iBase_EntitySetHandle fine_set, iBase_EntitySetHandle * intx_set, double * dep_coords, double radius,
 int nc, int nelem, MPI_Fint comm, int * ierr)

 */
#ifdef __cplusplus
extern "C" {
#endif

// this will initialize more, integrals / weights cells, density tag, etc
void initialize_area_and_tracer(iMesh_Instance instance,
    iBase_EntitySetHandle imesh_euler_set, double * area_vals, int * ierr) {

  EntityHandle eul_set = (EntityHandle) imesh_euler_set;

  moab::Interface * mb = MOABI;
  *ierr = 1;

  ErrorCode rval;

  Range eulQuads;
  rval = mb->get_entities_by_type(eul_set, MBQUAD, eulQuads);
  ERRORV(rval, "can't get eulerian quads");

  // area of the euler element is fixed, store it; it is used to recompute the averages at each
  // time step
  std::string tag_name4("Area");
  rval = mb->tag_get_handle(tag_name4.c_str(), 1, MB_TYPE_DOUBLE, areaTag,
      MB_TAG_DENSE | MB_TAG_CREAT);
  ERRORV(rval, "can't get area tag");
#ifndef NDEBUG
  std::cout << " num quads = " << eulQuads.size() << "\n";
#endif
  ERRORV(rval, "can't set area for each quad");

  rval = mb->tag_set_data(areaTag, eulQuads, &area_vals[0]); MB_CHK_ERR_RET(rval);
  // other tags density, bounds, linear appx for it

  rval = mb->tag_get_handle("Density", 1, MB_TYPE_DOUBLE, rhoTag,
        MB_TAG_CREAT | MB_TAG_DENSE);MB_CHK_ERR_RET(rval);

  // Create tag for cell barycenters in 3D Cartesian space
  rval = mb->tag_get_handle("CellBarycenter", 3, MB_TYPE_DOUBLE,
        barycenterTag, MB_TAG_CREAT | MB_TAG_DENSE); MB_CHK_ERR_RET(rval);

  // Create tag for cell centers of mass in 3D Cartesian space

  rval = mb->tag_get_handle("CellCenterOfMass", 3, MB_TYPE_DOUBLE,
        centerOfMassTag, MB_TAG_CREAT | MB_TAG_DENSE); MB_CHK_ERR_RET(rval);

  // Create tag for integrals of (x, y, x^2, y^2, xy) over Eulerian cells
  rval = mb->tag_get_handle("CellIntegral", 6, MB_TYPE_DOUBLE, cellIntTag,
        MB_TAG_CREAT | MB_TAG_DENSE); MB_CHK_ERR_RET(rval);

    // Create tag for cell density reconstruction coefficients

  rval = mb->tag_get_handle("LinearCoefRho", 3, MB_TYPE_DOUBLE, rhoCoefTag,
        MB_TAG_CREAT | MB_TAG_DENSE); MB_CHK_ERR_RET(rval);


  // Create tag for index of gnomonic plane for each cell
  rval = mb->tag_get_handle("GnomonicPlane", 1, MB_TYPE_INTEGER, planeTag,
        MB_TAG_CREAT | MB_TAG_DENSE); MB_CHK_ERR_RET(rval);

  // Create tag for density bounds
  rval = mb->tag_get_handle("DensityBounds", 2, MB_TYPE_DOUBLE,
        rhoBoundsTag, MB_TAG_CREAT | MB_TAG_DENSE); MB_CHK_ERR_RET(rval);

  // Create tag for intersection weights
  rval = mb->tag_get_handle("Weights", 6, moab::MB_TYPE_DOUBLE, weightsTag,
        MB_TAG_CREAT | MB_TAG_DENSE); MB_CHK_ERR_RET(rval);

  // get global id tag
  rval = mb->tag_get_handle(GLOBAL_ID_TAG_NAME, 1,
       MB_TYPE_INTEGER, gid, MB_TAG_DENSE); MB_CHK_ERR_RET(rval);

   // get cell plane
   get_gnomonic_plane(mb, eul_set, planeTag);

   // get Eulerian cell integrals over x, y, x^2, y^2, xy (and cell area)
   get_eul_cell_integrals(mb, eul_set, planeTag, areaTag, cellIntTag);

   // get cell barycenters
   get_barycenters(mb, eul_set, planeTag, cellIntTag, barycenterTag);

   if (pcomm)
   {
     rval = pcomm->exchange_tags(barycenterTag, allCells);MB_CHK_ERR_RET(rval);
   }
   // allCells have
  *ierr = 0;
  return;
}

void update_tracer_test(iMesh_Instance instance,
    iBase_EntitySetHandle imesh_euler_set,
    iBase_EntitySetHandle imesh_output_set, int numTracers,
    int transport_type, double * density_vals, double * tracer_vals, int * ierr) {

  moab::Interface * mb = MOABI;
  *ierr = 1;

  // some tags will depend on the number of tracers
  // Create tag for cell tracer
  Tag tauTag = 0;
  Tag tauCoefTag = 0;
  Tag tauBoundsTag = 0;

  ErrorCode rval = mb->tag_get_handle("Tracer", numTracers, MB_TYPE_DOUBLE, tauTag,
      MB_TAG_CREAT | MB_TAG_DENSE); MB_CHK_ERR_RET(rval);

  if (7==transport_type) // linear reconstruction
  {
    // really, this should be created only once; also linear coeffs, boundaries, etc
    // Create/retrieve tag for cell tracer reconstruction coefficients
    rval = mb->tag_get_handle("LinearCoefTau", 3*numTracers, moab::MB_TYPE_DOUBLE, tauCoefTag,
        MB_TAG_CREAT | MB_TAG_DENSE); MB_CHK_ERR_RET(rval);
    // Create/retrieve tag for tracer bounds
    rval = mb->tag_get_handle("TracerBounds", 2*numTracers, moab::MB_TYPE_DOUBLE,
          tauBoundsTag, moab::MB_TAG_CREAT | moab::MB_TAG_DENSE); MB_CHK_ERR_RET(rval);
  }

  EntityHandle eul_set = (EntityHandle) imesh_euler_set;
  EntityHandle output_set = (EntityHandle) imesh_output_set;

  Range eulQuads;
  rval = mb->get_entities_by_type(eul_set, MBQUAD, eulQuads);
  ERRORV(rval, "can't get eulerian quads");

  // area of the euler element is fixed, store it; it is used to recompute the averages at each
  // time step
#ifndef NDEBUG
  std::cout << " num quads = " << eulQuads.size() << "\n";
#endif
  rval = mb->tag_set_data(tauTag, eulQuads, &tracer_vals[0]);
  ERRORV(rval, "can't set tracer data");

  if (6==transport_type)
  {
    rval = pworker->update_tracer_data(output_set, tauTag, areaTag);
    ERRORV(rval, "can't update tracer ");
  }
  else
  {
    // if 7 == transport_type
    rval = mb->tag_set_data(rhoTag, eulQuads, &density_vals[0]); MB_CHK_ERR_RET(rval);

    std::vector<Tag> tags;
    tags.push_back(rhoTag); tags.push_back(tauTag);
    if (dprint)
    {
      num_update_calls++;
      if (num_update_calls%nlev==level_to_print)
      {
        std::vector<Tag> wtags;
        wtags.push_back(rhoTag);
        for (int k=1; k<=numTracers; k++)
        {
          std::stringstream tagName;
          tagName << "Tracer" << k ;
          Tag tracTag;
          rval = mb->tag_get_handle(tagName.str().c_str(), 1, moab::MB_TYPE_DOUBLE,
                    tracTag, moab::MB_TAG_CREAT | moab::MB_TAG_DENSE); MB_CHK_ERR_RET(rval);
          wtags.push_back(tracTag);
          std::vector<double> pickVals(eulQuads.size());
          for (size_t i=0; i<eulQuads.size(); i++)
          {
            pickVals[i] = tracer_vals[(k-1)+i*numTracers];
          }
          rval = mb->tag_set_data(tracTag, eulQuads, &pickVals[0]);MB_CHK_ERR_RET(rval);
        }

        int nstep=num_update_calls/nlev;
        std::stringstream fff;
        fff << "before0" << nstep << ".h5m";

        rval = mb->write_file(fff.str().c_str(), 0, "PARALLEL=WRITE_PART", &eul_set, 1, &wtags[0], wtags.size()); MB_CHK_ERR_RET(rval);
      }
    }
    rval = pcomm->exchange_tags(tags, tags, allCells);

    // start copy
    // get density bounds, filter the adjacent cells with allCells ()
    get_neighborhood_bounds(mb, eul_set, rhoTag, rhoBoundsTag, &allCells);

    // get tracer bounds
    get_neighborhood_bounds(mb, eul_set, tauTag, tauBoundsTag, &allCells);

       // get linear reconstruction coefficients for density
    get_linear_reconstruction(mb, eul_set, rhoTag, planeTag, barycenterTag,
           rhoCoefTag, &allCells); // also filter for adjacent cells

    // limit linear reconstruction coefficients for density
    limit_linear_reconstruction(mb, eul_set, rhoTag, planeTag, rhoBoundsTag,
           barycenterTag, rhoCoefTag);

       // get cell centers of mass
    get_centers_of_mass(mb, eul_set, planeTag, rhoCoefTag, cellIntTag,
           centerOfMassTag);

     // need to communicate center of mass tag for ghost elements, too
     // it is needed for linear reconstruction of tracer coeffs
    rval = pcomm->exchange_tags(centerOfMassTag, allCells); MB_CHK_ERR_RET(rval);
     // get linear reconstruction coefficients for tracer
    get_linear_reconstruction(mb, eul_set, tauTag, planeTag, centerOfMassTag,
         tauCoefTag, &allCells);

       // limit linear reconstruction coefficients for tracer
    limit_linear_reconstruction(mb, eul_set, tauTag, planeTag, tauBoundsTag,
           centerOfMassTag, tauCoefTag);

    // intersection weights (i.e. area, x integral, and y integral over cell intersections)
    get_intersection_weights(mb, eul_set, output_set, planeTag,
        weightsTag);

    // update the density and tracer
    rval = pworker->update_density_and_tracers(rhoTag, areaTag, rhoCoefTag,
            tauTag, tauCoefTag, weightsTag, planeTag); MB_CHK_ERR_RET(rval);

    rval = mb->tag_get_data(rhoTag, eulQuads, &density_vals[0]); MB_CHK_ERR_RET(rval);
    //end copy
  }

  rval = mb->tag_get_data(tauTag, eulQuads, &tracer_vals[0]);
  ERRORV(rval, "can't get tracer data");

  if (dprint && 7==transport_type)
  {

    if (num_update_calls%nlev==level_to_print)
    {
      std::vector<Tag> wtags;
      wtags.push_back(rhoTag);
      for (int k=1; k<=numTracers; k++)
      {
        std::stringstream tagName;
        tagName << "Tracer" << k ;
        Tag tracTag;
        rval = mb->tag_get_handle(tagName.str().c_str(), 1, moab::MB_TYPE_DOUBLE,
                  tracTag, moab::MB_TAG_CREAT | moab::MB_TAG_DENSE); MB_CHK_ERR_RET(rval);
        wtags.push_back(tracTag);
        std::vector<double> pickVals(eulQuads.size());
        for (size_t i=0; i<eulQuads.size(); i++)
        {
          pickVals[i] = tracer_vals[(k-1)+i*numTracers];
        }
        rval = mb->tag_set_data(tracTag, eulQuads, &pickVals[0]);MB_CHK_ERR_RET(rval);
      }
      std::vector<double> tcvals(3*numTracers*eulQuads.size());
      rval = mb->tag_get_data(tauBoundsTag, eulQuads, &tcvals[0]);MB_CHK_ERR_RET(rval);
      for (int k=1; k<=numTracers; k++)
      {
        std::stringstream tagName;
        tagName << "Bounds" << k ;
        Tag bTag;
        rval = mb->tag_get_handle(tagName.str().c_str(), 2, moab::MB_TYPE_DOUBLE,
                  bTag, moab::MB_TAG_CREAT | moab::MB_TAG_DENSE); MB_CHK_ERR_RET(rval);
        wtags.push_back(bTag);
        std::vector<double> pickVals(2*eulQuads.size());
        for (size_t i=0; i<eulQuads.size(); i++)
        {
          pickVals[2*i] = tcvals[2*(k-1)+i*numTracers*2];
          pickVals[2*i+1] = tcvals[2*(k-1)+1+i*numTracers*2];
        }
        rval = mb->tag_set_data(bTag, eulQuads, &pickVals[0]);MB_CHK_ERR_RET(rval);
      }
      rval = mb->tag_get_data(tauCoefTag, eulQuads, &tcvals[0]);MB_CHK_ERR_RET(rval);
      for (int k=1; k<=numTracers; k++)
      {
        std::stringstream tagName;
        tagName << "TagCoef" << k ;
        Tag coTag;
        rval = mb->tag_get_handle(tagName.str().c_str(), 3, moab::MB_TYPE_DOUBLE,
                  coTag, moab::MB_TAG_CREAT | moab::MB_TAG_DENSE); MB_CHK_ERR_RET(rval);
        wtags.push_back(coTag);
        std::vector<double> pickVals(3*eulQuads.size());
        for (size_t i=0; i<eulQuads.size(); i++)
        {
          pickVals[3*i] = tcvals[3*(k-1)+i*numTracers*3];
          pickVals[3*i+1] = tcvals[3*(k-1)+1+i*numTracers*3];
          pickVals[3*i+2] = tcvals[3*(k-1)+2+i*numTracers*3];
        }
        rval = mb->tag_set_data(coTag, eulQuads, &pickVals[0]);MB_CHK_ERR_RET(rval);
      }

      int nstep=num_update_calls/nlev;
      std::stringstream fff;
      fff << "after0" << nstep << ".h5m";

      rval = mb->write_file(fff.str().c_str(), 0, "PARALLEL=WRITE_PART", &eul_set, 1, &wtags[0], wtags.size()); MB_CHK_ERR_RET(rval);
    }
  }

  *ierr = 0;
  return;
}

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
  Tag tagArea2 = 0;
  std::string tag_name4("Area");
  rval = mb->tag_get_handle(tag_name4.c_str(), 1, MB_TYPE_DOUBLE, tagArea2,
      MB_TAG_DENSE | MB_TAG_CREAT);
  ERRORV(rval, "can't get area tag");

  rval = worker.update_tracer_data(outputSet, tagElem, tagArea2);
  ERRORV(rval, "can't update tracer ");

  // everything can be deleted now from intx data; polygons, etc.

  *ierr = 0;
  return;
}

ErrorCode create_coarse_mesh(Interface * mb, ParallelComm * pcomm1,
    EntityHandle coarseSet, double * coords, int * corners, int nc, int nelem,
    EntityHandle & start_vert, int & totalNumVertices, int & numCornerVertices,
    std::vector<double *> & coordv) {

  int rank = pcomm1->proc_config().proc_rank();

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

  numCornerVertices = (int) corn1.size();

  // these corners are actually dofs in spectral mesh, so they will vary from 0 to np*np*nel_global (at most)
  // they are shared on edges, corners of spectral elements, so number of dofs is lower than this estimate
  // still, the number of coarse vertices is close to number of cells + 2

  std::map<int, int> index_map;
  for (size_t i = 0; i < corn1.size(); i++)
    index_map[corn1[i]] = i;

  // estimate for the vertices that will get created, from euler formula
  // v-e+f = 2 for whole sphere
  // for one task (distributed) v-e+f = 1, for one connected region
  // if multiple connex reqions (out of HFSC distribution from homme), v-e+f = k, k is number of connectivity regs
  // so in any case, e = v+f -k, where k is at least 1 (1, 2, 3, etc)
  // so in any case number of coarse edges is at most e_max = v+f-1

  // this will give an estimate for the number of "fine" vertices
  int e_max = nelem + numCornerVertices - 1;
  // total number of extra vertices will be
  int numVertsOnEdges = (nc - 1) * e_max;
  int numExtraVerts = numVertsOnEdges + (nc - 1) * (nc - 1) * nelem; // internal fine vertices

  totalNumVertices = numCornerVertices + numExtraVerts; // this could be overestimated, because we are not sure
  // about the number of edges

  // used to determine the if the nodes are matching at corners of elements
  std::vector<int> index_used(numCornerVertices, 0);

  ReadUtilIface *read_iface;
  rval = mb->query_interface(read_iface);
  ERRORR(rval, "can't get query intf ");

  EntityHandle start_elem, *connect;
// create verts, num is 2(nquads+1) because they're in a 1d row; will initialize coords in loop over quads later
  rval = read_iface->get_node_coords(3, totalNumVertices, 0, start_vert,
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
// create quads; one quad for each spectral element; later, we will create fine quads
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

  // see how much we overestimated the number e_max
#ifndef NDEBUG
  std::cout << " on rank " << rank << " e_max is " << e_max
      << " actual number of edges: " << coarseEdges.size() << "\n";
#endif
  rval = pcomm1->resolve_shared_ents(coarseSet, 2, 1); // resolve vertices and edges
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

// start_v and coordv refer to all vertices, including the coarse ones
ErrorCode fill_coord_on_edges(Interface * mb, std::vector<double*> & coordv,
    double * coords, Range & edges, EntityHandle start_v, Range & coarseQuads,
    int nc, int numCornerVertices, Tag & fineVertOnEdgeTag) {
  ErrorCode rval = MB_SUCCESS;

  double * coordv2[3];
  for (int k = 0; k < 3; k++)
    coordv2[k] = coordv[k] + numCornerVertices; // they will start later

  assert(NC==nc);

  int edges_index[4][NC - 1]; // indices in the coords array for vertices, oriented positively
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
  for (int j = 1; j <= nc - 1; j++) {
    edges_index[0][j - 1] = j;                        // for nc = 3: 1, 2
    edges_index[1][j - 1] = (j + 1) * (nc + 1) - 1;         //             7, 11
    edges_index[2][j - 1] = (nc + 1) * (nc + 1) - j - 1;   //             14, 13
    edges_index[3][j - 1] = (nc - j) * (nc + 1);         //             8, 4
  }
  //int num_quads=(int)coarseQuads.size();
  int stride = (nc + 1) * (nc + 1);
  int indexv = 0;
  for (Range::iterator eit = edges.begin(); eit != edges.end(); eit++) {
    EntityHandle edge = *eit;
    std::vector<EntityHandle> faces;
    rval = mb->get_adjacencies(&edge, 1, 2, false, faces);
    ERRORR(rval, "can't get adjacent faces.");
    if (faces.size() < 1)
      return MB_FAILURE;
    int sense = 0, side_number = -1, offset = -1;
    EntityHandle quad = faces[0]; // just consider first quad
    rval = mb->side_number(quad, edge, side_number, sense, offset);
    ERRORR(rval, "can't get side number");
    int indexq = coarseQuads.index(quad);

    if (indexq == -1)
      return MB_FAILURE;

    EntityHandle firstRefinedV = start_v + numCornerVertices + indexv;
    rval = mb->tag_set_data(fineVertOnEdgeTag, &edge, 1, &firstRefinedV);
    ERRORR(rval, "can't set refined vertex tag");
    // copy the right coordinates from the coords array to coordv2 array

    double * start_quad = &coords[3 * stride * indexq];
    if (sense > 0) {
      for (int k = 1; k <= nc - 1; k++) {
        int index_in_quad = edges_index[side_number][k - 1] * 3;
        coordv2[0][indexv] = start_quad[index_in_quad];
        coordv2[1][indexv] = start_quad[index_in_quad + 1];
        coordv2[2][indexv] = start_quad[index_in_quad + 2];
        indexv++;
      }
    } else {
      // sense < 0, so we will traverse the edge in inverse sense
      for (int k = 1; k <= nc - 1; k++) {
        int index_in_quad = edges_index[side_number][nc - 1 - k] * 3;
        coordv2[0][indexv] = start_quad[index_in_quad];
        coordv2[1][indexv] = start_quad[index_in_quad + 1];
        coordv2[2][indexv] = start_quad[index_in_quad + 2];
        indexv++;
      }
    }
  }
  return rval;
}
/*
 ErrorCode resolve_interior_verts_on_bound_edges(Interface * mb, ParallelComm * pcomm,
 Range & edges)
 {
 // edges are coarse edges;
 ErrorCode rval;
 int rank=pcomm->proc_config().proc_rank();
 Range sharedCoarseEdges=edges;// filter the non shared ones
 rval = pcomm->filter_pstatus(sharedCoarseEdges, PSTATUS_SHARED, PSTATUS_AND);
 ERRORR(rval, "can't filter coarse edges  ");
 ParallelMergeMesh pmerge(pcomm, 0.0001);
 ErrorCode rval = pmerge.merge();

 return rval;
 }*/
ErrorCode create_fine_mesh(Interface * mb, ParallelComm * pcomm1,
    EntityHandle coarseSet, EntityHandle fine_set, double * coords, int nc,
    int nelem, int transport_type, EntityHandle start_vert, int numCornerVertices,
    std::vector<double *> & coordv) {
  int rank = pcomm1->proc_config().proc_rank();
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

  /*std::cout <<" local coarse mesh on rank " << rank << "  "<< coarseQuads.size() << " quads, "
   << edges.size() << " edges, " << verts.size() <<  " vertices.\n";*/

  int dum_id = -1;
  Tag partitionTag;
  mb->tag_get_handle(PARALLEL_PARTITION_TAG_NAME, 1, MB_TYPE_INTEGER,
      partitionTag, MB_TAG_SPARSE | MB_TAG_CREAT, &dum_id);
  // fine mesh, with all coordinates
  // std::vector<double *> coordv2;

  ReadUtilIface *read_iface;
  rval = mb->query_interface(read_iface);
  ERRORR(rval, "can't get query intf ");

  ;
  // create verts, (nc+1)*(nc+1)*nelem - verts.size()
  //
  /* int numVertsOnEdges = (nc-1)*(int)edges.size();
   int numExtraVerts = numVertsOnEdges + (nc-1)*(nc-1)*nelem; // internal fine vertices
   rval = read_iface->get_node_coords(3, numExtraVerts, 0,
   start_vert, coordv2);*/
  // ERRORR(rval, "can't get coords fine mesh");
  // fill coordinates for vertices on the edges, then vertices in the interior of coarse quads
  // we know that all quads are in order, their index corresponds to index in coords array
  rval = fill_coord_on_edges(mb, coordv, coords, edges, start_vert, coarseQuads,
      nc, numCornerVertices, fineVertTag);
  ERRORR(rval, "can't fill edges vertex coords on fine mesh");

  EntityHandle start_elem, *connect;
  rval = read_iface->get_element_connect(nelem * nc * nc, 4, MBQUAD, 0,
      start_elem, connect);
  ERRORR(rval, "can't create elements fine mesh");

  int iv = (nc - 1) * edges.size() + numCornerVertices; // iv is in the coordv array indices
  start_vert = start_vert + numCornerVertices + (nc - 1) * edges.size();

  /* // add a child to the mesh set, with the ordered vertices, as they come out in the list of coordinates
   EntityHandle vertSet;
   rval = mb->create_meshset(MESHSET_ORDERED, vertSet);
   ERRORR(rval, "can't create vertex set ");

   rval = mb->add_parent_child(fine_set, vertSet);
   ERRORR(rval, "can't create parent child relation between fine set and vertSet ");*/

  std::vector<EntityHandle> vertList;
  vertList.reserve(nelem * (nc + 1) * (nc + 1)); // will have a list of vertices, in order

  // now fill coordinates on interior nodes; also mark the start for each interior vertex
  //  in a coarse quad
  for (int ie = 0; ie < nelem; ie++) {
    // just fill coordinates for an array of (nc-1)*(nc-1) vertices
    EntityHandle firstVert = start_vert + (nc - 1) * (nc - 1) * ie;
    EntityHandle eh = coarseQuads[ie];
    rval = mb->tag_set_data(fineVertTag, &eh, 1, &firstVert);
    ERRORR(rval, "can't set refined vertex tag");

    int index_coords = stride * ie;
    for (int j = 1; j <= (nc - 1); j++) {
      for (int i = 1; i <= (nc - 1); i++) {
        int indx2 = 3 * (index_coords + (nc + 1) * j + i);
        coordv[0][iv] = coords[indx2];
        coordv[1][iv] = coords[indx2 + 1];
        coordv[2][iv] = coords[indx2 + 2];
        iv++;
      }
    }
  }

  Range quads3(start_elem, start_elem + nelem * nc * nc - 1);

  int ic = 0;
  for (int ie = 0; ie < nelem; ie++) {
    // just fill coordinates for an array of (nc-1)*(nc-1) vertices
    EntityHandle arr2[NC + 1][NC + 1]; //
    /*
     *     (nc,0)         (nc,nc)
     *
     *
     *     (1,0)           (1,nc)
     *     (0,0) (0,1)     (0,nc)
     */
    EntityHandle coarseQ = coarseQuads[ie];
    const EntityHandle * conn4 = NULL;
    int nnodes = 0;
    rval = mb->get_connectivity(coarseQ, conn4, nnodes);
    ERRORR(rval, "can't get conn of coarse quad");
    if (nnodes != 4)
      return MB_FAILURE;

    arr2[0][0] = conn4[0];
    arr2[nc][0] = conn4[1];
    arr2[nc][nc] = conn4[2];
    arr2[0][nc] = conn4[3];

    // get the coarse edges
    std::vector<EntityHandle> aedges;
    rval = mb->get_adjacencies(&coarseQ, 1, 1, false, aedges);
    ERRORR(rval, "can't get adje edges of coarse quad");
    assert((int )aedges.size() == 4);

    for (int k = 0; k < 4; k++) {
      EntityHandle edh = aedges[k];
      /*
       edges_index[0][j-1] = j;                        // for nc = 3: 1, 2
       edges_index[1][j-1] = (j+1)*(nc+1) - 1;         //             7, 11
       edges_index[2][j-1] = (nc+1)*(nc+1) - j - 1;    //             14, 13
       edges_index[3][j-1] = (nc-j) * (nc+1) ;         //             8, 4
       */
      int sense = 0, side_number = -1, offset = -1;
      rval = mb->side_number(coarseQ, edh, side_number, sense, offset);
      ERRORR(rval, "can't get side number");
      EntityHandle firstV; // first vertex on edge, if edge oriented positively
      rval = mb->tag_get_data(fineVertTag, &edh, 1, &firstV);
      ERRORR(rval, "can't get first vertex tag on edge");
      if (sense > 0) {
        if (0 == side_number) {
          for (int i = 1; i <= nc - 1; i++)
            arr2[i][0] = firstV + i - 1;
        } else if (1 == side_number) {
          for (int j = 1; j <= nc - 1; j++)
            arr2[nc][j] = firstV + j - 1;
        } else if (2 == side_number) {
          for (int i = nc - 1; i >= 1; i--)
            arr2[i][nc] = firstV + nc - 1 - i;
        } else if (3 == side_number) {
          for (int j = nc - 1; j >= 1; j--)
            arr2[0][j] = firstV + nc - 1 - j;
        }
      } else // if (sense<0)
      {
        if (0 == side_number) {
          for (int i = 1; i <= nc - 1; i++)
            arr2[i][0] = firstV + nc - 1 - i;
        } else if (1 == side_number) {
          for (int j = 1; j <= nc - 1; j++)
            arr2[nc][j] = firstV + nc - 1 - j;
        } else if (2 == side_number) {
          for (int i = nc - 1; i >= 1; i--)
            arr2[i][nc] = firstV + i - 1;
        } else if (3 == side_number) {
          for (int j = nc - 1; j >= 1; j--)
            arr2[0][j] = firstV + j - 1;
        }
      }
    }
    // fill the interior points matrix
    EntityHandle firstV; // first vertex on interior of coarse quad
    rval = mb->tag_get_data(fineVertTag, &coarseQ, 1, &firstV);
    ERRORR(rval, "can't get first vertex tag on coarse tag");
    int inc = 0;
    for (int j = 1; j <= nc - 1; j++) {
      for (int i = 1; i <= nc - 1; i++) {
        arr2[i][j] = firstV + inc;
        inc++;
      }
    }
    // fill now the matrix of quads; vertices are from 0 to (nc+1)*(nc+1)-1
    for (int j1 = 0; j1 < nc; j1++) {
      for (int i1 = 0; i1 < nc; i1++) {
        connect[ic++] = arr2[i1][j1]; // first one
        connect[ic++] = arr2[i1 + 1][j1]; //
        connect[ic++] = arr2[i1 + 1][j1 + 1]; // opp diagonal
        connect[ic++] = arr2[i1][j1 + 1]; //
      }
    }

    for (int j1 = 0; j1 <= nc; j1++) {
      for (int i1 = 0; i1 <= nc; i1++) {
        vertList.push_back(arr2[i1][j1]);
      }
    }
  }
  /*
   rval = mb->add_entities(vertSet, &vertList[0], (int)vertList.size());
   ERRORR(rval,"can't add to the vert set the list of ordered vertices");*/

  mb->add_entities(fine_set, quads3);
  // notify MOAB of the new elements
  rval = read_iface->update_adjacencies(start_elem, nelem * nc * nc, 4,
      connect);
  ERRORR(rval, "can't update adjacencies on fine quads");

  rval = mb->tag_set_data(partitionTag, &fine_set, 1, &rank);
  ERRORR(rval, "can't set partition tag on fine set");

  // the vertices on the boundary edges of the partition need to be shared and resolved
  ParallelMergeMesh pmerge(pcomm1, 0.0001);
  rval = pmerge.merge();
  ERRORR(rval, "can't resolve vertices on interior of boundary edges ");

  // delete the coarse mesh, except vertices
  pcomm1->delete_entities(coarseQuads);
  pcomm1->delete_entities(edges);

  rval = mb->get_connectivity(quads3, verts);
  ERRORR(rval, "can't get vertices ");

  Range owned_verts = verts;
  rval = pcomm1->filter_pstatus(owned_verts, PSTATUS_NOT_OWNED, PSTATUS_NOT);
  ERRORR(rval, "can't filter for owned vertices only");

  Range entities[4];
  entities[0] = owned_verts;
  entities[2] = quads3;
  // assign new ids only for owned entities
  // will eliminate gaps in global id space for vertices
  rval = pcomm1->assign_global_ids(entities, 2, 1, true, false);
  ERRORR(rval, "can't assign global ids for vertices ");

  /*ErrorCode ParallelComm::assign_global_ids( Range entities[],
   const int dimension,
   const int start_id,
   const bool parallel,
   const bool owned_only) */
#ifndef NDEBUG
  std::stringstream fff;
  fff << "fine0" << pcomm1->proc_config().proc_rank() << ".h5m";
  mb->write_mesh(fff.str().c_str(), &fine_set, 1);
#endif

  rval = mb->write_file("fine.h5m", 0, "PARALLEL=WRITE_PART", &fine_set, 1);
  ERRORR(rval, "can't write set 3, fine ");

  // exchange 1 layer of ghost elements for fine set
  // in the end, fine set will be euler set, I think

  if (7==transport_type)
  {
    rval = pcomm1->exchange_ghost_cells(2, // int ghost_dim
        0, // int bridge_dim
        1, // int num_layers
        0, // int addl_ents
        true); MB_CHK_ERR(rval);
    // these should be used only for ghosting
    rval = mb->get_entities_by_dimension(0, 2, allCells); MB_CHK_ERR(rval);
  }

   // more debugging stuff
#if 0
   fff.clear();
   fff << "fine_gh_" <<  pcomm1->proc_config().proc_rank() << ".h5m";
   EntityHandle ghSet;
   mb->create_meshset(MESHSET_SET, ghSet);
   mb->add_entities(ghSet, allCells);
   mb->write_mesh(fff.str().c_str(), &ghSet, 1);
#endif

  // we need to keep a mapping index, from the coords array to the vertex handles
  // so, for a given vertex entity handle, at what index in the coords array the vertex
  // coordinates are?
  // in the coords array, vertices are repeated 2 times if they are interior to a coarse edge, and
  // repeated 3 or 4 times if they are a corner vertex in a coarse quad

  numVertices = (int) verts.size();
  mapping_to_coords = new int[numVertices];
  for (int k = 0; k < numVertices; k++)
    mapping_to_coords[k] = -1; // it means it was not located yet in vertList
  // vertList is parallel to the coords and dep_coords array

  // now loop over vertsList, and see where
  // vertList has nelem * (nc+1)*(nc+1) vertices; loop over them, and see where are they located

  for (int kk = 0; kk < (int) vertList.size(); kk++) {
    EntityHandle v = vertList[kk];
    int index = verts.index(v);
    if (-1 == index) {
      std::cout << " can't locate vertex " << v << " in vertex Range \n";
      return MB_FAILURE;
    }
    if (mapping_to_coords[index] == -1) // it means the vertex v was not yet encountered in the vertList
        {
      mapping_to_coords[index] = kk;
    }
  }
  // check that every mapping has an index different from -1
  for (int k = 0; k < numVertices; k++) {
    if (mapping_to_coords[k] == -1) {
      {
        std::cout << " vertex at index " << k << " in vertex Range " << verts[k]
            << " is not mapped \n";
        return MB_FAILURE; //
      }
    }
  }

  return rval;
}

// this is called from Fortran 90
void create_mesh(iMesh_Instance instance,
    iBase_EntitySetHandle * imesh_euler_set,
    iBase_EntitySetHandle * imesh_departure_set,
    iBase_EntitySetHandle * imesh_intx_set, double * coords, int * corners,
    int nc, int nelem, int transport_type, MPI_Fint comm, int * ierr) {
  /* double * coords=(double*) icoords;
   int * corners = (int*) icorners;*/
  *ierr = 1;
  Interface * mb = MOABI;
  MPI_Comm mpicomm = MPI_Comm_f2c(comm);
  // instantiate parallel comm now or not?
  pcomm = new ParallelComm(mb, mpicomm);

  EntityHandle coarseSet;
  ErrorCode rval = mb->create_meshset(MESHSET_SET, coarseSet);
  ERRORV(rval, "can't create coarse set ");

  EntityHandle start_vert;
  int totalNumVertices;
  int numCornerVertices;
  std::vector<double *> coordv;
  rval = create_coarse_mesh(mb, pcomm, coarseSet, coords, corners, nc, nelem,
      start_vert, totalNumVertices, numCornerVertices, coordv);
  ERRORV(rval, "can't create coarse set ");

  EntityHandle fine_set;
  rval = mb->create_meshset(MESHSET_SET, fine_set);
  ERRORV(rval, "can't create fine set ");

  rval = create_fine_mesh(mb, pcomm, coarseSet, fine_set, coords, nc, nelem,
      transport_type, start_vert, numCornerVertices, coordv);
  ERRORV(rval, "can't create fine mesh set ");

  // now make euler set the fine set; the order of vertices will be the same after deep copy
  //*imesh_departure_set = (iBase_EntitySetHandle) fine_set;
  *imesh_euler_set = (iBase_EntitySetHandle) fine_set;
  EntityHandle euler_set = fine_set;

  EntityHandle lagr_set;
  rval = mb->create_meshset(MESHSET_SET, lagr_set);
  ERRORV(rval, "can't create moab lagr set ");
  // *imesh_euler_set = (iBase_EntitySetHandle) euler_set;
  *imesh_departure_set = (iBase_EntitySetHandle) lagr_set;

  // call in cslam utils
  // it will copy the second set from the first set
  rval = deep_copy_set_with_quads(mb, fine_set, lagr_set);
  ERRORV(rval, "can't populate lagrange set ");

  EntityHandle intx_set;
  rval = mb->create_meshset(MESHSET_SET, intx_set);
  ERRORV(rval, "can't create output set ");

  *imesh_intx_set = (iBase_EntitySetHandle) intx_set;

  pworker = new Intx2MeshOnSphere(mb);

  pworker->set_box_error(100 * gtol);
  Range local_verts;
  rval = pworker->build_processor_euler_boxes(euler_set, local_verts); // output also the local_verts
  ERRORV(rval, "can't compute euler boxes ");
  pworker->SetErrorTolerance(gtol);

  *ierr = 0;
  return;
}
ErrorCode set_departure_points_position(Interface * mb, EntityHandle lagrSet,
    double * dep_coords, double radius2) {

  // the departure quads are created in the same order as the fine quads
  // for each coarse element, there are nc*nc fine quads
  // their vertices are in order
  Range lagr_quads;
  ErrorCode rval = mb->get_entities_by_type(lagrSet, MBQUAD, lagr_quads);
  ERRORR(rval, "can't get lagrange quads");

  // get all vertices from lagr_quads
  Range lVerts;
  rval = mb->get_connectivity(lagr_quads, lVerts);
  ERRORR(rval, "can't get lagrangian vertices (departure)");

  // they are parallel to the verts Array, they must have the same number of vertices
  assert(numVertices == (int )lVerts.size());

  for (int i = 0; i < numVertices; i++) {
    EntityHandle v = lVerts[i];
    int index = mapping_to_coords[i];
    assert(-1 != index);

    SphereCoords sph;
    sph.R = radius2;
    sph.lat = dep_coords[2 * index];
    sph.lon = dep_coords[2 * index + 1];

    CartVect depPoint = spherical_to_cart(sph);
    rval = mb->set_coords(&v, 1, (double*) depPoint.array());
    ERRORR(rval, "can't set position of vertex");
  }

  return MB_SUCCESS;
}
void intersection_at_level(iMesh_Instance instance,
    iBase_EntitySetHandle fine_set, iBase_EntitySetHandle lagr_set,
    iBase_EntitySetHandle intx_set, double * dep_coords, double radius2,
    int * ierr) {
  *ierr = 1;
  Interface * mb = MOABI;
  //MPI_Comm mpicomm = MPI_Comm_f2c(comm);
  // instantiate parallel comm now or not?

  EntityHandle lagrMeshSet = (EntityHandle) lagr_set;

  //ParallelComm *pcomm = ParallelComm::get_pcomm(mb, 0);
  if (NULL == pcomm)
    return; // error is 1

  // set the departure tag on the fine mesh vertices
  ErrorCode rval = set_departure_points_position(mb, lagrMeshSet, dep_coords,
      radius2);
  if (dprint)
  {
    if (num_update_calls%nlev==level_to_print-1)
    {
      int nstep=num_update_calls/nlev;
      std::stringstream fff;
      fff << "depMesh0" << nstep << ".h5m";

      rval = mb->write_file(fff.str().c_str(), 0, "PARALLEL=WRITE_PART", &lagrMeshSet, 1); MB_CHK_ERR_RET(rval);
    }
  }
  ERRORV(rval, "can't set departure tag");
  if (debug) {
    std::stringstream fff;
    fff << "lagr0" << pcomm->proc_config().proc_rank() << ".vtk";
    rval = mb->write_mesh(fff.str().c_str(), &lagrMeshSet, 1);
    ERRORV(rval, "can't write covering set ");
  }

  // it should be done earlier
  pworker->SetRadius(radius);

  EntityHandle covering_set;
  rval = pworker->create_departure_mesh_3rd_alg(lagrMeshSet, covering_set);
  ERRORV(rval, "can't compute covering set ");

  if (debug) {
    std::stringstream fff;
    fff << "cover" << pcomm->proc_config().proc_rank() << ".vtk";
    rval = mb->write_mesh(fff.str().c_str(), &covering_set, 1);

    ERRORV(rval, "can't write covering set ");
  }
  EntityHandle intxSet = (EntityHandle) intx_set;
  rval = pworker->intersect_meshes(covering_set, (EntityHandle) fine_set,
      intxSet);
  ERRORV(rval, "can't intersect ");

  if (debug) {
    std::stringstream fff;
    fff << "intx0" << pcomm->proc_config().proc_rank() << ".vtk";
    rval = mb->write_mesh(fff.str().c_str(), &intxSet, 1);
    ERRORV(rval, "can't write covering set ");
  }

  return;
}
void cleanup_after_intersection(iMesh_Instance instance,
    iBase_EntitySetHandle lagr_set, iBase_EntitySetHandle euler_set,
    iBase_EntitySetHandle intx_set, int * ierr) {
  *ierr = 1;
  Interface * mb = MOABI;
  // delete elements
  // delete now the polygons and the elements of out_set
  // also, all verts that are not in euler set or lagr_set
  Range allVerts;
  ErrorCode rval = mb->get_entities_by_dimension(0, 0, allVerts);
  ERRORV(rval, "can't get all vertices");

  Range allElems;
  rval = mb->get_entities_by_dimension(0, 2, allElems);
  ERRORV(rval, "can't get all elems");

  Range polys=allCells; // first euler polys and ghosts, if created
  rval = mb->get_entities_by_dimension((EntityHandle) euler_set, 2, polys); // redundant, these are already here, from allCells FIXME

  // add to polys range the lagr polys
  rval = mb->get_entities_by_dimension((EntityHandle) lagr_set, 2, polys); // do not delete lagr set either, with its vertices
  ERRORV(rval, "can't get all polys from lagr set");
  // add to the connecVerts range all verts, from all initial polys
  Range vertsToStay;
  rval = mb->get_connectivity(polys, vertsToStay);
  ERRORV(rval, "get verts that stay");

  Range todeleteVerts = subtract(allVerts, vertsToStay);

  Range todeleteElem = subtract(allElems, polys); // this is coarse mesh too (if still here)

  // empty the out mesh set
  EntityHandle out_set = (EntityHandle) intx_set;
  rval = mb->clear_meshset(&out_set, 1);
  ERRORV(rval, "clear mesh set");

  rval = mb->delete_entities(todeleteElem);
  ERRORV(rval, "delete intx elements");
  rval = mb->delete_entities(todeleteVerts);
  ERRORV(rval, "failed to delete intx vertices");

  *ierr = 0;
  return;
}

void cleanup_after_simulation(int * ierr) {
  delete[] mapping_to_coords;
  numVertices = 0;
  *ierr = 0;
  return;
}
#ifdef __cplusplus
} // extern "C"
#endif

