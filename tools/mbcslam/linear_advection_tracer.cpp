
#include <iostream>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "Intx2MeshOnSphere.hpp"
#include <math.h>
#include "moab/ParallelComm.hpp"
#include "moab/ProgOptions.hpp"
#include "MBParallelConventions.h"
#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"
#include "TestUtil.hpp"
#include "CslamUtils.hpp"

std::string file_name("./uniform_30.g");
//std::string file_name("./uniform_120.g");

#ifdef MESHDIR
std::string TestDir( STRINGIFY(MESHDIR) );
#else
std::string TestDir(".");
#endif

using namespace moab;

moab::ErrorCode update_density(moab::Interface * mb, moab::EntityHandle euler_set,
    moab::EntityHandle lagr_set, moab::EntityHandle out_set, moab::Tag & rhoTag, 
    moab::Tag & areaTag, moab::Tag & rhoCoefTag, moab::Tag & weightsTag, moab::Tag & planeTag);

moab::ErrorCode update_tracer(moab::Interface * mb, moab::EntityHandle euler_set,
    moab::EntityHandle lagr_set, moab::EntityHandle out_set, moab::Tag & tauTag, 
    moab::Tag & areaTag, moab::Tag & rhoCoefTag, moab::Tag & tauCoefTag, moab::Tag & weightsTag, 
    moab::Tag & planeTag);

// covering set is output now (!)
moab::ErrorCode get_departure_grid(moab::Interface * mb, moab::EntityHandle euler_set,
                                   moab::EntityHandle lagr_set, moab::EntityHandle & covering_set, int tStep, moab::Range & connecVerts);

void get_barycenters(moab::Interface * mb, moab::EntityHandle set, moab::Tag &planeTag, moab::Tag &cellIntTag,
                     moab::Tag &barycenterTag);

void get_centers_of_mass(moab::Interface * mb, moab::EntityHandle set, moab::Tag &planeTag, moab::Tag &rhoCoefTag,
                         moab::Tag &cellIntTag, moab::Tag &centerOfMassTag);

void get_eul_cell_integrals(moab::Interface * mb, moab::EntityHandle set, moab::Tag &planeTag, moab::Tag &areaTag, moab::Tag &cellIntTag);

void get_gnomonic_plane(moab::Interface * mb, moab::EntityHandle set, moab::Tag &planeTag);

void get_linear_reconstruction(moab::Interface * mb, moab::EntityHandle set, moab::Tag &cellValTag, moab::Tag &planeTag, 
                               moab::Tag &centerTag, moab::Tag &linearCoefTag);

void get_neighborhood_bounds(moab::Interface * mb, moab::EntityHandle set, moab::Tag &cellValTag, moab::Tag &boundsTag); 

void limit_linear_reconstruction(moab::Interface * mb, moab::EntityHandle set, moab::Tag &cellValTag, moab::Tag &planeTag, 
                                 moab::Tag &boundsTag, moab::Tag &centerTag, moab::Tag &linearCoefTag);

void get_intersection_weights(moab::Interface * mb, moab::EntityHandle euler_set, moab::EntityHandle lagr_set, 
                              moab::EntityHandle intx_set, moab::Tag &planeTag, moab::Tag &weightsTag);

void set_initial_values(moab::Interface * mb, moab::EntityHandle euler_set, moab::Tag & centerTag, 
                        moab::Tag & valTag, int field_type);

// functions to compute departure point locations
void departure_point_swirl(moab::CartVect & arrival_point, double t, double delta_t, moab::CartVect & departure_point);
void departure_point_swirl_rot(moab::CartVect & arrival_point, double t, double delta_t, moab::CartVect & departure_point);
void departure_point_rotation(moab::CartVect & arrival_point, double t, double delta_t, moab::CartVect & departure_point);


double gtol = 1.e-9; // this is for geometry tolerance

double radius = 1.;

bool writeFiles = true;
bool parallelWrite = false;
bool velocity = false;

int numSteps = 100; // number of times with velocity displayed at points
double T = 5;

Intx2MeshOnSphere * pworker = NULL;

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);

  // set up MOAB interface and parallel communication
  moab::Core moab;
  moab::Interface& mb = moab;
  moab::ParallelComm pcomm(&mb, MPI_COMM_WORLD);

  //int rank = mb_pcomm->proc_config().proc_rank();
  int rank = pcomm.proc_config().proc_rank();

  // create meshset
  moab::EntityHandle euler_set;
  moab::ErrorCode rval = mb.create_meshset(MESHSET_SET, euler_set);

  std::stringstream opts;
  //opts << "PARALLEL=READ_PART;PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;GATHER_SET=0;PARTITION_METHOD=TRIVIAL_PARTITION;VARIABLE=";
  opts << "PARALLEL=READ_PART;PARTITION;PARALLEL_RESOLVE_SHARED_ENTS";
  //rval = mb.load_file(file_name.c_str(), &euler_set, opts.str().c_str());
  // read a homme file, partitioned in 4 so far
  std::string fileN = TestDir + "/mbcslam/fine4.h5m";
  const char *filename_mesh1 = fileN.c_str();

  rval = mb.load_file(filename_mesh1, &euler_set, opts.str().c_str());  MB_CHK_ERR(rval);
  /*int num_entities;
   rval = mb.get_number_entities_by_dimension(euler_set, 2, num_entities);MB_CHK_ERR(rval);*/

  rval = pcomm.exchange_ghost_cells(2, // int ghost_dim
      0, // int bridge_dim
      1, // int num_layers
      0, // int addl_ents
      true); MB_CHK_ERR(rval); // bool store_remote_handles

  // check if euler set is changed

  /*   rval = mb.get_number_entities_by_dimension(euler_set, 2, num_entities);MB_CHK_ERR(rval);*/
  // so indeed, euler set is not modified; ghost elements are created, and we verified below
  // Create tag for cell density

  moab::Range allCells; /// these will be all cells on current task, used for tag exchange
  rval = mb.get_entities_by_dimension(0, 2, allCells); MB_CHK_ERR(rval);
  moab::Tag rhoTag = 0;
  rval = mb.tag_get_handle("Density", 1, moab::MB_TYPE_DOUBLE, rhoTag,
      moab::MB_TAG_CREAT | moab::MB_TAG_DENSE); MB_CHK_ERR(rval);

  // Create tag for cell tracer
  moab::Tag tauTag = 0;
  rval = mb.tag_get_handle("Tracer", 1, moab::MB_TYPE_DOUBLE, tauTag,
      moab::MB_TAG_CREAT | moab::MB_TAG_DENSE); MB_CHK_ERR(rval);

  // Create tag for cell area
  moab::Tag areaTag = 0;
  rval = mb.tag_get_handle("Area", 1, moab::MB_TYPE_DOUBLE, areaTag,
      moab::MB_TAG_CREAT | moab::MB_TAG_DENSE); MB_CHK_ERR(rval);

  // Create tag for cell barycenters in 3D Cartesian space
  moab::Tag barycenterTag = 0;
  rval = mb.tag_get_handle("CellBarycenter", 3, moab::MB_TYPE_DOUBLE,
      barycenterTag, moab::MB_TAG_CREAT | moab::MB_TAG_DENSE); MB_CHK_ERR(rval);

  // Create tag for cell centers of mass in 3D Cartesian space
  moab::Tag centerOfMassTag = 0;
  rval = mb.tag_get_handle("CellCenterOfMass", 3, moab::MB_TYPE_DOUBLE,
      centerOfMassTag, moab::MB_TAG_CREAT | moab::MB_TAG_DENSE); MB_CHK_ERR(rval);

  // Create tag for integrals of (x, y, x^2, y^2, xy) over Eulerian cells
  moab::Tag cellIntTag = 0;
  rval = mb.tag_get_handle("CellIntegral", 6, moab::MB_TYPE_DOUBLE, cellIntTag,
      moab::MB_TAG_CREAT | moab::MB_TAG_DENSE); MB_CHK_ERR(rval);

  // Create tag for cell density reconstruction coefficients
  moab::Tag rhoCoefTag = 0;
  rval = mb.tag_get_handle("LinearCoefRho", 3, moab::MB_TYPE_DOUBLE, rhoCoefTag,
      moab::MB_TAG_CREAT | moab::MB_TAG_DENSE); MB_CHK_ERR(rval);

  // Create tag for cell tracer reconstruction coefficients
  moab::Tag tauCoefTag = 0;
  rval = mb.tag_get_handle("LinearCoefTau", 3, moab::MB_TYPE_DOUBLE, tauCoefTag,
      moab::MB_TAG_CREAT | moab::MB_TAG_DENSE); MB_CHK_ERR(rval);

  // Create tag for index of gnomonic plane for each cell
  moab::Tag planeTag = 0;
  rval = mb.tag_get_handle("GnomonicPlane", 1, moab::MB_TYPE_INTEGER, planeTag,
      moab::MB_TAG_CREAT | moab::MB_TAG_DENSE); MB_CHK_ERR(rval);

  // Create tag for density bounds
  moab::Tag rhoBoundsTag = 0;
  rval = mb.tag_get_handle("DensityBounds", 2, moab::MB_TYPE_DOUBLE,
      rhoBoundsTag, moab::MB_TAG_CREAT | moab::MB_TAG_DENSE); MB_CHK_ERR(rval);

  // Create tag for tracer bounds
  moab::Tag tauBoundsTag = 0;
  rval = mb.tag_get_handle("TracerBounds", 2, moab::MB_TYPE_DOUBLE,
      tauBoundsTag, moab::MB_TAG_CREAT | moab::MB_TAG_DENSE); MB_CHK_ERR(rval);

  // Create tag for intersection weights
  moab::Tag weightsTag = 0;
  rval = mb.tag_get_handle("Weights", 6, moab::MB_TYPE_DOUBLE, weightsTag,
      moab::MB_TAG_CREAT | moab::MB_TAG_DENSE); MB_CHK_ERR(rval);

  // get cell plane
  get_gnomonic_plane(&mb, euler_set, planeTag);

  // get Eulerian cell integrals over x, y, x^2, y^2, xy (and cell area)
  get_eul_cell_integrals(&mb, euler_set, planeTag, areaTag, cellIntTag);

  // get cell barycenters
  get_barycenters(&mb, euler_set, planeTag, cellIntTag, barycenterTag);

  // Set density distribution
  set_initial_values(&mb, euler_set, barycenterTag, rhoTag, 4);

  // Set tracer distribution
  // (1 - cosine bells, 2 - Gaussian hills, 3 - notched cylinders, 4 - constant)
  set_initial_values(&mb, euler_set, barycenterTag, tauTag, 2);

  // Get initial values for use in error computation
  moab::Range redEls;
  rval = mb.get_entities_by_dimension(euler_set, 2, redEls); MB_CHK_ERR(rval);
  std::vector<double> iniValsRho(redEls.size());
  rval = mb.tag_get_data(rhoTag, redEls, &iniValsRho[0]); MB_CHK_ERR(rval);
  std::vector<double> iniValsTau(redEls.size());
  rval = mb.tag_get_data(tauTag, redEls, &iniValsTau[0]); MB_CHK_ERR(rval);

  // Get Lagrangian set
  moab::EntityHandle out_set, lagrange_set, covering_set; // covering set created again; maybe it should be created only once
  rval = mb.create_meshset(MESHSET_SET, out_set); MB_CHK_ERR(rval);
  rval = mb.create_meshset(MESHSET_SET, lagrange_set); MB_CHK_ERR(rval);
  //rval = mb.create_meshset(MESHSET_SET, covering_set); MB_CHK_ERR(rval);

  rval = deep_copy_set(&mb, euler_set, lagrange_set); MB_CHK_ERR(rval);
  moab::EntityHandle dum = 0;
  moab::Tag corrTag;
  rval = mb.tag_get_handle(CORRTAGNAME, 1, MB_TYPE_HANDLE, corrTag,
      MB_TAG_DENSE | MB_TAG_CREAT, &dum); MB_CHK_ERR(rval);

  //Set up intersection of two meshes
  pworker = new Intx2MeshOnSphere(&mb);
  pworker->SetErrorTolerance(gtol);
  pworker->SetRadius(radius);
  pworker->set_box_error(100 * gtol);

  // these stay fixed for one run
  moab::Range local_verts;
  rval = pworker->build_processor_euler_boxes(euler_set, local_verts); MB_CHK_ERR(rval);// output also the local_verts

  // density, tracer tags will be used to compute the linear reconstruction using neighbor information
  // for that, we will need those values on the ghost cells too
  // exchange tags for ghosted elements
  std::vector<moab::Tag> tags;
  tags.push_back(rhoTag); tags.push_back(tauTag); tags.push_back(barycenterTag);

  rval = pcomm.exchange_tags(tags, tags, allCells); MB_CHK_ERR(rval);

  // write the part
  std::stringstream meshFile;
  meshFile << "mesh_" << rank << ".vtk";
  rval = mb.write_file(meshFile.str().c_str()); MB_CHK_ERR(rval);

  // loop over time to update density
  for (int ts = 1; ts < numSteps + 1; ts++) {

    if (ts == 1)  // output initial condition
        {
      std::stringstream newTracer;
      newTracer << "Tracer" << rank << "_" << ts - 1 << ".vtk";
      rval = mb.write_file(newTracer.str().c_str(), 0, 0, &euler_set, 1);
    }

    // get density bounds
    get_neighborhood_bounds(&mb, euler_set, rhoTag, rhoBoundsTag);

    // get tracer bounds
    get_neighborhood_bounds(&mb, euler_set, tauTag, tauBoundsTag);

    // get linear reconstruction coefficients for density
    get_linear_reconstruction(&mb, euler_set, rhoTag, planeTag, barycenterTag,
        rhoCoefTag);

    // limit linear reconstruction coefficients for density
    limit_linear_reconstruction(&mb, euler_set, rhoTag, planeTag, rhoBoundsTag,
        barycenterTag, rhoCoefTag);

    // get cell centers of mass
    get_centers_of_mass(&mb, euler_set, planeTag, rhoCoefTag, cellIntTag,
        centerOfMassTag);

    // need to communicate center of mass tag for ghost elements, too
    // it is needed for linear reconstruction of tracer coeffs
    rval = pcomm.exchange_tags(centerOfMassTag, allCells); MB_CHK_ERR(rval);
    // get linear reconstruction coefficients for tracer
    get_linear_reconstruction(&mb, euler_set, tauTag, planeTag, centerOfMassTag,
        tauCoefTag);

    // limit linear reconstruction coefficients for tracer
    limit_linear_reconstruction(&mb, euler_set, tauTag, planeTag, tauBoundsTag,
        centerOfMassTag, tauCoefTag);

    // get depature grid
    rval = get_departure_grid(&mb, euler_set, lagrange_set, covering_set, ts,
        local_verts);

    if (writeFiles) // so if write
    {
      std::stringstream coverFile;
      coverFile << "cover" << rank << "_" << ts << ".vtk";
      rval = mb.write_file(coverFile.str().c_str(), 0, 0, &covering_set, 1);
    }

    // intersect the meshes
    rval = pworker->intersect_meshes(covering_set, euler_set, out_set);

    // intersection weights (i.e. area, x integral, and y integral over cell intersections)
    get_intersection_weights(&mb, euler_set, covering_set, out_set, planeTag,
        weightsTag);

    // update the density
    rval = update_density(&mb, euler_set, lagrange_set, out_set, rhoTag,
        areaTag, rhoCoefTag, weightsTag, planeTag);

    // update the tracer
    rval = update_tracer(&mb, euler_set, lagrange_set, out_set, tauTag, areaTag,
        rhoCoefTag, tauCoefTag, weightsTag, planeTag);

    if (writeFiles && (ts % 2 == 0)) // so if write
        {
      std::stringstream newTracer;
      newTracer << "Tracer" << rank << "_" << ts << ".vtk";
      rval = mb.write_file(newTracer.str().c_str(), 0, 0, &euler_set, 1);
    }

    // delete the polygons and elements of out_set
    moab::Range allVerts;
    rval = mb.get_entities_by_dimension(0, 0, allVerts);

    moab::Range allElems;
    rval = mb.get_entities_by_dimension(0, 2, allElems);

    // get Eulerian and lagrangian cells
    moab::Range polys;
    rval = mb.get_entities_by_dimension(euler_set, 2, polys);
    rval = mb.get_entities_by_dimension(lagrange_set, 2, polys); // do not delete lagr set either, with its vertices

    // add to the connecVerts range all verts, from all initial polys
    moab::Range vertsToStay;
    rval = mb.get_connectivity(polys, vertsToStay);

    moab::Range todeleteVerts = subtract(allVerts, vertsToStay);

    moab::Range todeleteElem = subtract(allElems, polys);
    // empty the out mesh set
    rval = mb.clear_meshset(&out_set, 1);

    rval = mb.delete_entities(todeleteElem);
    rval = mb.delete_entities(todeleteVerts);
    if (rank == 0)
      std::cout << " step: " << ts << "\n";

  }

  //final vals and errors
  moab::Range::iterator iter = redEls.begin();
  double norm1rho = 0.;
  double norm2rho = 0.;
  double exact2rho = 0.;
  double exact1rho = 0.;
  double norm1tau = 0.;
  double norm2tau = 0.;
  double exact2tau = 0.;
  double exact1tau = 0.;
  int count = 0;
  void * data;
  int j = 0; // index in iniVals
  while (iter != redEls.end()) {
    rval = mb.tag_iterate(rhoTag, iter, redEls.end(), count, data);
    double * ptrDensity = (double*) data;

    rval = mb.tag_iterate(tauTag, iter, redEls.end(), count, data);
    double * ptrTracer = (double*) data;

    rval = mb.tag_iterate(areaTag, iter, redEls.end(), count, data);
    double * ptrArea = (double*) data;
    for (int i = 0; i < count; i++, iter++, ptrTracer++, ptrArea++, j++) {
      //double area = *ptrArea;
      norm1rho += fabs(*ptrDensity - iniValsRho[j]) * (*ptrArea);
      norm2rho += (*ptrDensity - iniValsRho[j]) * (*ptrDensity - iniValsRho[j])
          * (*ptrArea);
      norm1tau += fabs(*ptrTracer - iniValsTau[j]) * (*ptrArea);
      norm2tau += (*ptrTracer - iniValsTau[j]) * (*ptrTracer - iniValsTau[j])
          * (*ptrArea);
      exact1rho += (iniValsRho[j]) * (*ptrArea);
      exact2rho += (iniValsRho[j]) * (iniValsRho[j]) * (*ptrArea);
      exact1tau += (iniValsTau[j]) * (*ptrArea);
      exact2tau += (iniValsTau[j]) * (iniValsTau[j]) * (*ptrArea);
    }
  }

  double total_norm1rho = 0;
  double total_norm2rho = 0;
  double total_exact1rho = 0;
  double total_exact2rho = 0;
  double total_norm1tau = 0;
  double total_norm2tau = 0;
  double total_exact1tau = 0;
  double total_exact2tau = 0;
  int mpi_err = MPI_Reduce(&norm1rho, &total_norm1rho, 1, MPI_DOUBLE, MPI_SUM,
      0, MPI_COMM_WORLD);
  mpi_err = MPI_Reduce(&norm2rho, &total_norm2rho, 1, MPI_DOUBLE, MPI_SUM, 0,
      MPI_COMM_WORLD);
  mpi_err = MPI_Reduce(&exact1rho, &total_exact1rho, 1, MPI_DOUBLE, MPI_SUM, 0,
      MPI_COMM_WORLD);
  mpi_err = MPI_Reduce(&exact2rho, &total_exact2rho, 1, MPI_DOUBLE, MPI_SUM, 0,
      MPI_COMM_WORLD);
  mpi_err = MPI_Reduce(&norm1tau, &total_norm1tau, 1, MPI_DOUBLE, MPI_SUM, 0,
      MPI_COMM_WORLD);
  mpi_err = MPI_Reduce(&norm2tau, &total_norm2tau, 1, MPI_DOUBLE, MPI_SUM, 0,
      MPI_COMM_WORLD);
  mpi_err = MPI_Reduce(&exact1tau, &total_exact1tau, 1, MPI_DOUBLE, MPI_SUM, 0,
      MPI_COMM_WORLD);
  mpi_err = MPI_Reduce(&exact2tau, &total_exact2tau, 1, MPI_DOUBLE, MPI_SUM, 0,
      MPI_COMM_WORLD);
  if (0 == rank)
    std::cout << " numSteps: " << numSteps << " 1-norm rho:"
        << total_norm1rho / total_exact1rho << " 2-norm rho:"
        << std::sqrt(total_norm2rho / total_exact2rho) << "\n";
  std::cout << "                1-norm tau:" << total_norm1tau / total_exact1tau
      << " 2-norm tau:" << std::sqrt(total_norm2tau / total_exact2tau) << "\n";

  MPI_Finalize();
  return 0;
}

void get_gnomonic_plane(moab::Interface * mb, moab::EntityHandle set, moab::Tag &planeTag)
{
  // get all entities of dimension 2
  moab::Range cells;
  ErrorCode rval = mb->get_entities_by_dimension(set, 2, cells);
  if (MB_SUCCESS != rval)
    return;

  for (Range::iterator it = cells.begin(); it != cells.end(); it++)
  {
    moab::EntityHandle icell = *it;

    // get the nodes
    const moab::EntityHandle * verts;
    int num_nodes;
    rval = mb->get_connectivity(icell, verts, num_nodes);
    if (MB_SUCCESS != rval)
      return;

    // get coordinates
    std::vector<double> coords(3 * num_nodes);
    rval = mb->get_coords(verts, num_nodes, &coords[0]);
    if (MB_SUCCESS != rval)
      return;

    // get cell center
     double centerx = 0;
     double centery = 0;
     double centerz = 0;
     for (int inode = 0; inode < num_nodes; inode++){
        centerx += coords[inode*3]/num_nodes;
        centery += coords[inode*3+1]/num_nodes;
        centerz += coords[inode*3+2]/num_nodes;
     }
     double rad = std::sqrt(centerx*centerx + centery*centery + centerz*centerz);
     centerx = centerx/rad;
     centery = centery/rad;
     centerz = centerz/rad;
     moab::CartVect center(centerx,centery,centerz);

    // define gnomonic plane based on cell center coordinates
     int plane = 0;
     decide_gnomonic_plane(center,plane);

     rval = mb->tag_set_data(planeTag, &icell, 1, &plane);
  }
  return;
}


void get_barycenters(moab::Interface * mb, moab::EntityHandle set, moab::Tag &planeTag, 
                     moab::Tag &cellIntTag, moab::Tag &barycenterTag)
{

  // get all entities of dimension 2
  moab::Range cells;
  ErrorCode rval = mb->get_entities_by_dimension(set, 2, cells);
  if (MB_SUCCESS != rval)
    return;

  // set sphere radius to 1
   double R = 1.0;

  for (Range::iterator it = cells.begin(); it != cells.end(); it++)
  {
    moab::EntityHandle icell = *it;

    // get the nodes
    const moab::EntityHandle * verts;
    int num_nodes;
    rval = mb->get_connectivity(icell, verts, num_nodes);
    if (MB_SUCCESS != rval)
      return;

    // get gnomonic plane 
     int plane = 0;
     rval = mb->tag_get_data(planeTag, &icell, 1, &plane );
     if (MB_SUCCESS != rval)
       return;

    // get cell integral values
     std::vector<double> cellInt(6);
     rval = mb->tag_get_data(cellIntTag, &icell, 1, &cellInt[0]);
     if (MB_SUCCESS != rval)
       return;

    // barycenter in gnomonic coordinates
     double bary_x = 0;
     double bary_y = 0;
     bary_x = cellInt[1] / cellInt[0];
     bary_y = cellInt[2] / cellInt[0];

    // barycenter in Cartesian X,Y,Z coordinates
     moab::CartVect barycent;
     reverse_gnomonic_projection(bary_x, bary_y, R, plane, barycent);

    // set barycenter
     std::vector<double> barycenter(3);
     barycenter[0] = barycent[0];
     barycenter[1] = barycent[1];
     barycenter[2] = barycent[2];
     rval = mb->tag_set_data(barycenterTag, &icell, 1, &barycenter[0]);

  }
  return;
}

void get_eul_cell_integrals(moab::Interface * mb, moab::EntityHandle set, moab::Tag &planeTag, 
                            moab::Tag &areaTag, moab::Tag &cellIntTag)
{
  // get all entities of dimension 2
  moab::Range cells;
  ErrorCode rval = mb->get_entities_by_dimension(set, 2, cells);
  if (MB_SUCCESS != rval)
    return;

  // set sphere radius to 1
   double R = 1.0;

  double total_area = 0;
  for (Range::iterator it = cells.begin(); it != cells.end(); it++)
  {
    moab::EntityHandle icell = *it;

    // get the nodes
    const moab::EntityHandle * verts;
    int num_nodes;
    rval = mb->get_connectivity(icell, verts, num_nodes);
    if (MB_SUCCESS != rval)
      return;

    // get coordinates
    std::vector<double> coords(3 * num_nodes);
    rval = mb->get_coords(verts, num_nodes, &coords[0]);
    if (MB_SUCCESS != rval)
      return;

    // get gnomonic plane 
     int plane = 0;
     rval = mb->tag_get_data(planeTag, &icell, 1, &plane );
     if (MB_SUCCESS != rval)
       return;

    // get vertex coordinates and project onto gnomonic plane
     std::vector<double> x(num_nodes);
     std::vector<double> y(num_nodes);
     double area   = 0;
     double int_x  = 0;
     double int_y  = 0;
     double int_x2 = 0;
     double int_y2 = 0;
     double int_xy = 0;
     for (int inode = 0; inode < num_nodes; inode++){
         double rad = sqrt(coords[inode*3]*coords[inode*3] + coords[inode*3+1]*coords[inode*3+1] + coords[inode*3+2]*coords[inode*3+2]);
         CartVect xyzcoord(coords[inode*3]/rad,coords[inode*3+1]/rad,coords[inode*3+2]/rad);
         gnomonic_projection(xyzcoord, R, plane, x[inode],y[inode]);
        
     }

    // integrate over cell 
     for (int inode = 0; inode < num_nodes; inode++){
         int inode2 = inode+1;
         if (inode2 >= num_nodes) inode2 = 0;
         double xmid = 0.5*(x[inode] + x[inode2]);
         double ymid = 0.5*(y[inode] + y[inode2]);
         double r1 = sqrt(1 + x[inode]*x[inode] + y[inode]*y[inode]);
         double rm = sqrt(1 + xmid*xmid + ymid*ymid);
         double r2 = sqrt(1 + x[inode2]*x[inode2] + y[inode2]*y[inode2]);
         double hx=x[inode2]-x[inode];

         // int dA
         area  += -hx*(y[inode]/(r1*(1 + x[inode]*x[inode]))
                            + 4.0*ymid/(rm*(1+xmid*xmid))
                            + y[inode2]/(r2*(1 + x[inode2]*x[inode2])))/6.0;

         // int x dA
         int_x  += -hx*(x[inode]*y[inode]/(r1*(1 + x[inode]*x[inode]))
                            + 4.0*xmid*ymid/(rm*(1+xmid*xmid))
                            + x[inode2]*y[inode2]/(r2*(1 + x[inode2]*x[inode2])))/6.0;

         // int y dA
         int_y  += hx*(1.0/r1 + 4.0/rm + 1.0/r2)/6.0;

         // int x^2 dA
         int_x2 += -hx*(x[inode]*x[inode]*y[inode]/(r1*(1.0 + x[inode]*x[inode]))
                            + 4.0*xmid*xmid*ymid/(rm*(1+xmid*xmid))
                            + x[inode2]*x[inode2]*y[inode2]/(r2*(1 + x[inode2]*x[inode2])))/6.0;

         // int y^2 dA
         int_y2  += hx*((y[inode]/r1 - asinh(y[inode]/sqrt(1.0 + x[inode]*x[inode])))
                          + 4.0*(ymid/rm - asinh(ymid/sqrt(1.0 + xmid*xmid)))
                          + (y[inode2]/r2 - asinh(y[inode2]/sqrt(1.0 + x[inode2]*x[inode2]))))/6.0;

         // int xy dA
         int_xy  += hx*(x[inode]/r1 + 4.0*xmid/rm + x[inode2]/r2)/6.0;
     }

     total_area+=area;

    // set integrals
     std::vector<double> ints(6);
     ints[0] = area;
     ints[1] = int_x;
     ints[2] = int_y;
     ints[3] = int_x2;
     ints[4] = int_y2;
     ints[5] = int_xy;
     rval = mb->tag_set_data(cellIntTag, &icell, 1, &ints[0]);

     // set area
     rval = mb->tag_set_data(areaTag, &icell, 1, &area);
  }

   std::cout << "total Eulerian cell area = " << total_area << "\n";

  return;
}

void get_centers_of_mass(moab::Interface * mb, moab::EntityHandle set, moab::Tag &planeTag, moab::Tag &rhoCoefTag,
                         moab::Tag &cellIntTag, moab::Tag &centerOfMassTag)
{
  // get all entities of dimension 2
  moab::Range cells;
  ErrorCode rval = mb->get_entities_by_dimension(set, 2, cells);
  if (MB_SUCCESS != rval)
    return;

  // set sphere radius to 1
   double R = 1.0;

  for (Range::iterator it = cells.begin(); it != cells.end(); it++)
  {
    moab::EntityHandle icell = *it;

    // get the nodes
    const moab::EntityHandle * verts;
    int num_nodes;
    rval = mb->get_connectivity(icell, verts, num_nodes);
    if (MB_SUCCESS != rval)
      return;

    // get gnomonic plane 
     int plane = 0;
     rval = mb->tag_get_data(planeTag, &icell, 1, &plane );
     if (MB_SUCCESS != rval)
       return;

    // get Eulerian cell integrals
     std::vector<double> cellInt(6);
     rval = mb->tag_get_data(cellIntTag, &icell, 1, &cellInt[0]);
     if (MB_SUCCESS != rval)
       return;

    // get linear rho coefficients
     std::vector<double> rhoCoefs(3);
     rval = mb->tag_get_data(rhoCoefTag, &icell, 1, &rhoCoefs[0]);
     if (MB_SUCCESS != rval)
       return;

     // get cell mass
     double mass = rhoCoefs[0]*cellInt[0] + rhoCoefs[1]*cellInt[1] + rhoCoefs[2]*cellInt[2];

     // center of mass defaults to barycenter
     double com_x = cellInt[1]/cellInt[0];
     double com_y = cellInt[2]/cellInt[0];

     // if mass is nonzero, compute actual center of mass
     if (mass > 1.0e-10) {
        // com_x = \int x rho^h
        com_x = rhoCoefs[0]*cellInt[1] + rhoCoefs[1]*cellInt[3] + rhoCoefs[2]*cellInt[5];
        // com_y = \int y rho^h
        com_y = rhoCoefs[0]*cellInt[2] + rhoCoefs[1]*cellInt[5] + rhoCoefs[2]*cellInt[4];
        com_x = com_x/mass;
        com_y = com_y/mass;
     }

    // center of mass in Cartesian X,Y,Z coordinates
     moab::CartVect com_XYZ;
     reverse_gnomonic_projection(com_x, com_y, R, plane, com_XYZ);

    // set center of mass
     std::vector<double> center_of_mass(3);
     center_of_mass[0] = com_XYZ[0];
     center_of_mass[1] = com_XYZ[1];
     center_of_mass[2] = com_XYZ[2];
     rval = mb->tag_set_data(centerOfMassTag, &icell, 1, &center_of_mass[0]);

  }
  return;
}

void set_initial_values(moab::Interface * mb, moab::EntityHandle euler_set, moab::Tag & barycenterTag, 
                        moab::Tag & rhoTag, int field_type)
{
  moab::ErrorCode rval = MB_SUCCESS;

  // get cells
  moab::Range cells;
  rval = mb->get_entities_by_dimension(euler_set, 2, cells);
  if (MB_SUCCESS != rval)
    return;

  // get barycenters
  std::vector<double> cell_barys(3*cells.size());
  rval = mb->tag_get_data(barycenterTag, cells, &cell_barys[0]);

  // loop over cells 
  int cell_ind = 0;
  for (Range::iterator it = cells.begin(); it != cells.end(); it++)
  {
    moab::EntityHandle icell = *it;
   
     // convert barycenter from 3-D Cartesian to lat/lon
      moab::CartVect bary_xyz(cell_barys[cell_ind*3],cell_barys[cell_ind*3+1],cell_barys[cell_ind*3+2]);
      moab::SphereCoords sphCoord = cart_to_spherical(bary_xyz);

      if (field_type == 1)  // cosine bells
      {
        //                 lon1,        lat1  lon2    lat2   b    c  hmax  r
        double params[] = { 5*M_PI/6.0, 0.0, 7*M_PI/6, 0.0, 0.1, 0.9, 1., 0.5};

        double rho_barycent = quasi_smooth_field(sphCoord.lon, sphCoord.lat, params);
        rval = mb->tag_set_data(rhoTag, &icell, 1, &rho_barycent);
      }

      if (field_type == 2)  // Gaussian hills
      {
        moab::CartVect p1, p2;
        moab::SphereCoords spr;
        spr.R = 1;
        //spr.lat = M_PI/3;
        //spr.lon= M_PI;
        spr.lat = 0.0;
        spr.lon= 5*M_PI/6.0;
        p1 = spherical_to_cart(spr);
        //spr.lat = -M_PI/3;
        spr.lon= 7*M_PI/6.0;
        p2 = spherical_to_cart(spr);
        // X1, Y1, Z1, X2, Y2, Z2, ?, ?
        double params[] = { p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], 1,    5.};

        double rho_barycent = smooth_field(sphCoord.lon, sphCoord.lat, params);
        rval = mb->tag_set_data(rhoTag, &icell, 1, &rho_barycent);
      }

      if (field_type == 3)  // Zalesak cylinders
      {
        //                   lon1,      lat1,    lon2,   lat2, b,   c,   r
        double params[] = { 5*M_PI/6.0, 0.0, 7*M_PI/6.0, 0.0, 0.1, 0.9, 0.5};

        double rho_barycent = slotted_cylinder_field(sphCoord.lon, sphCoord.lat, params);
        rval = mb->tag_set_data(rhoTag, &icell, 1, &rho_barycent);
      }

      if (field_type == 4)  // constant
      {
        double rho_barycent = 1.0;
        rval = mb->tag_set_data(rhoTag, &icell, 1, &rho_barycent);
      }

     cell_ind++;
  }

}

void get_linear_reconstruction(moab::Interface * mb, moab::EntityHandle set, moab::Tag &cellValTag, 
                               moab::Tag &planeTag, moab::Tag &centerTag, moab::Tag &linearCoefTag)
{
  // get all entities of dimension 2
  Range cells;
  ErrorCode rval = mb->get_entities_by_dimension(set, 2, cells);
  if (MB_SUCCESS != rval)
    return;

  // Get coefficients for reconstruction (in cubed-sphere coordinates)
  for (Range::iterator it = cells.begin(); it != cells.end(); it++)
  {
    moab::EntityHandle icell = *it;

    // get the nodes, then the coordinates const moab::EntityHandle * verts;
    const moab::EntityHandle * verts;
    int num_nodes;
    rval = mb->get_connectivity(icell, verts, num_nodes);
    if (MB_SUCCESS != rval)
      return;

    moab::Range adjacentEdges;
    rval = mb->get_adjacencies(&icell, 1, 1, true, adjacentEdges);

    // get adjacent cells from edges
    moab::Range adjacentCells;
    rval = mb->get_adjacencies(adjacentEdges, 2, true, adjacentCells, Interface::UNION);

    // get gnomonic plane 
     int plane = 0;
     rval = mb->tag_get_data(planeTag, &icell, 1, &plane );

    std::vector<double> dx(adjacentCells.size() - 1);
    std::vector<double> dy(adjacentCells.size() - 1);
    std::vector<double> dr(adjacentCells.size() - 1);
    double center_x;
    double center_y;

    // get center of cell where reconstruction occurs
     double rad = 1;
     std::vector<double> cent(3);
     rval = mb->tag_get_data(centerTag, &icell, 1, &cent[0] );
     CartVect cellcenter(cent[0],cent[1],cent[2]);
     double cellx = 0;
     double celly = 0;
     gnomonic_projection(cellcenter, rad, plane, cellx, celly);

    // get cell average value
     double cellVal = 0;
     rval = mb->tag_get_data(cellValTag, &icell, 1, &cellVal );

    // get centers of surrounding cells 
     std::vector<double> cell_cents(3*adjacentCells.size());
     rval = mb->tag_get_data(centerTag, adjacentCells, &cell_cents[0]);

    // get density of surrounding cells 
     std::vector<double> adjCellVals(adjacentCells.size());
     rval = mb->tag_get_data(cellValTag, adjacentCells, &adjCellVals[0]);

     std::size_t jind = 0;
     for (std::size_t i=0; i< adjacentCells.size(); i++){
        
         if (adjacentCells[i] != icell) {

            CartVect cent_xyz(cell_cents[i*3],cell_cents[i*3+1],cell_cents[i*3+2]);
            gnomonic_projection(cent_xyz, rad, plane, center_x, center_y);

            dx[jind] = center_x - cellx;
            dy[jind] = center_y - celly;
            dr[jind] = adjCellVals[i] - cellVal;
           
            jind++;
         }
     }

     std::vector<double> linearCoef(3);
     if (adjacentCells.size() == 5) {

       // compute normal equations matrix
        double N11 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] + dx[3]*dx[3];
        double N22 = dy[0]*dy[0] + dy[1]*dy[1] + dy[2]*dy[2] + dy[3]*dy[3];
        double N12 = dx[0]*dy[0] + dx[1]*dy[1] + dx[2]*dy[2] + dx[3]*dy[3];

       // rhs
        double Rx = dx[0]*dr[0] + dx[1]*dr[1] + dx[2]*dr[2] + dx[3]*dr[3];
        double Ry = dy[0]*dr[0] + dy[1]*dr[1] + dy[2]*dr[2] + dy[3]*dr[3];

       // determinant
        double Det = N11*N22 - N12*N12;

       // solution
        linearCoef[1] = (Rx*N22 - Ry*N12)/Det;
        linearCoef[2] = (Ry*N11 - Rx*N12)/Det;
        linearCoef[0] = cellVal - linearCoef[1]*cellx - linearCoef[2]*celly;

     }
     else
     {
        // default to first order
        linearCoef[0] = cellVal;
        linearCoef[1] = 0.0;
        linearCoef[2] = 0.0;
        std::cout<< "Need 4 adjacent cells for linear reconstruction! \n";
     }
     
     rval = mb->tag_set_data(linearCoefTag, &icell, 1, &linearCoef[0]);

  }
  return;
}

void limit_linear_reconstruction(moab::Interface * mb, moab::EntityHandle set, moab::Tag &cellValTag, 
                                 moab::Tag &planeTag, moab::Tag & boundsTag, moab::Tag & centerTag,
                                 moab::Tag &linearCoefTag)
{
  // get all entities of dimension 2
  Range cells;
  ErrorCode rval = mb->get_entities_by_dimension(set, 2, cells);
  if (MB_SUCCESS != rval)
    return;

  // Get coefficients for reconstruction (in cubed-sphere coordinates)
  for (Range::iterator it = cells.begin(); it != cells.end(); it++)
  {
    moab::EntityHandle icell = *it;

    double R = 1.0;

    // get the nodes
    const moab::EntityHandle * verts;
    int num_nodes;
    rval = mb->get_connectivity(icell, verts, num_nodes);
    if (MB_SUCCESS != rval)
      return;

    // get coordinates
    std::vector<double> coords(3 * num_nodes);
    rval = mb->get_coords(verts, num_nodes, &coords[0]);
    if (MB_SUCCESS != rval)
      return;

    // get gnomonic plane 
     int plane = 0;
     rval = mb->tag_get_data(planeTag, &icell, 1, &plane );

    // get vertex coordinates and project onto gnomonic plane
     std::vector<double> x(num_nodes);
     std::vector<double> y(num_nodes);
     for (int inode = 0; inode < num_nodes; inode++){
         double rad = sqrt(coords[inode*3]*coords[inode*3] + coords[inode*3+1]*coords[inode*3+1] + coords[inode*3+2]*coords[inode*3+2]);
         CartVect xyzcoord(coords[inode*3]/rad,coords[inode*3+1]/rad,coords[inode*3+2]/rad);
         gnomonic_projection(xyzcoord, R, plane, x[inode],y[inode]);
     }

     // get Eulerian cell reconstruction coefficients
     std::vector<double>  linearCoefs(3);
     rval = mb->tag_get_data(linearCoefTag, &icell, 1, &linearCoefs[0]);
     if (MB_SUCCESS != rval)
         return;

     // get min/max value of linear function in cell (occurs at endpoints)
     double valMin = 9999.0;
     double valMax = 0.0;
     for (int inode = 0; inode < num_nodes; inode++){
        valMin = std::min(valMin,(linearCoefs[0] + linearCoefs[1]*x[inode] + linearCoefs[2]*y[inode]));
        valMax = std::max(valMax,(linearCoefs[0] + linearCoefs[1]*x[inode] + linearCoefs[2]*y[inode]));
     }

    // get center of cell where reconstruction occurs
     double rad = 1;
     std::vector<double> cent(3);
     rval = mb->tag_get_data(centerTag, &icell, 1, &cent[0] );
     CartVect cellcenter(cent[0],cent[1],cent[2]);
     double cellx = 0;
     double celly = 0;
     gnomonic_projection(cellcenter, rad, plane, cellx, celly);

    // get cell average value
     double cellVal = 0;
     rval = mb->tag_get_data(cellValTag, &icell, 1, &cellVal );

    // get cell bounds 
     std::vector<double> bounds(2);
     rval = mb->tag_get_data(boundsTag, &icell, 1, &bounds[0] );

     // calculate limiting coefficient
      double alphaMin = 1.0;
      double alphaMax = 1.0;

      if (std::abs(valMin - cellVal) > 1.0e-10) {
        alphaMin = std::max(0.0, (bounds[0] - cellVal)/(valMin - cellVal)); 
      }

      if (std::abs(valMax - cellVal) > 1.0e-10) {
        alphaMax = std::max(0.0, (bounds[1] - cellVal)/(valMax - cellVal)); 
      }

      double alpha  = std::min(std::min(1.0, alphaMin), alphaMax);

      if (alpha<1.0)
      {
        linearCoefs[1] = alpha*linearCoefs[1];
        linearCoefs[2] = alpha*linearCoefs[2];
        linearCoefs[0] = cellVal - linearCoefs[1]*cellx - linearCoefs[2]*celly;
      }

      rval = mb->tag_set_data(linearCoefTag, &icell, 1, &linearCoefs[0]);

  }
  return;
}

void get_neighborhood_bounds(moab::Interface * mb, moab::EntityHandle set, moab::Tag &cellValTag, 
                             moab::Tag &boundsTag)
{
  // get all entities of dimension 2
  Range cells;
  ErrorCode rval = mb->get_entities_by_dimension(set, 2, cells);
  if (MB_SUCCESS != rval)
    return;

  // Get min/max value from neighborhood of cell
  for (Range::iterator it = cells.begin(); it != cells.end(); it++)
  {
    moab::EntityHandle icell = *it;

    // get the nodes
    const moab::EntityHandle * verts;
    int num_nodes;
    rval = mb->get_connectivity(icell, verts, num_nodes);
    if (MB_SUCCESS != rval)
      return;

    // get adjacent cells from nodes
     moab::Range adjacentCells;
     rval = mb->get_adjacencies(verts, num_nodes, 2, true, adjacentCells, Interface::UNION);

    // get cell value of middle cell
     double cellVal = 0;
     rval = mb->tag_get_data(cellValTag, &icell, 1, &cellVal );
     double minVal = cellVal;
     double maxVal = cellVal;

    // get values of surrounding cells and take min/max
     std::vector<double> adjCellVals(adjacentCells.size());
     rval = mb->tag_get_data(cellValTag, adjacentCells, &adjCellVals[0]);

     for (std::size_t i=0; i< adjacentCells.size(); i++){
        
         if (adjacentCells[i] != icell) {
            minVal = std::min(minVal,adjCellVals[i]);
            maxVal = std::max(maxVal,adjCellVals[i]);
         }
     }

     std::vector<double> bounds(2);
     bounds[0] = minVal;
     bounds[1] = maxVal;
     rval = mb->tag_set_data(boundsTag, &icell, 1, &bounds[0]);

  }
  return;
}


moab::ErrorCode get_departure_grid(moab::Interface * mb, moab::EntityHandle euler_set,
    moab::EntityHandle lagr_set, moab::EntityHandle & covering_set, int tStep, Range & connecVerts)
{
  ErrorCode rval = MB_SUCCESS;

  EntityHandle dum=0;
  Tag corrTag;
  mb->tag_get_handle(CORRTAGNAME, 1, MB_TYPE_HANDLE, corrTag,
                                             MB_TAG_DENSE, &dum);

  double t = tStep * T / numSteps; // numSteps is global; so is T
  double delta_t = T / numSteps; // this is global too, actually
  // double delta_t = 0.0001;
  // double t = delta_t;

  Range polys;
  rval = mb->get_entities_by_dimension(euler_set, 2, polys);

  // change coordinates of lagr mesh vertices
  for (Range::iterator vit = connecVerts.begin(); vit != connecVerts.end();
      vit++)
  {
    moab::EntityHandle oldV = *vit;
    CartVect posi;
    rval = mb->get_coords(&oldV, 1, &(posi[0]));
    // cslam utils, case 1
    CartVect newPos;
    departure_point_swirl_rot(posi, t, delta_t, newPos);
    newPos = radius * newPos; // do we need this? the radius should be 1
    moab::EntityHandle new_vert;
    rval = mb->tag_get_data(corrTag, &oldV, 1, &new_vert);
    // set the new position for the new vertex
    rval = mb->set_coords(&new_vert, 1, &(newPos[0]));
  }

  // if in parallel, we have to move some elements to another proc, and receive other cells
  // from other procs
  rval = pworker->create_departure_mesh_3rd_alg(lagr_set, covering_set);


  return rval;
}

// !!! For now serial !!!
moab::ErrorCode update_density(moab::Interface * mb, moab::EntityHandle euler_set,
    moab::EntityHandle lagr_set, moab::EntityHandle out_set, moab::Tag & rhoTag,
    moab::Tag & areaTag, moab::Tag & rhoCoefTag, moab::Tag & weightsTag,
    moab::Tag & planeTag)
{

//  moab::ParallelComm * parcomm = ParallelComm::get_pcomm(mb, 0);

  ErrorCode rval = MB_SUCCESS;

  double R = 1.0;

  moab::EntityHandle dum=0;
  Tag corrTag;
  mb->tag_get_handle(CORRTAGNAME, 1, MB_TYPE_HANDLE, corrTag,
                                             MB_TAG_DENSE, &dum);

  moab::Tag gid;
  rval = mb->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid, MB_TAG_DENSE);
    if (MB_SUCCESS != rval)
      return rval;

  // get all polygons out of out_set; then see where are they coming from
  moab::Range polys;
  rval = mb->get_entities_by_dimension(out_set, 2, polys);
    if (MB_SUCCESS != rval)
      return rval;

  // get all Lagrangian cells 
  moab::Range rs1;
  rval = mb->get_entities_by_dimension(lagr_set, 2, rs1);
    if (MB_SUCCESS != rval)
      return rval;

  // get all Eulerian cells 
  moab::Range rs2;
  rval = mb->get_entities_by_dimension(euler_set, 2, rs2);
    if (MB_SUCCESS != rval)
      return rval;

  // get gnomonic plane for Eulerian cells
  std::vector<int>  plane(rs2.size());
  rval = mb->tag_get_data(planeTag, rs2, &plane[0]);
    if (MB_SUCCESS != rval)
      return rval;

  // get Eulerian cell reconstruction coefficients
  std::vector<double>  rhoCoefs(3*rs2.size());
  rval = mb->tag_get_data(rhoCoefTag, rs2, &rhoCoefs[0]);
    if (MB_SUCCESS != rval)
      return rval;

  // get intersection weights 
  std::vector<double>  weights(6*polys.size());
  rval = mb->tag_get_data(weightsTag, polys, &weights[0]);
    if (MB_SUCCESS != rval)
      return rval;

  // Initialize the new values
  std::vector<double> newValues(rs2.size(), 0.);// initialize with 0 all of them

  // For each polygon get red/blue parent
  moab::Tag redParentTag;
  moab::Tag blueParentTag;
  rval = mb->tag_get_handle("RedParent", 1, MB_TYPE_INTEGER, redParentTag, MB_TAG_DENSE);
    if (MB_SUCCESS != rval)
      return rval;
  rval = mb->tag_get_handle("BlueParent", 1, MB_TYPE_INTEGER, blueParentTag, MB_TAG_DENSE);
    if (MB_SUCCESS != rval)
      return rval;

  // mass_lagr = (\sum_intx \int rho^h(x,y) dV)
  // rho_eul^n+1 = mass_lagr/area_eul
  double check_intx_area = 0.;
  int polyIndex = 0;
  for (Range::iterator it= polys.begin(); it!=polys.end(); it++)
  {

    moab::EntityHandle poly=*it;
    int blueIndex, redIndex;
    rval =  mb->tag_get_data(blueParentTag, &poly, 1, &blueIndex);
    moab::EntityHandle blue = rs1[blueIndex];
       
    rval = mb->tag_get_data(redParentTag, &poly, 1, &redIndex);

    moab::EntityHandle redArr;
    rval = mb->tag_get_data(corrTag, &blue, 1, &redArr);
    int arrRedIndex = rs2.index(redArr);

    // sum into new density values
    newValues[arrRedIndex] += rhoCoefs[redIndex*3]*weights[polyIndex*6] + rhoCoefs[redIndex*3+1]*weights[polyIndex*6+1] 
                                + rhoCoefs[redIndex*3+2]*weights[polyIndex*6+2];

    check_intx_area += weights[polyIndex*6];

    polyIndex++;

  }


 // now divide by red area (current)
  int j=0;
  Range::iterator iter = rs2.begin();
  void * data=NULL; //used for stored area
  int count =0;
  double total_mass_local=0.;
  while (iter != rs2.end())
  {
    rval = mb->tag_iterate(areaTag, iter, rs2.end(), count, data);
    double * ptrArea=(double*)data;
    for (int i=0; i<count; i++, iter++, j++, ptrArea++)
    {
      total_mass_local+=newValues[j];
      newValues[j]/= (*ptrArea);
    }
  }

  rval = mb->tag_set_data(rhoTag, rs2, &newValues[0]);

  std::cout <<"total mass now:" << total_mass_local << "\n";
  std::cout <<"check: total intersection area: (4 * M_PI * R^2): "  << 4 * M_PI * R*R << " " << check_intx_area << "\n";

  return MB_SUCCESS;
}

// !!! For now serial !!!
moab::ErrorCode update_tracer(moab::Interface * mb, moab::EntityHandle euler_set,
    moab::EntityHandle lagr_set, moab::EntityHandle out_set, moab::Tag & tauTag, 
    moab::Tag & areaTag, moab::Tag & rhoCoefTag, moab::Tag & tauCoefTag, moab::Tag & weightsTag, 
    moab::Tag & planeTag)
{

//  moab::ParallelComm * parcomm = ParallelComm::get_pcomm(mb, 0);

  ErrorCode rval = MB_SUCCESS;

  double R = 1.0;

  moab::EntityHandle dum=0;
  Tag corrTag;
  mb->tag_get_handle(CORRTAGNAME, 1, MB_TYPE_HANDLE, corrTag,
                                             MB_TAG_DENSE, &dum);

  moab::Tag gid;
  rval = mb->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid, MB_TAG_DENSE);
    if (MB_SUCCESS != rval)
      return rval;

  // get all polygons out of out_set; then see where are they coming from
  moab::Range polys;
  rval = mb->get_entities_by_dimension(out_set, 2, polys);
    if (MB_SUCCESS != rval)
      return rval;

  // get all Lagrangian cells 
  moab::Range rs1;
  rval = mb->get_entities_by_dimension(lagr_set, 2, rs1);
    if (MB_SUCCESS != rval)
      return rval;

  // get all Eulerian cells 
  moab::Range rs2;
  rval = mb->get_entities_by_dimension(euler_set, 2, rs2);
    if (MB_SUCCESS != rval)
      return rval;

  // get gnomonic plane for Eulerian cells
  std::vector<int>  plane(rs2.size());
  rval = mb->tag_get_data(planeTag, rs2, &plane[0]);
    if (MB_SUCCESS != rval)
      return rval;

  // get Eulerian cell reconstruction coefficients
  std::vector<double>  rhoCoefs(3*rs2.size());
  rval = mb->tag_get_data(rhoCoefTag, rs2, &rhoCoefs[0]);
    if (MB_SUCCESS != rval)
      return rval;

  // get Eulerian cell reconstruction coefficients
  std::vector<double>  tauCoefs(3*rs2.size());
  rval = mb->tag_get_data(tauCoefTag, rs2, &tauCoefs[0]);
    if (MB_SUCCESS != rval)
      return rval;

  // get intersection weights 
  std::vector<double>  weights(6*polys.size());
  rval = mb->tag_get_data(weightsTag, polys, &weights[0]);
    if (MB_SUCCESS != rval)
      return rval;

  // Initialize the new values
  std::vector<double> newValues(rs2.size(), 0.);// initialize with 0 all of them
  std::vector<double> newMass(rs2.size(), 0.);// initialize with 0 all of them

  // For each polygon get red/blue parent
  moab::Tag redParentTag;
  moab::Tag blueParentTag;
  rval = mb->tag_get_handle("RedParent", 1, MB_TYPE_INTEGER, redParentTag, MB_TAG_DENSE);
    if (MB_SUCCESS != rval)
      return rval;
  rval = mb->tag_get_handle("BlueParent", 1, MB_TYPE_INTEGER, blueParentTag, MB_TAG_DENSE);
    if (MB_SUCCESS != rval)
      return rval;

  // tracer_mass_lagr = (\sum_intx \int rho^h(x,y) tau^h(x,y) dV)
  // tau_eul^n+1 = tracer_mass_lagr/mass_lagr
  double check_intx_area = 0.;
  int polyIndex = 0;
  for (Range::iterator it= polys.begin(); it!=polys.end(); it++)
  {

    moab::EntityHandle poly=*it;
    int blueIndex, redIndex;
    rval =  mb->tag_get_data(blueParentTag, &poly, 1, &blueIndex);
    moab::EntityHandle blue = rs1[blueIndex];
       
    rval = mb->tag_get_data(redParentTag, &poly, 1, &redIndex);

    moab::EntityHandle redArr;
    rval = mb->tag_get_data(corrTag, &blue, 1, &redArr);
    int arrRedIndex = rs2.index(redArr);

    // sum into new values
    newMass[arrRedIndex] += rhoCoefs[redIndex*3]*weights[polyIndex*6] + rhoCoefs[redIndex*3+1]*weights[polyIndex*6+1] 
                                + rhoCoefs[redIndex*3+2]*weights[polyIndex*6+2];

    double A = rhoCoefs[redIndex*3]*tauCoefs[redIndex*3];
    double B = rhoCoefs[redIndex*3+1]*tauCoefs[redIndex*3] + rhoCoefs[redIndex*3]*tauCoefs[redIndex*3+1];
    double C = rhoCoefs[redIndex*3+2]*tauCoefs[redIndex*3] + rhoCoefs[redIndex*3]*tauCoefs[redIndex*3+2];
    double D = rhoCoefs[redIndex*3+1]*tauCoefs[redIndex*3+1];
    double E = rhoCoefs[redIndex*3+2]*tauCoefs[redIndex*3+2];
    double F = rhoCoefs[redIndex*3+1]*tauCoefs[redIndex*3+2] + rhoCoefs[redIndex*3+2]*tauCoefs[redIndex*3+1];
    newValues[arrRedIndex] += A*weights[polyIndex*6] + B*weights[polyIndex*6+1] + C*weights[polyIndex*6+2]
                              + D*weights[polyIndex*6+3] + E*weights[polyIndex*6+4] + F*weights[polyIndex*6+5];

    check_intx_area += weights[polyIndex*6];

    polyIndex++;

  }


 // now divide by red area (current)
  int j=0;
  Range::iterator iter = rs2.begin();
  void * data=NULL; //used for stored area
  int count =0;
  double total_mass_local=0.;
  while (iter != rs2.end())
  {
    rval = mb->tag_iterate(areaTag, iter, rs2.end(), count, data);
    double * ptrArea=(double*)data;
    for (int i=0; i<count; i++, iter++, j++, ptrArea++)
    {
      total_mass_local+=newValues[j];
      if (newMass[j] > 1.0e-12)
      {
         newValues[j]/= newMass[j];
      }
      else
      {
         newValues[j] = 0.0;
      }
    }
  }

  rval = mb->tag_set_data(tauTag, rs2, &newValues[0]);

  std::cout <<"total tracer mass now:" << total_mass_local << "\n";
  std::cout <<"check: total intersection area: (4 * M_PI * R^2): "  << 4 * M_PI * R*R << " " << check_intx_area << "\n";

  return MB_SUCCESS;
}



/*
 *  Deformational flow 
 */
void departure_point_swirl(moab::CartVect & arrival_point, double t, double delta_t, moab::CartVect & departure_point)
{

  // always assume radius is 1 here?
  moab::SphereCoords sph = cart_to_spherical(arrival_point);
  double k = 2.4; //flow parameter
  /*     radius needs to be within some range   */
  double  sl2 = sin(sph.lon/2);
  double pit = M_PI * t / T;
  double omega = M_PI/T;
  double costheta = cos(sph.lat);
  //double u = k * sl2*sl2 * sin(2*sph.lat) * cos(pit);
  double v = k * sin(sph.lon) * costheta * cos(pit);
  //double psi = k * sl2 * sl2 *costheta * costheta * cos(pit);
  double u_tilda = 2*k*sl2*sl2*sin(sph.lat)*cos(pit);

  // formula 35, page 8
  // this will approximate dep point using a Taylor series with up to second derivative
  // this will be O(delta_t^3) exact.
  double lon_dep = sph.lon - delta_t * u_tilda -delta_t*delta_t * k * sl2 *
      ( sl2 * sin (sph.lat) * sin(pit) * omega
          - u_tilda * sin(sph.lat) * cos(pit) * cos (sph.lon/2)
          - v * sl2 * costheta * cos(pit)   );
  // formula 36, page 8 again
  double lat_dep = sph.lat - delta_t*v - delta_t * delta_t/4* k *
      ( sin(sph.lon)* cos(sph.lat) * sin(pit) * omega
          - u_tilda * cos(sph.lon) * cos(sph.lat) * cos(pit)
          + v * sin(sph.lon) * sin(sph.lat) * cos(pit)  );
  moab::SphereCoords sph_dep;
  sph_dep.R = 1.; // radius
  sph_dep.lat = lat_dep;
  sph_dep.lon = lon_dep;

  departure_point = spherical_to_cart(sph_dep);
  return;
}

/*
 *  Deformational flow with rotation
 */
void departure_point_swirl_rot(moab::CartVect & arrival_point, double t, double delta_t, moab::CartVect & departure_point)
{

  moab::SphereCoords sph = cart_to_spherical(arrival_point);
  double omega = M_PI/T;
  double gt = cos(M_PI*t/T);

  double lambda = sph.lon - 2.0*omega*t;
  double u_tilda = 4.0*sin(lambda)*sin(lambda)*sin(sph.lat)*gt + 2.0*omega;
  double v = 2.0*sin(2.0*lambda)*cos(sph.lat)*gt;

  double lon_dep = sph.lon - delta_t*u_tilda - delta_t*delta_t*2.0*sin(lambda) *
                   (   sin(lambda)*sin(sph.lat)*sin(omega*t)*omega
                     - sin(lambda)*cos(sph.lat)*cos(omega*t)*v
                     - 2.0*cos(lambda)*sin(sph.lat)*cos(omega*t)*u_tilda);

  double lat_dep = sph.lat - delta_t*v - delta_t*delta_t*2.0*
                   (  cos(sph.lat)*sin(omega*t)*omega*sin(lambda)*cos(lambda)
                    - 2.0*u_tilda*cos(sph.lat)*cos(omega*t)*cos(lambda)*cos(lambda)
                    + u_tilda*cos(sph.lat)*cos(omega*t) 
                    + v*sin(sph.lat)*cos(omega*t)*sin(lambda)*cos(lambda));

  moab::SphereCoords sph_dep;
  sph_dep.R = 1.; // radius
  sph_dep.lat = lat_dep;
  sph_dep.lon = lon_dep;

  departure_point = spherical_to_cart(sph_dep);
  return;
}

/*
 *  Zonal flow
 */
void departure_point_rotation(CartVect & arrival_point, double t, double delta_t, CartVect & departure_point)
{

  // rotation angle (0 - around equator, pi/2 - over poles)
  double alpha = 0.0;

  //radius = 1
  moab::SphereCoords sph = cart_to_spherical(arrival_point);

  // angular velocity
  double omega = 2.0*M_PI;

  // lat/lon of rotated pole
  double lon_p = M_PI;
  double lat_p = M_PI/2.0 - alpha;

  // rotate spherical coordinates so that rotation is along new equator
   double coslat = cos(sph.lat);
   double sinlat = sin(sph.lat);
   double coslatp = cos(lat_p);
   double sinlatp = sin(lat_p);
   double lon_rot = atan2(coslat*sin(sph.lon - lon_p), coslat*sinlatp*cos(sph.lon-lon_p) - coslatp*sinlat);
   double lat_rot = asin(sinlat*sinlatp + coslat*coslatp*cos(sph.lon - lon_p));

  // update position in rotated coords (NOTE: for rotation in these coords only longitude changes)
   lon_rot = lon_rot - delta_t*omega;

  // convert back to standard spherical coords
   double lon_temp = atan2(cos(lat_rot)*sin(lon_rot), sin(lat_rot)*coslatp + cos(lat_rot)*cos(lon_rot)*sinlatp);
   double lon_dep = lon_temp + lon_p; 
   double lat_dep = asin(sin(lat_rot)*sinlatp - cos(lat_rot)*coslatp*cos(lon_rot));

   SphereCoords sph_dep;
   sph_dep.R = 1.; // radius
   sph_dep.lat = lat_dep;
   sph_dep.lon = lon_dep;

   departure_point = spherical_to_cart(sph_dep);

  return;
}

void get_intersection_weights(moab::Interface * mb, moab::EntityHandle euler_set, moab::EntityHandle lagr_set, 
                              moab::EntityHandle intx_set, moab::Tag &planeTag, moab::Tag &weightsTag)
{
  // get all intersection polygons
  moab::Range polys;
  ErrorCode rval = mb->get_entities_by_dimension(intx_set, 2, polys);
  if (MB_SUCCESS != rval)
    return;

  // get all Eulerian cells
  moab::Range eul_cells;
  rval = mb->get_entities_by_dimension(euler_set, 2, eul_cells);
  if (MB_SUCCESS != rval)
    return;

  // get all Lagrangian cells
  moab::Range lagr_cells;
  rval = mb->get_entities_by_dimension(lagr_set, 2, lagr_cells);
  if (MB_SUCCESS != rval)
    return;

  // get tag for Eulerian parent cell of intersection polygon
  moab::Tag redParentTag;
  rval = mb->tag_get_handle("RedParent", 1, MB_TYPE_INTEGER, redParentTag, MB_TAG_DENSE);

  // get tag for Lagrangian parent cell of intersection polygon
  moab::Tag blueParentTag;
  rval = mb->tag_get_handle("BlueParent", 1, MB_TYPE_INTEGER, blueParentTag, MB_TAG_DENSE);

  // get gnomonic plane for Eulerian cells
  std::vector<int>  plane(eul_cells.size());
  rval = mb->tag_get_data(planeTag, eul_cells, &plane[0]);
    if (MB_SUCCESS != rval)
      return;

  double total_area = 0.;
  for (moab::Range::iterator it = polys.begin(); it != polys.end(); it++)
  {
    moab::EntityHandle poly = *it;

    // get the nodes 
    const moab::EntityHandle * verts;
    int num_nodes;
    rval = mb->get_connectivity(poly, verts, num_nodes);
    if (MB_SUCCESS != rval)
      return;

    // get coordinates
    std::vector<double> coords(3 * num_nodes);
    rval = mb->get_coords(verts, num_nodes, &coords[0]);
    if (MB_SUCCESS != rval)
      return;

    // get index of Eulerian parent cell for polygon
    int redIndex;
    rval = mb->tag_get_data(redParentTag, &poly, 1, &redIndex);

    std::vector<double> x(num_nodes);
    std::vector<double> y(num_nodes);
    double poly_area  = 0;
    double poly_intx  = 0;
    double poly_inty  = 0;
    double poly_intx2 = 0;
    double poly_inty2 = 0;
    double poly_intxy = 0;
    double R = 1.0;
    for (int inode = 0; inode < num_nodes; inode++){
         double rad = sqrt(coords[inode*3]*coords[inode*3] + coords[inode*3+1]*coords[inode*3+1] + coords[inode*3+2]*coords[inode*3+2]);
         moab::CartVect xyzcoord(coords[inode*3]/rad,coords[inode*3+1]/rad,coords[inode*3+2]/rad);
         gnomonic_projection(xyzcoord, R, plane[redIndex], x[inode],y[inode]);
    }

    std::vector<double> weights(6);
    for (int inode = 0; inode < num_nodes; inode++){
        int inode2 = inode+1;
        if (inode2 >= num_nodes) inode2 = 0;
        double xmid = 0.5*(x[inode] + x[inode2]);
        double ymid = 0.5*(y[inode] + y[inode2]);
        double r1 = sqrt(1 + x[inode]*x[inode] + y[inode]*y[inode]);
        double rm = sqrt(1 + xmid*xmid + ymid*ymid);
        double r2 = sqrt(1 + x[inode2]*x[inode2] + y[inode2]*y[inode2]);
        double hx=x[inode2]-x[inode];

        poly_area  += -hx*(y[inode]/(r1*(1.0 + x[inode]*x[inode]))
                            + 4.0*ymid/(rm*(1.0 + xmid*xmid))
                            + y[inode2]/(r2*(1.0 + x[inode2]*x[inode2])))/6.0;

        poly_intx  += -hx*(x[inode]*y[inode]/(r1*(1.0 + x[inode]*x[inode]))
                            + 4.0*xmid*ymid/(rm*(1.0 +xmid*xmid))
                            + x[inode2]*y[inode2]/(r2*(1.0 + x[inode2]*x[inode2])))/6.0;

        poly_inty  += hx*(1.0/r1 + 4.0/rm + 1.0/r2)/6.0;

        poly_intx2 += -hx*(x[inode]*x[inode]*y[inode]/(r1*(1.0 + x[inode]*x[inode]))
                            + 4.0*xmid*xmid*ymid/(rm*(1+xmid*xmid))
                            + x[inode2]*x[inode2]*y[inode2]/(r2*(1 + x[inode2]*x[inode2])))/6.0;

        poly_inty2 += hx*((y[inode]/r1 - asinh(y[inode]/sqrt(1.0 + x[inode]*x[inode])))
                         + 4.0*(ymid/rm - asinh(ymid/sqrt(1.0 + xmid*xmid)))
                         + (y[inode2]/r2 - asinh(y[inode2]/sqrt(1.0 + x[inode2]*x[inode2]))))/6.0;

        poly_intxy += hx*(x[inode]/r1 + 4.0*xmid/rm + x[inode2]/r2)/6.0;
    }
     weights[0] = poly_area;
     weights[1] = poly_intx;
     weights[2] = poly_inty;
     weights[3] = poly_intx2;
     weights[4] = poly_inty2;
     weights[5] = poly_intxy;

     total_area += poly_area;

     rval = mb->tag_set_data(weightsTag, &poly, 1, &weights[0]);

  }

      std::cout << "polygon area = " << total_area << "\n";
  return;
}
