#include <iostream>
#include <sstream>
#include <iomanip>
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


#ifdef MESHDIR
std::string TestDir( STRINGIFY(MESHDIR) );
#else
std::string TestDir(".");
#endif

bool debugflag = false;
using namespace moab;

// covering set is output now (!)
moab::ErrorCode get_departure_grid(moab::Interface * mb,
    moab::EntityHandle euler_set, moab::EntityHandle lagr_set,
    moab::EntityHandle & covering_set, int tStep, moab::Range & connecVerts);

double gtol = 1.e-12; // this is for geometry tolerance

double radius = 1.;

bool writeFiles = false;

int numSteps = 200; // number of times with velocity displayed at points
double T = 5;

Intx2MeshOnSphere * pworker = NULL;

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);

  ProgOptions opts;

//  std::string firstDoubleTag;
  //  opts.addOpt<string>("double_tag_cell,d", "add double tag on cells", &firstDoubleTag);

  std::string fileIn = TestDir + "/mbcslam/fine4.h5m";

  opts.addOpt<std::string>("inpFile,i",
      "Specify the input file name string (default fine4.h5m)", &fileIn);
  opts.addOpt<int>(std::string("timeSteps,t")," number of time steps (default 200)",&numSteps);

  opts.addOpt<void>("writeFiles,w", "write result files (default false) ", &writeFiles);

  opts.parseCommandLine(argc, argv);

  // set up MOAB interface and parallel communication
  moab::Core moab;
  moab::Interface& mb = moab;
  moab::ParallelComm pcomm(&mb, MPI_COMM_WORLD);

  //int rank = mb_pcomm->proc_config().proc_rank();
  int rank = pcomm.proc_config().proc_rank();

  // create meshset
  moab::EntityHandle euler_set;
  moab::ErrorCode rval = mb.create_meshset(MESHSET_SET, euler_set);

  std::stringstream opts1;
  //opts << "PARALLEL=READ_PART;PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;GATHER_SET=0;PARTITION_METHOD=TRIVIAL_PARTITION;VARIABLE=";
  opts1 << "PARALLEL=READ_PART;PARTITION;PARALLEL_RESOLVE_SHARED_ENTS";
  //rval = mb.load_file(file_name.c_str(), &euler_set, opts.str().c_str());
  // read a homme file, partitioned in 4 so far

  const char * filename_mesh1 = fileIn.c_str();
  rval = mb.load_file(filename_mesh1, &euler_set, opts1.str().c_str()); MB_CHK_ERR(rval);

  /*int num_entities;
   rval = mb.get_number_entities_by_dimension(euler_set, 2, num_entities);MB_CHK_ERR(rval);*/

  rval = pcomm.exchange_ghost_cells(2, // int ghost_dim
      0, // int bridge_dim
      1, // int num_layers
      0, // int addl_ents
      true); MB_CHK_ERR(rval);


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
      moab::MB_TAG_CREAT | moab::MB_TAG_DENSE);  MB_CHK_ERR(rval);

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
      moab::MB_TAG_CREAT | moab::MB_TAG_DENSE);  MB_CHK_ERR(rval);

  // Create tag for cell density reconstruction coefficients
  moab::Tag rhoCoefTag = 0;
  rval = mb.tag_get_handle("LinearCoefRho", 3, moab::MB_TYPE_DOUBLE, rhoCoefTag,
      moab::MB_TAG_CREAT | moab::MB_TAG_DENSE);  MB_CHK_ERR(rval);

  // Create tag for cell tracer reconstruction coefficients
  moab::Tag tauCoefTag = 0;
  rval = mb.tag_get_handle("LinearCoefTau", 3, moab::MB_TYPE_DOUBLE, tauCoefTag,
      moab::MB_TAG_CREAT | moab::MB_TAG_DENSE);
  MB_CHK_ERR(rval);

  // Create tag for index of gnomonic plane for each cell
  moab::Tag planeTag = 0;
  rval = mb.tag_get_handle("GnomonicPlane", 1, moab::MB_TYPE_INTEGER, planeTag,
      moab::MB_TAG_CREAT | moab::MB_TAG_DENSE);
  MB_CHK_ERR(rval);

  // Create tag for density bounds
  moab::Tag rhoBoundsTag = 0;
  rval = mb.tag_get_handle("DensityBounds", 2, moab::MB_TYPE_DOUBLE,
      rhoBoundsTag, moab::MB_TAG_CREAT | moab::MB_TAG_DENSE);
  MB_CHK_ERR(rval);

  // Create tag for tracer bounds
  moab::Tag tauBoundsTag = 0;
  rval = mb.tag_get_handle("TracerBounds", 2, moab::MB_TYPE_DOUBLE,
      tauBoundsTag, moab::MB_TAG_CREAT | moab::MB_TAG_DENSE);
  MB_CHK_ERR(rval);

  // Create tag for intersection weights
  moab::Tag weightsTag = 0;
  rval = mb.tag_get_handle("Weights", 6, moab::MB_TYPE_DOUBLE, weightsTag,
      moab::MB_TAG_CREAT | moab::MB_TAG_DENSE);
  MB_CHK_ERR(rval);
  moab::Tag gid;
  rval = mb.tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid,
      MB_TAG_DENSE);
  MB_CHK_ERR(rval);

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
  rval = mb.get_entities_by_dimension(euler_set, 2, redEls);
  MB_CHK_ERR(rval);
  std::vector<double> iniValsRho(redEls.size());
  rval = mb.tag_get_data(rhoTag, redEls, &iniValsRho[0]);
  MB_CHK_ERR(rval);
  std::vector<double> iniValsTau(redEls.size());
  rval = mb.tag_get_data(tauTag, redEls, &iniValsTau[0]);
  MB_CHK_ERR(rval);

  // Get Lagrangian set
  moab::EntityHandle out_set, lagrange_set, covering_set, extended_set; // covering set created again; maybe it should be created only once
  rval = mb.create_meshset(MESHSET_SET, out_set);
  MB_CHK_ERR(rval);
  rval = mb.create_meshset(MESHSET_SET, lagrange_set);
  MB_CHK_ERR(rval);
  rval = mb.create_meshset(MESHSET_SET, extended_set);
  MB_CHK_ERR(rval);
  // add all cells, which include ghosts
  rval = mb.add_entities(extended_set, allCells);
  //rval = mb.create_meshset(MESHSET_SET, covering_set); MB_CHK_ERR(rval);

  rval = deep_copy_set(&mb, euler_set, lagrange_set);
  MB_CHK_ERR(rval);
  moab::EntityHandle dum = 0;
  moab::Tag corrTag;
  rval = mb.tag_get_handle(CORRTAGNAME, 1, MB_TYPE_HANDLE, corrTag,
      MB_TAG_DENSE | MB_TAG_CREAT, &dum);
  MB_CHK_ERR(rval);

  //Set up intersection of two meshes
  pworker = new Intx2MeshOnSphere(&mb);
  pworker->SetErrorTolerance(gtol);
  pworker->SetRadius(radius);
  pworker->set_box_error(100 * gtol);

  // these stay fixed for one run
  moab::Range local_verts;
  rval = pworker->build_processor_euler_boxes(euler_set, local_verts);
  MB_CHK_ERR(rval); // output also the local_verts

  // density, tracer tags will be used to compute the linear reconstruction using neighbor information
  // for that, we will need those values on the ghost cells too
  // exchange tags for ghosted elements
  std::vector<moab::Tag> tags;
  tags.push_back(rhoTag);
  tags.push_back(tauTag);
  tags.push_back(barycenterTag);

  rval = pcomm.exchange_tags(tags, tags, allCells);
  MB_CHK_ERR(rval);
  tags.pop_back(); // barycenter tag needs to be updated earlier in the loop

  // write the part
  std::stringstream meshFile;
  meshFile << "mesh_" << rank << ".vtk";
  rval = mb.write_file(meshFile.str().c_str(), 0, 0, &extended_set, 1);
  MB_CHK_ERR(rval);

  // loop over time to update density
  for (int ts = 1; ts < numSteps + 1; ts++) {

    if (ts == 1)  // output initial condition
    {
      std::stringstream newTracer2;
      newTracer2 << "Tracer_00" << ".h5m";
      rval = mb.write_file(newTracer2.str().c_str(), 0, "PARALLEL=WRITE_PART",
          &euler_set, 1);MB_CHK_ERR(rval);
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
    rval = pcomm.exchange_tags(centerOfMassTag, allCells);
    MB_CHK_ERR(rval);
    // get linear reconstruction coefficients for tracer
    get_linear_reconstruction(&mb, euler_set, tauTag, planeTag, centerOfMassTag,
        tauCoefTag);

    // limit linear reconstruction coefficients for tracer
    limit_linear_reconstruction(&mb, euler_set, tauTag, planeTag, tauBoundsTag,
        centerOfMassTag, tauCoefTag);

    // get depature grid
    rval = get_departure_grid(&mb, euler_set, lagrange_set, covering_set, ts,
        local_verts);
    MB_CHK_ERR(rval);

    // intersect the meshes
    rval = pworker->intersect_meshes(covering_set, euler_set, out_set);
    MB_CHK_ERR(rval);

    if (writeFiles && rank < 10 && ts<2) // so if write
    {
      std::stringstream lagrFile;
      lagrFile << "lagr" << rank << "_" << ts << ".vtk";
      rval = mb.write_file(lagrFile.str().c_str(), 0, 0, &lagrange_set, 1); MB_CHK_ERR(rval);
      std::stringstream coverFile;
      coverFile << "cover" << rank << "_" << ts << ".vtk";
      rval = mb.write_file(coverFile.str().c_str(), 0, 0, &covering_set, 1);
      MB_CHK_ERR(rval);
      std::stringstream intxFile;
      intxFile << "intx" << rank << "_" << ts << ".vtk";
      rval = mb.write_file(intxFile.str().c_str(), 0, 0, &out_set, 1);
      MB_CHK_ERR(rval);
    }

    // intersection weights (i.e. area, x integral, and y integral over cell intersections)
    get_intersection_weights(&mb, euler_set, out_set, planeTag, weightsTag);

    // update the density and tracer
    rval = pworker->update_density_and_tracers(rhoTag, areaTag, rhoCoefTag,
        tauTag, tauCoefTag, weightsTag, planeTag);
    MB_CHK_ERR(rval);
    /*rval = update_density(&mb, euler_set, lagrange_set, out_set, rhoTag,
     areaTag, rhoCoefTag, weightsTag, planeTag);

     // update the tracer
     rval = update_tracer(&mb, euler_set, lagrange_set, out_set, tauTag, areaTag,
     rhoCoefTag, tauCoefTag, weightsTag, planeTag);*/

    rval = pcomm.exchange_tags(tags, tags, allCells);
    MB_CHK_ERR(rval);
    if (writeFiles) // so if write
    {
      std::stringstream newTracer2;
      newTracer2 << "Tracer_0" << ts << ".h5m";
      rval = mb.write_file(newTracer2.str().c_str(), 0, "PARALLEL=WRITE_PART",
          &euler_set, 1);
      MB_CHK_ERR(rval);
    }

    // debug
    if (debugflag) {
      std::vector<double> rhov(redEls.size());
      std::vector<int> gids(redEls.size());
      std::vector<double> tauv(redEls.size());
      mb.tag_get_data(rhoTag, redEls, &rhov[0]);
      mb.tag_get_data(gid, redEls, &gids[0]);
      mb.tag_get_data(tauTag, redEls, &tauv[0]);
      for (int p = 0; p < (int) pcomm.size(); p++) {
        if (rank == p) {
          for (size_t i = 0; i < redEls.size(); i++)
            std::cout << " gid:" << gids[i] << std::setprecision(14) << " dens:"
                << rhov[i] << " tau:" << tauv[i] << "\n";
        }
        MPI_Barrier(pcomm.comm());
      }

    }
    // delete the polygons and elements of out_set
    moab::Range allVerts;
    rval = mb.get_entities_by_dimension(0, 0, allVerts);
    MB_CHK_ERR(rval);

    moab::Range allElems;
    rval = mb.get_entities_by_dimension(0, 2, allElems);
    MB_CHK_ERR(rval);

    // get Eulerian and lagrangian cells
    moab::Range polys = allCells; // these are ghosted cells, too, in euler set, we need them
    rval = mb.get_entities_by_dimension(lagrange_set, 2, polys);
    MB_CHK_ERR(rval); // do not delete lagr set either, with its vertices

    // add to the connecVerts range all verts, from all initial polys
    moab::Range vertsToStay;
    rval = mb.get_connectivity(polys, vertsToStay);
    MB_CHK_ERR(rval);

    moab::Range todeleteVerts = subtract(allVerts, vertsToStay);

    moab::Range todeleteElem = subtract(allElems, polys);
    // empty the out mesh set
    rval = mb.clear_meshset(&out_set, 1);
    MB_CHK_ERR(rval);

    rval = mb.delete_entities(todeleteElem);
    MB_CHK_ERR(rval);
    rval = mb.delete_entities(todeleteVerts);
    MB_CHK_ERR(rval);
    if (rank == 0)
      std::cout << " step: " << ts << "\n";
    // temporary, stop here
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
    MB_CHK_ERR(rval);
    double * ptrDensity = (double*) data;

    rval = mb.tag_iterate(tauTag, iter, redEls.end(), count, data);
    MB_CHK_ERR(rval);
    double * ptrTracer = (double*) data;

    rval = mb.tag_iterate(areaTag, iter, redEls.end(), count, data);
    MB_CHK_ERR(rval);
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
  MPI_Reduce(&norm1rho, &total_norm1rho, 1, MPI_DOUBLE, MPI_SUM, 0,
      MPI_COMM_WORLD);
  MPI_Reduce(&norm2rho, &total_norm2rho, 1, MPI_DOUBLE, MPI_SUM, 0,
      MPI_COMM_WORLD);
  MPI_Reduce(&exact1rho, &total_exact1rho, 1, MPI_DOUBLE, MPI_SUM, 0,
      MPI_COMM_WORLD);
  MPI_Reduce(&exact2rho, &total_exact2rho, 1, MPI_DOUBLE, MPI_SUM, 0,
      MPI_COMM_WORLD);
  MPI_Reduce(&norm1tau, &total_norm1tau, 1, MPI_DOUBLE, MPI_SUM, 0,
      MPI_COMM_WORLD);
  MPI_Reduce(&norm2tau, &total_norm2tau, 1, MPI_DOUBLE, MPI_SUM, 0,
      MPI_COMM_WORLD);
  MPI_Reduce(&exact1tau, &total_exact1tau, 1, MPI_DOUBLE, MPI_SUM, 0,
      MPI_COMM_WORLD);
  MPI_Reduce(&exact2tau, &total_exact2tau, 1, MPI_DOUBLE, MPI_SUM, 0,
      MPI_COMM_WORLD);

  if (0 == rank) {
    std::cout << " numSteps: " << numSteps << " 1-norm rho:"
        << total_norm1rho / total_exact1rho << " 2-norm rho:"
        << std::sqrt(total_norm2rho / total_exact2rho) << "\n";
    std::cout << "                1-norm tau:"
        << total_norm1tau / total_exact1tau << " 2-norm tau:"
        << std::sqrt(total_norm2tau / total_exact2tau) << "\n";
  }
  MPI_Finalize();
  return 0;
}

moab::ErrorCode get_departure_grid(moab::Interface * mb,
    moab::EntityHandle euler_set, moab::EntityHandle lagr_set,
    moab::EntityHandle & covering_set, int tStep, Range & connecVerts) {
  ErrorCode rval = MB_SUCCESS;

  EntityHandle dum = 0;
  Tag corrTag;
  mb->tag_get_handle(CORRTAGNAME, 1, MB_TYPE_HANDLE, corrTag, MB_TAG_DENSE,
      &dum);

  double t = tStep * T / numSteps; // numSteps is global; so is T
  double delta_t = T / numSteps; // this is global too, actually
  // double delta_t = 0.0001;
  // double t = delta_t;

  Range polys;
  rval = mb->get_entities_by_dimension(euler_set, 2, polys);

  // change coordinates of lagr mesh vertices
  for (Range::iterator vit = connecVerts.begin(); vit != connecVerts.end();
      vit++) {
    moab::EntityHandle oldV = *vit;
    CartVect posi;
    rval = mb->get_coords(&oldV, 1, &(posi[0]));
    // cslam utils, case 1
    CartVect newPos;
    departure_point_swirl_rot(posi, t, T, delta_t, newPos);
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
