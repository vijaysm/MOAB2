
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

#ifdef MESHDIR
std::string TestDir( STRINGIFY(MESHDIR) );
#else
std::string TestDir(".");
#endif

//std::string file_name("./uniform_30.g");
//std::string file_name("./uniform_120.g");
//std::string file_name("./eulerHomme.vtk");

using namespace moab;

moab::ErrorCode update_density(moab::Interface * mb, moab::EntityHandle euler_set,
    moab::EntityHandle lagr_set, moab::EntityHandle out_set, moab::Tag & rhoTag, 
    moab::Tag & areaTag, moab::Tag & rhoCoefsTag, moab::Tag & weightsTag, moab::Tag & planeTag);

moab::ErrorCode get_departure_grid(moab::Interface * mb, moab::EntityHandle euler_set,
                                   moab::EntityHandle lagr_set, moab::EntityHandle covering_set, int tStep, moab::Range & connecVerts);

double gtol = 1.e-9; // this is for geometry tolerance

double radius = 1.;

bool writeFiles = true;
bool parallelWrite = false;
bool velocity = false;

int numSteps = 200; // number of times with velocity displayed at points
double T = 5;

Intx2MeshOnSphere * pworker = NULL;

int main(int argc, char *argv[]) {

   MPI_Init(&argc, &argv);

  // set up MOAB interface and parallel communication
   moab::Core moab;
   moab::Interface& mb = moab;
   moab::ParallelComm mb_pcomm(&mb, MPI_COMM_WORLD);

   //int rank = mb_pcomm->proc_config().proc_rank();
   int rank = mb_pcomm.proc_config().proc_rank();

  // create meshset
   moab::EntityHandle euler_set;
   moab::ErrorCode rval = mb.create_meshset(MESHSET_SET, euler_set);

   std::stringstream opts;
   //opts << "PARALLEL=READ_PART;PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;GATHER_SET=0;PARTITION_METHOD=TRIVIAL_PARTITION;VARIABLE=";
   //opts << "PARALLEL=READ_PART;PARTITION;PARALLEL_RESOLVE_SHARED_ENTS";
   //rval = mb.load_file(file_name.c_str(), &euler_set, opts.str().c_str());
   std::string fileN = TestDir + "/mbcslam/fine4.h5m";

   rval = mb.load_file(fileN.c_str(), &euler_set);
   CHECK_ERR(rval);

   // Create tag for cell density
    moab::Tag rhoTag = 0;
    rval=mb.tag_get_handle("Density",1,moab::MB_TYPE_DOUBLE,rhoTag,moab::MB_TAG_CREAT|moab::MB_TAG_DENSE);
    CHECK_ERR(rval);

   // Create tag for cell area
    moab::Tag areaTag = 0;
    rval=mb.tag_get_handle("Area",1,moab::MB_TYPE_DOUBLE,areaTag,moab::MB_TAG_CREAT|moab::MB_TAG_DENSE);
    CHECK_ERR(rval);
   // Create tag for cell barycenters in 3D Cartesian space
    moab::Tag barycenterTag = 0;
    rval=mb.tag_get_handle("CellBarycenter",3,moab::MB_TYPE_DOUBLE,barycenterTag,moab::MB_TAG_CREAT|moab::MB_TAG_DENSE);
    CHECK_ERR(rval);
   // Create tag for cell density reconstruction coefficients
    moab::Tag rhoCoefTag = 0;
    rval=mb.tag_get_handle("LinearCoefRho",3,moab::MB_TYPE_DOUBLE,rhoCoefTag,moab::MB_TAG_CREAT|moab::MB_TAG_DENSE);
    CHECK_ERR(rval);
   // Create tag for index of gnomonic plane for each cell
    moab::Tag planeTag = 0;
    rval=mb.tag_get_handle("gnomonicPlane",1,moab::MB_TYPE_INTEGER,planeTag,moab::MB_TAG_CREAT|moab::MB_TAG_DENSE);
    CHECK_ERR(rval);
   // Create tag for intersection weights
    moab::Tag weightsTag = 0;
    rval=mb.tag_get_handle("Weights",3,moab::MB_TYPE_DOUBLE,weightsTag,moab::MB_TAG_CREAT|moab::MB_TAG_DENSE);
    CHECK_ERR(rval);
   // get cell plane
    get_gnomonic_plane(&mb, euler_set, planeTag);

   // get cell barycenters (and cell area)
    get_barycenters(&mb, euler_set, planeTag, areaTag, barycenterTag);

   // Set density distributions
    set_density(&mb,euler_set, barycenterTag, rhoTag, 1);

   // Get initial values for use in error computation
    moab::Range redEls;
    rval = mb.get_entities_by_dimension(euler_set, 2, redEls);
    std::vector<double> iniValsRho(redEls.size());
    rval = mb.tag_get_data(rhoTag, redEls, &iniValsRho[0]);
    CHECK_ERR(rval);
   // Get Lagrangian set
    moab::EntityHandle out_set, lagrange_set, covering_set;  
    rval = mb.create_meshset(MESHSET_SET, out_set);
    CHECK_ERR(rval);
    rval = mb.create_meshset(MESHSET_SET, lagrange_set);
    CHECK_ERR(rval);
    rval = mb.create_meshset(MESHSET_SET, covering_set);
    CHECK_ERR(rval);
    rval = deep_copy_set(&mb, euler_set, lagrange_set);
    CHECK_ERR(rval);
    moab::EntityHandle dum = 0;
    moab::Tag corrTag;
    rval = mb.tag_get_handle(CORRTAGNAME,1, MB_TYPE_HANDLE, corrTag,
                                           MB_TAG_DENSE|MB_TAG_CREAT, &dum);
    CHECK_ERR(rval);
   //Set up intersection of two meshes

/*
    moab::Intx2MeshOnSphere worker(&mb);
    worker.SetRadius(radius);
    worker.SetErrorTolerance(gtol);
*/


    // pworker is global here; maybe we should pass it around
    pworker = new Intx2MeshOnSphere(&mb);
    pworker->SetErrorTolerance(gtol);
    pworker->SetRadius(radius);
    pworker->set_box_error(100*gtol);


    // these stay fixed for one run
    moab::Range local_verts;
    rval = pworker->build_processor_euler_boxes(euler_set, local_verts);// output also the local_verts
    //rval = worker.build_processor_euler_boxes(euler_set, local_verts);// output also the local_verts

    // loop over time to update density
    for (int ts=1; ts < numSteps + 1; ts++){
 
       if (ts  == 1)  // output initial condition
       {
           std::stringstream newDensity;
           newDensity << "Density" << rank << "_" << ts-1 << ".vtk";
           rval = mb.write_file(newDensity.str().c_str(), 0, 0, &euler_set, 1);
       }

      // get linear reconstruction coefficients
       get_linear_reconstruction(&mb, euler_set, rhoTag, planeTag, barycenterTag, rhoCoefTag);

      // get depature grid
       rval =  get_departure_grid(&mb, euler_set, lagrange_set, covering_set, ts, local_verts);

      // intersect the meshes
       rval = pworker->intersect_meshes(lagrange_set, euler_set, out_set);

      // intersection weights (i.e. area, x integral, and y integral over cell intersections)
       get_intersection_weights3(&mb, euler_set, lagrange_set, out_set, planeTag, weightsTag);

      // update the density
       rval = update_density(&mb, euler_set, lagrange_set, out_set, rhoTag, areaTag,
                              rhoCoefTag, weightsTag, planeTag);

       if (writeFiles && (ts % 5 == 0)) // so if write
       {
           std::stringstream newDensity;
           newDensity << "Density" << rank << "_" << ts << ".vtk";
           rval = mb.write_file(newDensity.str().c_str(), 0, 0, &euler_set, 1);
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
        if (rank==0)
            std::cout << " step: " << ts << "\n";

    }

    //final vals and errors
    moab::Range::iterator iter = redEls.begin();
    double norm1 = 0.;
    double norm2 = 0.;
    double exact2 = 0.;
    double exact1 = 0.;
    int count =0;
    void * data;
    int j=0;// index in iniVals
    while (iter != redEls.end())
    {
       rval = mb.tag_iterate(rhoTag, iter, redEls.end(), count, data);
       double * ptrTracer=(double*)data;

       rval = mb.tag_iterate(areaTag, iter, redEls.end(), count, data);
       double * ptrArea=(double*)data;
       for (int i=0; i<count; i++, iter++, ptrTracer++, ptrArea++, j++)
       {
          //double area = *ptrArea;
          norm1+=fabs(*ptrTracer - iniValsRho[j])* (*ptrArea);
          norm2+=(*ptrTracer - iniValsRho[j])*(*ptrTracer - iniValsRho[j])* (*ptrArea);
          exact1+=(iniValsRho[j])* (*ptrArea);
          exact2+=(iniValsRho[j])*(iniValsRho[j])* (*ptrArea);
       }
    }

   double total_norm1=0;
   double total_norm2=0;
   double total_exact1=0;
   double total_exact2=0;
   int mpi_err = MPI_Reduce(&norm1, &total_norm1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if (mpi_err)
     std::cout <<" error in MPI_reduce:" << mpi_err << "\n";
   mpi_err = MPI_Reduce(&norm2, &total_norm2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   mpi_err = MPI_Reduce(&exact1, &total_exact1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   mpi_err = MPI_Reduce(&exact2, &total_exact2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if (0==rank)
     std::cout << " numSteps:" << numSteps << " 1-norm:" << total_norm1/total_exact1 << " 2-norm:" << total_norm2/total_exact2 << "\n";

   MPI_Finalize();
   return 0;
}

moab::ErrorCode get_departure_grid(moab::Interface * mb, moab::EntityHandle euler_set,
    moab::EntityHandle lagr_set, moab::EntityHandle covering_set, int tStep, Range & connecVerts)
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

// !!! For now serial !!!
moab::ErrorCode update_density(moab::Interface * mb, moab::EntityHandle euler_set,
    moab::EntityHandle lagr_set, moab::EntityHandle out_set, moab::Tag & rhoTag,
    moab::Tag & areaTag, moab::Tag & rhoCoefsTag, moab::Tag & weightsTag,
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
  rval = mb->tag_get_data(rhoCoefsTag, rs2, &rhoCoefs[0]);
    if (MB_SUCCESS != rval)
      return rval;

  // get intersection weights 
  std::vector<double>  weights(3*polys.size());
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
    newValues[arrRedIndex] += rhoCoefs[redIndex*3]*weights[polyIndex*3] + rhoCoefs[redIndex*3+1]*weights[polyIndex*3+1] 
                                + rhoCoefs[redIndex*3+2]*weights[polyIndex*3+2];

    check_intx_area += weights[polyIndex*3+2];

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


