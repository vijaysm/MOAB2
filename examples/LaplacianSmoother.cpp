/** @example LaplacianSmoother.cpp \n
 * \brief Perform Laplacian relaxation on a mesh and its dual \n
 * <b>To run</b>: mpiexec -np <np> LaplacianSmoother [filename]\n
 *
 * Briefly, Laplacian relaxation is a technique to smooth out a mesh.  The centroid of each cell is computed from its
 * vertex positions, then vertices are placed at the average of their connected cells' centroids.
 *
 * In the parallel algorithm, an extra ghost layer of cells is exchanged.  This allows us to compute the centroids
 * for boundary cells on each processor where they appear; this eliminates the need for one round of data exchange
 * (for those centroids) between processors.  New vertex positions must be sent from owning processors to processors
 * sharing those vertices.  Convergence is measured as the maximum distance moved by any vertex.  
 * 
 * In this implementation, a fixed number of iterations is performed.  The final mesh is output to 'laplacianfinal.h5m'
 * in the current directory (H5M format must be used since the file is written in parallel).
 *
 * Usage: mpiexec -n 2 valgrind ./LaplacianSmoother -f input/surfrandomtris-64part.h5m -r 2 -p 2 -n 25
 */

#include <iostream>
#include <sstream>
#include "moab/Core.hpp"
#ifdef MOAB_HAVE_MPI
#  include "moab/ParallelComm.hpp"
#  include "MBParallelConventions.h"
#endif
#include "moab/Skinner.hpp"
#include "moab/CN.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CartVect.hpp"
#include "moab/NestedRefine.hpp"

using namespace moab;
using namespace std;

#ifdef MESH_DIR
string test_file_name = string(MESH_DIR) + string("/surfrandomtris-4part.h5m");
#else
string test_file_name = string("input/surfrandomtris-4part.h5m");
#endif

#define RC MB_CHK_ERR(rval)
#define dbgprint(MSG)                                     \
  do {                                                    \
      if (!global_rank) std::cerr << MSG << std::endl;    \
  } while(false)

ErrorCode perform_laplacian_smoothing(Core *mb, Range &verts, int dim, Tag fixed, 
                                   bool use_hc=false, bool use_acc=false, int num_its=10, 
                                   double rel_eps=1e-5, double alpha=0.0, double beta=0.5, int report_its=1);

ErrorCode hcFilter(Core* mb, moab::Range& verts, int dim, Tag fixed, 
                    std::vector<double>& verts_o, std::vector<double>& verts_n, 
                    double alpha, double beta);

ErrorCode laplacianFilter(Core* mb, moab::Range& verts, int dim, Tag fixed,
                            std::vector<double>& verts_o, std::vector<double>& verts_n);

int main(int argc, char **argv)
{
  int num_its = 10;
  int num_ref = 1;
  int num_dim = 2;
  int report_its = 1;
  int num_degree = 2;
  bool use_hc=false;
  bool use_acc=false;
  double alpha = 0.5, beta = 0.0;
  double rel_eps = 1e-5;
  const int nghostrings = 1;
  ProgOptions opts;
  ErrorCode rval;
  std::stringstream sstr;
  int global_rank;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &global_rank );

  // Decipher program options from user
  opts.addOpt<int>(std::string("niter,n"),
      std::string("Number of Laplacian smoothing iterations (default=10)"), &num_its);
  opts.addOpt<double>(std::string("eps,e"),
      std::string("Tolerance for the Laplacian smoothing error (default=1e-5)"), &rel_eps);
  opts.addOpt<double>(std::string("alpha,a"),
      std::string("Tolerance for the Laplacian smoothing error (default=0.0)"), &alpha);
  opts.addOpt<double>(std::string("beta,b"),
      std::string("Tolerance for the Laplacian smoothing error (default=0.5)"), &beta);
  opts.addOpt<int>(std::string("dim,d"),
      std::string("Topological dimension of the mesh (default=2)"), &num_dim);
  opts.addOpt<std::string>(std::string("file,f"),
      std::string("Input mesh file to smoothen (default=input/surfrandomtris-4part.h5m)"), &test_file_name);
  opts.addOpt<int>(std::string("nrefine,r"),
      std::string("Number of uniform refinements to perform and apply smoothing cycles (default=1)"), &num_ref);
  opts.addOpt<int>(std::string("ndegree,p"),
      std::string("Degree of uniform refinement (default=2)"), &num_degree);
  opts.addOpt<void>(std::string("humphrey,c"),
      std::string("Use Humphrey’s Classes algorithm to reduce shrinkage of Laplacian smoother (default=false)"), &use_hc);
  opts.addOpt<void>(std::string("aitken,d"),
      std::string("Use Aitken \\delta^2 acceleration to improve convergence of Lapalace smoothing algorithm (default=false)"), &use_acc);

  opts.parseCommandLine(argc, argv);

  // get MOAB and ParallelComm instances
  Core *mb = new Core;
  if (NULL == mb) return 1;
  int global_size = 1;

  // get the ParallelComm instance
  ParallelComm *pcomm = new ParallelComm(mb, MPI_COMM_WORLD);
  global_size = pcomm->size();

  string roptions,woptions;
  if (global_size > 1) { // if reading in parallel, need to tell it how
    sstr.str("");
    sstr << "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=" 
         << num_dim << ".0." << nghostrings
         << ";DEBUG_IO=0;DEBUG_PIO=0";
    roptions = sstr.str();
    woptions = "PARALLEL=WRITE_PART";
  }

  // read the file
  moab::EntityHandle fileset,currset;
  rval = mb->create_meshset(MESHSET_SET, fileset); RC;
  currset = fileset;
  rval = mb->load_file(test_file_name.c_str(), &fileset, roptions.c_str()); RC;

  std::vector<EntityHandle> hsets(num_ref+1,fileset);
  if (num_ref) {
    // Perform uniform refinement of the smoothed mesh
#ifdef MOAB_HAVE_MPI
    NestedRefine *uref = new NestedRefine(mb, pcomm, currset);
#else
    NestedRefine *uref = new NestedRefine(mb, 0, currset);
#endif

    std::vector<int> num_degrees(num_ref, num_degree);
    rval = uref->generate_mesh_hierarchy(num_ref, &num_degrees[0], hsets); RC;

    // Now exchange 1 layer of ghost elements, using vertices as bridge
    // (we could have done this as part of reading process, using the PARALLEL_GHOSTS read option)
    rval = uref->exchange_ghosts(hsets, nghostrings); RC;

    delete uref;
  }

  for (int iref=0; iref <= num_ref; ++iref) {

    // specify which set we are currently working on
    currset = hsets[iref];

    // make tag to specify fixed vertices, since it's input to the algorithm; use a default value of non-fixed
    // so we only need to set the fixed tag for skin vertices
    Tag fixed;
    int def_val = 0;
    rval = mb->tag_get_handle("fixed", 1, MB_TYPE_INTEGER, fixed, MB_TAG_CREAT | MB_TAG_DENSE, &def_val); RC;

    // get all vertices and faces
    Range verts, faces, skin_verts;
    rval = mb->get_entities_by_type(currset, MBVERTEX, verts); RC;
    rval = mb->get_entities_by_dimension(currset, num_dim, faces); RC;
    dbgprint ( "Found " << verts.size() << " vertices and " << faces.size() << " elements" );
    
    // get the skin vertices of those faces and mark them as fixed; we don't want to fix the vertices on a
    // part boundary, but since we exchanged a layer of ghost faces, those vertices aren't on the skin locally
    // ok to mark non-owned skin vertices too, I won't move those anyway
    // use MOAB's skinner class to find the skin
    Skinner skinner(mb);
    rval = skinner.find_skin(currset, faces, true, skin_verts); RC; // 'true' param indicates we want vertices back, not faces

    std::vector<int> fix_tag(skin_verts.size(), 1); // initialized to 1 to indicate fixed
    rval = mb->tag_set_data(fixed, skin_verts, &fix_tag[0]); RC;
    // exchange tags on owned verts for fixed points
    if (global_size > 1) {
#ifdef MOAB_HAVE_MPI
      rval = pcomm->exchange_tags(fixed, skin_verts); RC;
#endif
    }

    // now perform the Laplacian relaxation
    rval = perform_laplacian_smoothing(mb, verts, num_dim, fixed, use_hc, use_acc, num_its, rel_eps, alpha, beta, report_its); RC;
    
    // output file, using parallel write
    sstr.str("");
    sstr << "output/LaplacianSmoother_" << iref << ".h5m";
    rval = mb->write_file(sstr.str().c_str(), "H5M", woptions.c_str(), &currset, 1); RC;

    // delete fixed tag, since we created it here
    rval = mb->tag_delete(fixed); RC;
  }

  delete pcomm;
  // delete MOAB instance
  delete mb;
  
  MPI_Finalize();
  return 0;
}

ErrorCode perform_laplacian_smoothing(Core *mb, Range &verts, int dim, Tag fixed, 
                                       bool use_hc, bool use_acc, int num_its, double rel_eps, 
                                       double alpha, double beta, int report_its) 
{
  ErrorCode rval;
  int global_rank = 0, global_size = 1;
  std::vector<double> verts_acc1, verts_acc2, verts_acc3;

#ifdef MOAB_HAVE_MPI
  ParallelComm *pcomm = ParallelComm::get_pcomm(mb, 0);
  global_rank = pcomm->rank();
  global_size = pcomm->size();
#endif
  
  dbgprint ( "-- Starting smoothing cycle --" );
  // perform Laplacian relaxation:
  // 1. setup: set vertex centroids from vertex coords; filter to owned verts; get fixed tags

  // get all verts coords into tag; don't need to worry about filtering out fixed verts, 
  // we'll just be setting to their fixed coords
  std::vector<double> verts_o, verts_n;
  verts_o.resize(3*verts.size());
  verts_n.resize(3*verts.size());
  rval = mb->get_coords(verts, &verts_o[0]); RC;
  const int nbytes = sizeof(double)*verts_o.size();

  Tag errt;
  double def_val = 0.0;
  rval = mb->tag_get_handle("error", 1, MB_TYPE_DOUBLE, errt, MB_TAG_CREAT | MB_TAG_DENSE, &def_val); RC;

  if (use_acc) {
    verts_acc1.resize( verts_o.size(), 0.0 );
    verts_acc2.resize( verts_o.size(), 0.0 );
    verts_acc3.resize( verts_o.size(), 0.0 );
    memcpy( &verts_acc1[0], &verts_o[0], nbytes );
    memcpy( &verts_acc2[0], &verts_o[0], nbytes );
    memcpy( &verts_acc3[0], &verts_o[0], nbytes );
  }

  double mxdelta = 0.0, global_max = 0.0;
  // 2. for num_its iterations:
  for (int nit = 0; nit < num_its; nit++) {
    
    mxdelta = 0.0;

    // 2a. Apply Laplacian Smoothing Filter to Mesh
    if (use_hc) {
      rval = hcFilter(mb, verts, dim, fixed, verts_o, verts_n, alpha, beta); RC;
    }
    else {
      rval = laplacianFilter(mb, verts, dim, fixed, verts_o, verts_n); RC;
    }

    if (use_acc) {

      if ( nit > 10 ) {
        for(unsigned i=0; i < verts_n.size(); ++i) {

          verts_acc1[i] = verts_acc2[i];
          verts_acc2[i] = verts_acc3[i];
          verts_acc3[i] = verts_n[i];

          double num = ( verts_acc3[i] * verts_acc1[i] - verts_acc2[i] * verts_acc2[i] );
          double den = ( verts_acc3[i] - 2.0 * verts_acc2[i] + verts_acc1[i] );
          verts_n[i] = (den > 1e-8 ?  num / den : verts_n[i]);

          // double num = ( verts_acc2[i] - verts_acc1[i] );
          // double den = ( verts_acc3[i] - 2.0 * verts_acc2[i] + verts_acc1[i] );
          // verts_n[i] = ( den > 1e-8 ? verts_acc1[i] - num * num / den : verts_n[i] );
        }
      }
      else {
        memcpy( &verts_acc1[0], &verts_acc2[0], nbytes );
        memcpy( &verts_acc2[0], &verts_acc3[0], nbytes );
        memcpy( &verts_acc3[0], &verts_n[0],    nbytes );
      }
    }

    // 2b. foreach owned vertex: compute change in coordinate norm
    for (unsigned v = 0; v < verts.size(); ++v) {
      double delta = (CartVect(&verts_n[3*v])-CartVect(&verts_o[3*v])).length();
      mxdelta = std::max(delta, mxdelta);
      EntityHandle vh = verts[v];
      // dbgprint( "\t\tError " << vh << ": Max delta = " << mxdelta << "." );
      rval = mb->tag_set_data(errt, &vh, 1, &mxdelta); RC;
    }

    // global reduce for maximum delta, then report it
    global_max = mxdelta;
#ifdef MOAB_HAVE_MPI
    if (global_size > 1)
      MPI_Allreduce(&mxdelta, &global_max, 1, MPI_DOUBLE, MPI_MAX, pcomm->comm());
#endif

    if (!(nit%report_its)) {
      dbgprint( "\tIterate " << nit << ": Global Max delta = " << global_max << "." );
    }

    if (global_max < rel_eps) break;
    else {
      verts_o.swap(verts_n);
    }
  }

  // write the tag back onto vertex coordinates
  rval = mb->set_coords(verts, &verts_n[0]); RC;

  dbgprint( "-- Final iterate error = " << global_max << ".\n" );
  return MB_SUCCESS;
}

/*
  Standard Laplacian Smooth Filter

  Additional references: http://www.doc.ic.ac.uk/~gr409/thesis-MSc.pdf
*/
ErrorCode laplacianFilter(Core* mb, moab::Range& verts, int dim, Tag fixed, std::vector<double>& verts_o, std::vector<double>& verts_n)
{
  ErrorCode rval;
  int index=0;
  std::vector<int> fix_tag(verts.size());
  double *data;
  bool use_updated=true;

  memcpy( &verts_n[0], &verts_o[0], sizeof(double)*verts_o.size() );

  if (use_updated)
    data = &verts_n[0];
  else
    data = &verts_o[0];

  rval = mb->tag_get_data(fixed, verts, &fix_tag[0]); RC;

  for (Range::const_iterator vit = verts.begin(); vit != verts.end(); ++vit, index+=3)
  {
    // if !fixed
    if (fix_tag[index/3]) {
      continue;
    }

    moab::Range adjverts, adjelems;
    // Find the neighboring vertices (1-ring neighborhood)
    rval = mb->get_adjacencies(&(*vit), 1, dim, false, adjelems); RC;
    rval = mb->get_connectivity(adjelems, adjverts); RC;
    adjverts.erase(*vit);

    const int nadjs = adjverts.size();
    if (nadjs)
    {
      double delta[3] = {0.0,0.0,0.0};

      // Add the vertices and divide by the number of vertices
      for (int j=0; j<nadjs; ++j)
      {
        const int joffset = verts.index(adjverts[j])*3;
        delta[0] += data[joffset+0];
        delta[1] += data[joffset+1];
        delta[2] += data[joffset+2];
      }

      verts_n[index+0] = delta[0] / nadjs;
      verts_n[index+1] = delta[1] / nadjs;
      verts_n[index+2] = delta[2] / nadjs;
      
      // dbgprint( "-- Element = " << index/3 << " adjacencies = " << nadjs );
      // dbgprint( "\t Old_Coords = " << verts_o[index] << ", New_Coords = " << verts_n[index]  );
      // dbgprint( "\t Old_Coords = " << verts_o[index+1] << ", New_Coords = " << verts_n[index+1]  );
      // dbgprint( "\t Old_Coords = " << verts_o[index+2] << ", New_Coords = " << verts_n[index+2]  );
    }
  }

  return MB_SUCCESS;
}

/*
  HC (Humphrey’s Classes) Smooth Algorithm - Reduces Shrinkage of Laplacian Smoother

  Link: http://informatikbuero.com/downloads/Improved_Laplacian_Smoothing_of_Noisy_Surface_Meshes.pdf

  Where sv - original points
      pv - previous points,
      alpha [0..1] influences previous points pv, e.g. 0
      beta  [0..1] e.g. > 0.5
*/
ErrorCode hcFilter(Core* mb, moab::Range& verts, int dim, Tag fixed, std::vector<double>& verts_o, std::vector<double>& verts_n, double alpha, double beta)
{
  ErrorCode rval;
  std::vector<double> verts_hc(verts_o.size());
  std::vector<int> fix_tag(verts.size());

  rval = mb->tag_get_data(fixed, verts, &fix_tag[0]); RC;

  // Perform Laplacian Smooth
  rval = laplacianFilter(mb, verts, dim, fixed, verts_o, verts_n); RC;

  // Compute Differences
  unsigned index=0;
  for (index = 0; index < verts_o.size(); ++index)
  {
    verts_hc[index] = verts_n[index] - (alpha * verts_o[index] + ( 1.0 - alpha ) * verts_o[index] );
    // const int offset = index*3;
    // verts_hc[offset+0] = verts_n[offset+0] - (alpha * verts_o[offset+0] + ( 1 - alpha ) * verts_o[offset+0] );
    // verts_hc[offset+1] = verts_n[offset+1] - (alpha * verts_o[offset+1] + ( 1 - alpha ) * verts_o[offset+1] );
    // verts_hc[offset+2] = verts_n[offset+2] - (alpha * verts_o[offset+2] + ( 1 - alpha ) * verts_o[offset+2] );
  }

  index=0;
  for (Range::const_iterator vit = verts.begin(); vit != verts.end(); ++vit, index+=3)
  {
    // if !fixed
    if (fix_tag[index/3]) {
      // verts_n[index+0] = verts_o[index+0];
      // verts_n[index+1] = verts_o[index+1];
      // verts_n[index+2] = verts_o[index+2];
      continue;
    }

    moab::Range adjverts, adjelems;
    // Find the neighboring vertices (1-ring neighborhood)
    rval = mb->get_adjacencies(&(*vit), 1, dim, false, adjelems); RC;
    rval = mb->get_connectivity(adjelems, adjverts); RC;
    adjverts.erase(*vit);

    const int nadjs = adjverts.size();
    double delta[3] = {0.0,0.0,0.0};

    for (int j=0; j<nadjs; ++j)
    {
      const int joffset = verts.index(adjverts[j])*3;
      delta[0] += verts_hc[joffset+0];
      delta[1] += verts_hc[joffset+1];
      delta[2] += verts_hc[joffset+2];
    }

    verts_n[index+0] -= beta * verts_hc[index+0] + ((1.0 - beta) / nadjs) * delta[0];
    verts_n[index+1] -= beta * verts_hc[index+1] + ((1.0 - beta) / nadjs) * delta[1];
    verts_n[index+2] -= beta * verts_hc[index+2] + ((1.0 - beta) / nadjs) * delta[2];
  }

  return MB_SUCCESS;
}

