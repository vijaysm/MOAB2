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
#include "moab/VerdictWrapper.hpp"
#include "matrix.h"

using namespace moab;
using namespace std;

#define WRITE_DEBUG_FILES

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

ErrorCode perform_laplacian_smoothing(Core *mb, Range& cells, Range &verts, int dim, Tag fixed, 
                                   bool use_hc=false, bool use_acc=false, int acc_method=1, int num_its=10, 
                                   double rel_eps=1e-5, double alpha=0.0, double beta=0.5, int report_its=1);

ErrorCode hcFilter(Core* mb, moab::ParallelComm* pcomm, moab::Range& verts, int dim, Tag fixed, 
                    std::vector<double>& verts_o, std::vector<double>& verts_n, 
                    double alpha, double beta);

ErrorCode laplacianFilter(Core* mb, moab::ParallelComm* pcomm, moab::Range& verts, int dim, Tag fixed,
                            std::vector<double>& verts_o, std::vector<double>& verts_n, bool use_updated=true);

int main(int argc, char **argv)
{
  int num_its = 10;
  int num_ref = 0;
  int num_dim = 2;
  int report_its = 1;
  int num_degree = 2;
  bool use_hc=false;
  bool use_acc=false;
  int acc_method=1;
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
  opts.addOpt<double>(std::string("alpha"),
      std::string("Tolerance for the Laplacian smoothing error (default=0.0)"), &alpha);
  opts.addOpt<double>(std::string("beta"),
      std::string("Tolerance for the Laplacian smoothing error (default=0.5)"), &beta);
  opts.addOpt<int>(std::string("dim,d"),
      std::string("Topological dimension of the mesh (default=2)"), &num_dim);
  opts.addOpt<std::string>(std::string("file,f"),
      std::string("Input mesh file to smoothen (default=surfrandomtris-4part.h5m)"), &test_file_name);
  opts.addOpt<int>(std::string("nrefine,r"),
      std::string("Number of uniform refinements to perform and apply smoothing cycles (default=1)"), &num_ref);
  opts.addOpt<int>(std::string("ndegree,p"),
      std::string("Degree of uniform refinement (default=2)"), &num_degree);
  opts.addOpt<void>(std::string("humphrey,c"),
      std::string("Use Humphrey’s Classes algorithm to reduce shrinkage of Laplacian smoother (default=false)"), &use_hc);
  opts.addOpt<void>(std::string("aitken,d"),
      std::string("Use Aitken \\delta^2 acceleration to improve convergence of Lapalace smoothing algorithm (default=false)"), &use_acc);
  opts.addOpt<int>(std::string("acc,a"),
      std::string("Type of vector Aitken process to use for acceleration (default=1)"), &acc_method);

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

    // get all vertices and cells
    Range verts, cells, skin_verts;
    rval = mb->get_entities_by_type(currset, MBVERTEX, verts); RC;
    rval = mb->get_entities_by_dimension(currset, num_dim, cells); RC;
    dbgprint ( "Found " << verts.size() << " vertices and " << cells.size() << " elements" );
    
    // get the skin vertices of those cells and mark them as fixed; we don't want to fix the vertices on a
    // part boundary, but since we exchanged a layer of ghost cells, those vertices aren't on the skin locally
    // ok to mark non-owned skin vertices too, I won't move those anyway
    // use MOAB's skinner class to find the skin
    Skinner skinner(mb);
    rval = skinner.find_skin(currset, cells, true, skin_verts); RC; // 'true' param indicates we want vertices back, not cells

    std::vector<int> fix_tag(skin_verts.size(), 1); // initialized to 1 to indicate fixed
    rval = mb->tag_set_data(fixed, skin_verts, &fix_tag[0]); RC;
    // exchange tags on owned verts for fixed points
    if (global_size > 1) {
#ifdef MOAB_HAVE_MPI
      rval = pcomm->exchange_tags(fixed, skin_verts); RC;
#endif
    }

    // now perform the Laplacian relaxation
    rval = perform_laplacian_smoothing(mb, cells, verts, num_dim, fixed, use_hc, use_acc, acc_method, num_its, rel_eps, alpha, beta, report_its); RC;
    
    // output file, using parallel write
    sstr.str("");
    sstr << "LaplacianSmoother_" << iref << ".h5m";
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

ErrorCode perform_laplacian_smoothing(Core *mb, Range& cells, Range &verts, int dim, Tag fixed, 
                                       bool use_hc, bool use_acc, int acc_method,
                                       int num_its, double rel_eps, 
                                       double alpha, double beta, int report_its) 
{
  ErrorCode rval;
  int global_rank = 0, global_size = 1;
  int nacc = 2; /* nacc_method: 1 = Method 2 from [1], 2 = Method 3 from [1] */
  std::vector<double> verts_acc1, verts_acc2, verts_acc3;
  double rat_theta=rel_eps, rat_alpha=rel_eps, rat_alphaprev=rel_eps;
#ifdef MOAB_HAVE_MPI
  const char *woptions = "PARALLEL=WRITE_PART";
#else
  const char *woptions = "";
#endif
  std::vector<int> fix_tag(verts.size());

  rval = mb->tag_get_data(fixed, verts, &fix_tag[0]); RC;

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
  verts_o.resize(3*verts.size(), 0.0);
  verts_n.resize(3*verts.size(), 0.0);
  // std::vector<const void*> vdata(1);
  // vdata[0] = &verts_n[0];
  // const int vdatasize = verts_n.size();
  void* vdata = &verts_n[0];
  rval = mb->get_coords(verts, &verts_o[0]); RC;
  const int nbytes = sizeof(double)*verts_o.size();

  Tag errt,vpost;
  double def_val[3] = {0.0,0.0,0.0};
  rval = mb->tag_get_handle("error", 1, MB_TYPE_DOUBLE, errt, MB_TAG_CREAT | MB_TAG_DENSE, def_val); RC;
  rval = mb->tag_get_handle("vpos", 3, MB_TYPE_DOUBLE, vpost, MB_TAG_CREAT | MB_TAG_DENSE, def_val); RC;
  // rval = mb->tag_set_by_ptr(vpost, verts, vdata); RC;

  if (use_acc) {
    verts_acc1.resize( verts_o.size(), 0.0 );
    verts_acc2.resize( verts_o.size(), 0.0 );
    verts_acc3.resize( verts_o.size(), 0.0 );
    memcpy( &verts_acc1[0], &verts_o[0], nbytes );
    memcpy( &verts_acc2[0], &verts_o[0], nbytes );
    memcpy( &verts_acc3[0], &verts_o[0], nbytes );
  }

  // Filter verts down to owned ones and get fixed tag for them
  Range owned_verts, shared_owned_verts;
  if (global_size > 1) {
#ifdef MOAB_HAVE_MPI
    rval = pcomm->filter_pstatus(verts, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &owned_verts);MB_CHK_ERR(rval);
#endif
  }
  else
    owned_verts = verts;

#ifdef MOAB_HAVE_MPI
  // Get shared owned verts, for exchanging tags
  rval = pcomm->get_shared_entities(-1, shared_owned_verts, 0, false, true);MB_CHK_ERR(rval);
  // Workaround: if no shared owned verts, put a non-shared one in the list, to prevent exchanging tags
  // for all shared entities
  if (shared_owned_verts.empty()) shared_owned_verts.insert(*verts.begin());
#endif

#ifdef WRITE_DEBUG_FILES
  {
    // output file, using parallel write
    std::stringstream sstr;
    sstr << "LaplacianSmootherIterate_0.h5m";
    rval = mb->write_file(sstr.str().c_str(), NULL, woptions); RC;
  }
#endif

  bool check_metrics = true;
  VerdictWrapper vw(mb);
  vw.set_size(1.); // for relative size measures; maybe it should be user input

  double mxdelta = 0.0, global_max = 0.0;
  // 2. for num_its iterations:
  for (int nit = 0; nit < num_its; nit++) {
    
    mxdelta = 0.0;

    // 2a. Apply Laplacian Smoothing Filter to Mesh
    if (use_hc) {
      rval = hcFilter(mb, pcomm, verts, dim, fixed, verts_o, verts_n, alpha, beta); RC;
    }
    else {
      rval = laplacianFilter(mb, pcomm, verts, dim, fixed, verts_o, verts_n); RC;
    }

    if (check_metrics) {
      bool smooth_success = true;
      std::vector<double> coords;
      double jac = 0.0;
      int num_nodes = 0;
      for (unsigned ic=0; ic < cells.size(); ++ic) {
        EntityHandle ecell = cells[ic];
        EntityType etype = mb->type_from_handle(ecell);
        Range everts;
        rval = mb->get_connectivity(&ecell, 1, everts, num_nodes); RC;
        coords.resize(num_nodes*3, 0.0);
        for (int iv=0; iv < num_nodes; ++iv) {
          const int offset = mb->id_from_handle(everts[iv])*3;
          for (int ii=0; ii < 3; ++ii)
            coords[iv*3+ii] = verts_n[offset+ii];
        }
        rval = vw.quality_measure(ecell, MB_SHAPE, jac, num_nodes, etype, &coords[0]); RC;
        if (jac < 1e-8) {
          dbgprint ( "Inverted element " << ic << " with jacobian = " << jac << "." ) ;
          smooth_success = false;
          break;
        }
      }

      if (!smooth_success) {
        verts_o.swap(verts_n);
        dbgprint ( "Found element inversions during smoothing. Stopping iterations." ) ;
        break;
      }
    }

    if (use_acc) {

      // if (acc_method < 5 || nit <= 3) {
      //   memcpy( &verts_acc1[0], &verts_acc2[0], nbytes );
      //   memcpy( &verts_acc2[0], &verts_acc3[0], nbytes );
      //   memcpy( &verts_acc3[0], &verts_n[0],    nbytes );
      // }

      rat_alphaprev = rat_alpha;
      for(unsigned i=0; i < verts_n.size(); ++i) {
        rat_alpha = std::max( rat_alpha, std::abs( (verts_acc3[i] - verts_acc2[i]) * (verts_acc2[i] - verts_acc1[i]) )/ ( (verts_acc2[i] - verts_acc1[i]) * (verts_acc2[i] - verts_acc1[i]) ) );
      }
      rat_theta = std::abs( rat_alpha / rat_alphaprev - 1.0 );

      if ( nit > 3 && (nit % nacc) && rat_theta < 1.0 ) {

        if (acc_method == 1) { /* Method 2 from ACCELERATION OF VECTOR SEQUENCES: http://onlinelibrary.wiley.com/doi/10.1002/cnm.1630020409/pdf */
          double vnorm = 0.0, den, acc_alpha = 0.0, acc_gamma = 0.0;
          for(unsigned i=0; i < verts_n.size(); ++i) {
            den = ( verts_acc3[i] - 2.0 * verts_acc2[i] + verts_acc1[i] );
            vnorm += den*den;
            acc_alpha += (verts_acc3[i] - verts_acc2[i]) * (verts_acc3[i] - verts_acc2[i]);
            acc_gamma += (verts_acc2[i] - verts_acc1[i]) * (verts_acc2[i] - verts_acc1[i]);
          }
          for(unsigned i=0; i < verts_n.size(); ++i) {
            verts_n[i] = verts_acc2[i] + ( acc_gamma * (verts_acc3[i] - verts_acc2[i]) - acc_alpha * (verts_acc2[i] - verts_acc1[i]) ) / vnorm;
          }
        }
        else if (acc_method == 2) { /* Method 3 from ACCELERATION OF VECTOR SEQUENCES: http://onlinelibrary.wiley.com/doi/10.1002/cnm.1630020409/pdf */
          double vnorm = 0.0, num = 0.0, den = 0.0, mu = 0.0;
          for(unsigned i=0; i < verts_n.size(); ++i) {
            num += ( verts_acc3[i] - verts_acc2[i] ) * ( verts_acc3[i] - 2.0 * verts_acc2[i] + verts_acc1[i] );
            den = ( verts_acc3[i] - 2.0 * verts_acc2[i] + verts_acc1[i] );
            vnorm += den*den;
          }
          mu = num / vnorm;
          for(unsigned i=0; i < verts_n.size(); ++i) {
            verts_n[i] = verts_acc3[i] + mu * ( verts_acc2[i] - verts_acc3[i] );
          }
        }
        else if (acc_method == 3) { /* Method 5 from ACCELERATION OF VECTOR SEQUENCES: http://onlinelibrary.wiley.com/doi/10.1002/cnm.1630020409/pdf */
          double num = 0.0, den = 0.0, mu = 0.0;
          for(unsigned i=0; i < verts_n.size(); ++i) {
            num += ( verts_acc3[i] - verts_acc2[i] ) * ( verts_acc2[i] - verts_acc1[i] );
            den += ( verts_acc2[i] - verts_acc1[i] ) * ( verts_acc3[i] - 2.0 * verts_acc2[i] + verts_acc1[i] );
          }
          mu = num / den;
          for(unsigned i=0; i < verts_n.size(); ++i) {
            verts_n[i] = verts_acc3[i] - mu * ( verts_acc3[i] - verts_acc2[i] );
          }
        }
        else if (acc_method == 4) { /* Method 8 from ACCELERATION OF VECTOR SEQUENCES: http://onlinelibrary.wiley.com/doi/10.1002/cnm.1630020409/pdf */
          double num = 0.0, den = 0.0, lambda = 0.0;
          for(unsigned i=0; i < verts_n.size(); ++i) {
            num += ( verts_acc3[i] - verts_acc2[i] ) * ( verts_acc3[i] - verts_acc2[i] );
            den += ( verts_acc2[i] - verts_acc1[i] ) * ( verts_acc2[i] - verts_acc1[i] );
          }
          lambda = std::sqrt(num / den);
          for(unsigned i=0; i < verts_n.size(); ++i) {
            verts_n[i] = verts_acc3[i] - lambda / (lambda - 1.0) * ( verts_acc3[i] - verts_acc2[i] );
          }
        }
        else if (acc_method == 5) { /* Minimum polynomial extrapolation of vector sequenes : https://en.wikipedia.org/wiki/Minimum_polynomial_extrapolation */
          /* Pseudo-code 
                U=x(:,2:end-1)-x(:,1:end-2);
                c=-pinv(U)*(x(:,end)-x(:,end-1));
                c(end+1,1)=1;
                s=(x(:,2:end)*c)/sum(c);
          */
          Matrix U(verts_n.size(), 2);
          Vector res(verts_n.size());
          for (unsigned ir=0; ir < verts_n.size(); ir++) {
            U(ir,0) = verts_acc2[ir]-verts_acc1[ir];
            U(ir,1) = verts_acc3[ir]-verts_acc2[ir];
            res[ir] = -(verts_n[ir]-verts_acc3[ir]);
          }
          // U.print();
          // Vector acc = QR(U).solve(res);
          Vector acc = solve(U,res);
          double accsum = acc[0]+acc[1]+1.0;
          for (unsigned ir=0; ir < verts_n.size(); ir++) {
            verts_n[ir] = (verts_acc1[ir]*acc[0] + verts_acc2[ir]*acc[1] + verts_acc3[ir]) / accsum;
          }

          memcpy( &verts_acc1[0], &verts_acc2[0], nbytes );
          memcpy( &verts_acc2[0], &verts_acc3[0], nbytes );
          memcpy( &verts_acc3[0], &verts_n[0],    nbytes );

        }
        else {
          int offset=0;
          for (Range::const_iterator vit = verts.begin(); vit != verts.end(); ++vit, offset+=3)
          {
            // if !fixed
            if (fix_tag[offset/3]) continue;

            CartVect num1 = ( CartVect(&verts_acc3[offset]) - CartVect(&verts_acc2[offset]) );
            CartVect num2 = ( CartVect(&verts_acc3[offset]) - 2.0 * CartVect(&verts_acc2[offset]) + CartVect(&verts_acc1[offset]) );

            num1.scale(num2);
            const double mu = num1.length_squared() / num2.length_squared();

            verts_n[offset+0] = verts_acc3[offset+0] + mu * (verts_acc2[offset+0] - verts_acc3[offset+0]);
            verts_n[offset+1] = verts_acc3[offset+1] + mu * (verts_acc2[offset+1] - verts_acc3[offset+1]);
            verts_n[offset+2] = verts_acc3[offset+2] + mu * (verts_acc2[offset+2] - verts_acc3[offset+2]);
          }
        }
      }
      memcpy( &verts_acc1[0], &verts_acc2[0], nbytes );
      memcpy( &verts_acc2[0], &verts_acc3[0], nbytes );
      memcpy( &verts_acc3[0], &verts_n[0],    nbytes );
    }



    // 2b. foreach owned vertex: compute change in coordinate norm
    for (unsigned v = 0; v < verts.size(); ++v) {
      double delta = (CartVect(&verts_n[3*v])-CartVect(&verts_o[3*v])).length();
      mxdelta = std::max(delta, mxdelta);
      EntityHandle vh = verts[v];
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

#ifdef WRITE_DEBUG_FILES
    {
      // write the tag back onto vertex coordinates
      rval = mb->set_coords(verts, &verts_n[0]); RC;
      // output VTK file for debugging purposes
      std::stringstream sstr;
      sstr << "LaplacianSmootherIterate_" << nit+1 << ".h5m";
      rval = mb->write_file(sstr.str().c_str(), NULL, woptions); RC;
    }
#endif

#ifdef MOAB_HAVE_MPI
      // 2c. exchange tags on owned verts
    if (global_size > 1) {
      rval = mb->tag_set_data(vpost, owned_verts, vdata); RC;
      rval = pcomm->exchange_tags(vpost, shared_owned_verts); RC;
    }
#endif

    if (global_max < rel_eps) break;
    else {
      std::copy(verts_n.begin(), verts_n.end(), verts_o.begin());
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
ErrorCode laplacianFilter(Core* mb, moab::ParallelComm* pcomm, moab::Range& verts, int dim, Tag fixed, std::vector<double>& verts_o, std::vector<double>& verts_n, bool use_updated)
{
  ErrorCode rval;
  std::vector<int> fix_tag(verts.size());
  double *data;

  memcpy( &verts_n[0], &verts_o[0], sizeof(double)*verts_o.size() );

  if (use_updated)
    data = &verts_n[0];
  else
    data = &verts_o[0];

    // filter verts down to owned ones and get fixed tag for them
  Range owned_verts;
  if (pcomm->size() > 1) {
#ifdef MOAB_HAVE_MPI
    rval = pcomm->filter_pstatus(verts, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &owned_verts);
    if (rval != MB_SUCCESS) return rval;
#endif
  }
  else
    owned_verts = verts;

  rval = mb->tag_get_data(fixed, owned_verts, &fix_tag[0]); RC;

  int vindex=0;
  for (Range::const_iterator vit = owned_verts.begin(); vit != owned_verts.end(); ++vit, vindex++)
  {
    // if !fixed
    if (fix_tag[vindex]) continue;

    const int index = verts.index(*vit) * 3;

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
ErrorCode hcFilter(Core* mb, moab::ParallelComm* pcomm, moab::Range& verts, int dim, Tag fixed, std::vector<double>& verts_o, std::vector<double>& verts_n, double alpha, double beta)
{
  ErrorCode rval;
  std::vector<double> verts_hc(verts_o.size());
  std::vector<int> fix_tag(verts.size());

  // Perform Laplacian Smooth
  rval = laplacianFilter(mb, pcomm, verts, dim, fixed, verts_o, verts_n); RC;

  // Compute Differences
  for (unsigned index = 0; index < verts_o.size(); ++index)
  {
    verts_hc[index] = verts_n[index] - (alpha * verts_n[index] + ( 1.0 - alpha ) * verts_o[index] );
  }

    // filter verts down to owned ones and get fixed tag for them
  Range owned_verts;
#ifdef MOAB_HAVE_MPI
  if (pcomm->size() > 1) {
    rval = pcomm->filter_pstatus(verts, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &owned_verts);
    if (rval != MB_SUCCESS) return rval;
  }
  else
#endif
    owned_verts = verts;

  rval = mb->tag_get_data(fixed, owned_verts, &fix_tag[0]); RC;

  int vindex=0;
  for (Range::const_iterator vit = owned_verts.begin(); vit != owned_verts.end(); ++vit, vindex++)
  {
    // if !fixed
    if (fix_tag[vindex]) continue;

    const int index = verts.index(*vit) * 3;

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

