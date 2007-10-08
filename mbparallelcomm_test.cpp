/** test of MBParallelComm functionality
 *
 * To run:
 *
 * mpirun -np <#procs> mbparallelcomm_test
 *
 */

#include "MBParallelComm.hpp"
#include "MBParallelConventions.h"
#include "MBTagConventions.hpp"
#include "MBCore.hpp"
#include "ScdVertexSeq.hpp"
#include "ScdElementSeq.hpp"
#include "EntitySequenceManager.hpp"
#include "mpi.h"
#include <iostream>

#define REALTFI 1

#define ERROR(a, b) {std::cerr << a << std::endl; return b;}

#define PRINT_LAST_ERROR {\
        std::string last_error;\
        result = mbImpl->get_last_error(last_error);\
        if (last_error.empty()) std::cerr << "(none)" << std::endl;\
        else std::cerr << last_error << std::endl;\
        }
MBErrorCode create_linear_mesh(MBInterface *mbImpl,
                               int N, int M, int &nshared);

MBErrorCode create_scd_mesh(MBInterface *mbImpl,
                            int IJK, int &nshared);

MBErrorCode read_file(MBInterface *mbImpl, const char *filename);

int main(int argc, char **argv) 
{
    // need to init MPI first, to tell how many procs and rank
  int err = MPI_Init(&argc, &argv);

  int nprocs, rank;
  err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // create MOAB instance based on that
  MBInterface *mbImpl = new MBCore(rank, nprocs);
  if (NULL == mbImpl) return 1;
  
  MBErrorCode result;

    // each interior proc has a vector of N+M vertices, sharing
    // M vertices each with lower- and upper-rank processors, except
    // procs on the end

    // get N, M from command line
  int N, M;
  if (argc < 3) {
    std::cerr << "No arguments passed; assuming N=10, M=2." << std::endl;
    N = 10;
    M = 2;
  }
  else {
    N = atoi(argv[1]);
    M = atoi(argv[2]);
  }

  int max_iter = 2;
  if (argc > 3) max_iter = atoi(argv[3]);

  for (int i = 0; i < max_iter; i++) {
    int nshared;
    MBErrorCode tmp_result = MB_SUCCESS;
    if (0 == i) {
      tmp_result = create_linear_mesh(mbImpl, N, M, nshared);
      if (MB_SUCCESS != tmp_result) {
        result = tmp_result;
        std::cerr << "Couldn't create linear mesh; error message:." 
                  << std::endl;
        PRINT_LAST_ERROR
        continue;
      }
    }
    else if (1 == i) {
      tmp_result = create_scd_mesh(mbImpl, N, nshared);
      if (MB_SUCCESS != tmp_result) {
        result = tmp_result;
        std::cerr << "Couldn't create structured mesh; error message:" 
                  << std::endl;
        PRINT_LAST_ERROR
        continue;
      }
    }
    else if (2 == i && argc > 4) {
        // read a file in parallel from the filename on the command line

      tmp_result = read_file(mbImpl, argv[4]);
      if (MB_SUCCESS != tmp_result) {
        result = tmp_result;
        std::cerr << "Couldn't read mesh; error message:" << std::endl;
        PRINT_LAST_ERROR
        continue;
      }
      nshared = -1;
    }

    if (MB_SUCCESS == tmp_result) {
        // now figure out which vertices are shared
      double wtime;
      if (0 == rank) wtime = MPI_Wtime();
      MBParallelComm *pcomm = new MBParallelComm(mbImpl);
      tmp_result = pcomm->resolve_shared_ents();
      if (MB_SUCCESS != tmp_result) {
        std::cerr << "Couldn't resolve shared entities; error message:" << std::endl;
        PRINT_LAST_ERROR
        result = tmp_result;
        continue;
      }
      
      if (0 == rank) wtime = MPI_Wtime() - wtime;

      MBRange shared_ents;
      tmp_result = pcomm->get_shared_entities(0, shared_ents);
      
      if (MB_SUCCESS != tmp_result) {
        std::cerr << "get_shared_entities returned error on proc " 
                  << rank << "; message: " << std::endl;
        PRINT_LAST_ERROR
        result = tmp_result;
      }
  
        // check # shared entities
      else if (0 <= nshared && nshared != (int) shared_ents.size()) {
        std::cerr << "Didn't get correct number of shared vertices on "
                  << "processor " << rank << std::endl;
        result = MB_FAILURE;
      }
      else if (i < 2) {
        std::cerr << "Proc " << rank;
        if (0 == i) std::cerr << " linear mesh succeeded." << std::endl;
        else std::cerr << " structured mesh succeeded." << std::endl;
        if (0 == rank) std::cerr << "   Time = " << wtime << "." << std::endl;
      }
      else {
        std::cerr << "Proc " << rank << " " << shared_ents.size()
                  << " shared entities." << std::endl;
      }
  
      delete pcomm;
      tmp_result = mbImpl->delete_mesh();
      if (MB_SUCCESS != tmp_result) {
        result = tmp_result;
        std::cerr << "Couldn't delete mesh on rank " << rank
                  << "; error message: " << std::endl;
        PRINT_LAST_ERROR
      }
    }
  }
  
  
  err = MPI_Finalize();

  if (MB_SUCCESS == result)
    std::cerr << "Proc " << rank << ": Success." << std::endl;
    
  return (MB_SUCCESS == result ? 0 : 1);
}

MBErrorCode read_file(MBInterface *mbImpl, const char *filename) 
{
  std::string options = "PARALLEL=BCAST_DELETE;PARTITION=MATERIAL_SET";
  MBEntityHandle file_set;
  MBErrorCode result = mbImpl->load_file(filename, file_set, 
                                         options.c_str());
  return result;
}

#define RR(a, b) {std::cerr << a; return b;}

MBErrorCode create_linear_mesh(MBInterface *mbImpl,
                               int N, int M, int &nshared) 
{
    /* create a mesh where each processor owns N vertices and shares
     * M vertices with processors on either end; for non-root procs,
     * that works out to N+M vertices.
     *
     * Number of vertices shared should be M on end procs, 2M on others
     */
  int my_rank = mbImpl->proc_rank();
  int my_size = mbImpl->proc_size();
  
  int nverts = N + M;
  if (0 == my_rank) nverts = N;

    // create my vertices and give them the right global ids
  MBRange my_verts;
  std::vector<double> coords(3*(nverts));
  std::fill(coords.begin(), coords.end(), 0.0);
  MBErrorCode result = mbImpl->create_vertices(&coords[0], nverts,
                                               my_verts);
  if (MB_SUCCESS != 0)
    RR("Failed to create vertices.", MB_FAILURE);
  
  std::vector<int> global_ids(N+M);
  for (int i = 0; i < nverts; i++)
    global_ids[i] = my_rank*N - (nverts-N) + i;
  
  int def_val = -1;
  MBTag gid_tag;
  result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE,
                              MB_TYPE_INTEGER, gid_tag,
                              &def_val, true);
  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) 
    RR("Failed to create tag.", MB_FAILURE);
  
  result = mbImpl->tag_set_data(gid_tag, my_verts, &global_ids[0]);
  if (MB_SUCCESS != result) RR("Failed to set global_id tag.", MB_FAILURE);

  nshared = ((my_rank == 0 || my_rank == my_size-1) ? M : 2*M);
  
  return MB_SUCCESS;
}

double LENGTH = 1.0;

void build_coords(const int nelem, double *&coords);

MBErrorCode create_scd_mesh(MBInterface *mbImpl,
                            int IJK, int &nshared) 
{
    /* Create a 3d mesh of hexes, parameterized in i, j, k, where
     * on processor of rank R the K parameterization starts at 
     * i=0, j=0, k=RK, and there are I+1, J+1, K+1 vertices in the
     * i, j, k directions, resp.; for now, make I, J, K equal
     *
     * Should share (I+1)(J+1) vertices on end procs, or twice that
     * on interior procs.
     */
  int my_rank = mbImpl->proc_rank();
  
    // make a 3d block of vertices
  MBEntitySequence *dum_seq = NULL;
  ScdVertexSeq *vseq = NULL;
  ScdElementSeq *eseq = NULL;
  EntitySequenceManager *seq_mgr = 
    dynamic_cast<MBCore*>(mbImpl)->sequence_manager();
  HomCoord vseq_minmax[2] = {HomCoord(0,0,0), 
                             HomCoord(IJK, IJK, IJK)};
  MBEntityHandle vstart, estart;
  
  MBErrorCode result = 
    seq_mgr->create_scd_sequence(vseq_minmax[0], vseq_minmax[1],
                                 MBVERTEX, 1, vstart, dum_seq);
  if (NULL != dum_seq) vseq = dynamic_cast<ScdVertexSeq*>(dum_seq);
  assert (MB_FAILURE != result && vstart != 0 && dum_seq != NULL && 
          vseq != NULL);
    // now the element sequence
  result = seq_mgr->create_scd_sequence(vseq_minmax[0], vseq_minmax[1], 
                                        MBHEX, 1, estart, dum_seq);
  if (NULL != dum_seq) eseq = dynamic_cast<ScdElementSeq*>(dum_seq);
  assert (MB_FAILURE != result && estart != 0 && dum_seq != NULL && 
          eseq != NULL);
  
    // only need to add one vseq to this, unity transform
    // trick: if I know it's going to be unity, just input 3 sets of equivalent points
  result = eseq->add_vsequence(vseq, vseq_minmax[0], vseq_minmax[0], 
                               vseq_minmax[0], 
                               vseq_minmax[0], vseq_minmax[0], vseq_minmax[0]);
  assert(MB_SUCCESS == result);

    // set the coordinates of the vertices
  double *coords = NULL;
  build_coords(IJK, coords);

    // offset coords by starting z distance
  double deltaz = my_rank*LENGTH;
  int num_verts = (IJK + 1)*(IJK + 1)*(IJK + 1);
  for (int k = 0; k < num_verts; k++)
    coords[3*k+2] += deltaz;
  
  MBRange vrange(vstart, vstart+num_verts-1);
  result = mbImpl->set_coords(vrange, coords);
  if (MB_SUCCESS != result) RR("Couldn't build coords array.", MB_FAILURE);
  
    // set global ids for vertices; reuse coords space
  int *gids = (int*) coords;
  int start_gid = 1 + my_rank * (IJK)*(IJK+1)*(IJK+1);
  for (int i = 0; i < num_verts; i++)
    gids[i] = start_gid++;
  
  int def_val = -1;
  MBTag gid_tag;
  result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE,
                              MB_TYPE_INTEGER, gid_tag,
                              &def_val, true);
  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) 
    RR("Failed to create tag.", MB_FAILURE);
  
  result = mbImpl->tag_set_data(gid_tag, vrange, gids);
  if (MB_SUCCESS != result) RR("Failed to set global_id tag.", MB_FAILURE);

  nshared = (IJK+1) * (IJK+1);
  
  return result;
}

void compute_edge(double *start, const int nelem,  const double xint,
                  const int stride) 
{
  for (int i = 1; i < nelem; i++) {
    start[i*stride] = start[0]+i*xint;
    start[nelem+1+i*stride] = start[nelem+1]+i*xint;
    start[2*(nelem+1)+i*stride] = start[2*(nelem+1)]+i*xint;
  }
}

void compute_face(double *a, const int nelem,  const double xint,
                  const int stride1, const int stride2) 
{
    // 2D TFI on a face starting at a, with strides stride1 in ada and stride2 in tse
  for (int j = 1; j < nelem; j++) {
    double tse = j * xint;
    for (int i = 1; i < nelem; i++) {
      double ada = i * xint;
      
      a[i*stride1+j*stride2] = (1.0 - ada)*a[i*stride1]
        + ada*a[i*stride1+nelem*stride2]
        + (1.0 - tse)*a[j*stride2]
        + tse*a[j*stride2+nelem*stride1]
        - (1.0 - tse)*(1.0 - ada)*a[0]
        - (1.0 - tse)*ada*a[nelem*stride1]
        - tse*(1.0 - ada)*a[nelem*stride2]
        - tse*ada*a[nelem*(stride1+stride2)];
      a[nelem+1+i*stride1+j*stride2] = (1.0 - ada)*a[nelem+1+i*stride1]
        + ada*a[nelem+1+i*stride1+nelem*stride2]
        + (1.0 - tse)*a[nelem+1+j*stride2]
        + tse*a[nelem+1+j*stride2+nelem*stride1]
        - (1.0 - tse)*(1.0 - ada)*a[nelem+1+0]
        - (1.0 - tse)*ada*a[nelem+1+nelem*stride1]
        - tse*(1.0 - ada)*a[nelem+1+nelem*stride2]
        - tse*ada*a[nelem+1+nelem*(stride1+stride2)];
      a[2*(nelem+1)+i*stride1+j*stride2] = (1.0 - ada)*a[2*(nelem+1)+i*stride1]
        + ada*a[2*(nelem+1)+i*stride1+nelem*stride2]
        + (1.0 - tse)*a[2*(nelem+1)+j*stride2]
        + tse*a[2*(nelem+1)+j*stride2+nelem*stride1]
        - (1.0 - tse)*(1.0 - ada)*a[2*(nelem+1)+0]
        - (1.0 - tse)*ada*a[2*(nelem+1)+nelem*stride1]
        - tse*(1.0 - ada)*a[2*(nelem+1)+nelem*stride2]
        - tse*ada*a[2*(nelem+1)+nelem*(stride1+stride2)];
    }
  }
}

void build_coords(const int nelem, double *&coords) 
{
    // allocate the memory
  int numv = nelem+1;
  int numv_sq = numv*numv;
  int tot_numv = numv*numv*numv;
  coords = new double[3*tot_numv];

// use FORTRAN-like indexing
#define VINDEX(i,j,k) (i + (j*numv) + (k*numv_sq))
  int idx;
  double scale1, scale2, scale3;
    // use these to prevent optimization on 1-scale, etc (real map wouldn't have
    // all these equal)
  scale1 = LENGTH/nelem;
  scale2 = LENGTH/nelem;
  scale3 = LENGTH/nelem;

#ifdef REALTFI
    // use a real TFI xform to compute coordinates
    // compute edges
    // i (stride=1)
  compute_edge(&coords[VINDEX(0,0,0)], nelem, scale1, 1);
  compute_edge(&coords[VINDEX(0,nelem,0)], nelem, scale1, 1);
  compute_edge(&coords[VINDEX(0,0,nelem)], nelem, scale1, 1);
  compute_edge(&coords[VINDEX(0,nelem,nelem)], nelem, scale1, 1);
    // j (stride=numv)
  compute_edge(&coords[VINDEX(0,0,0)], nelem, scale1, numv);
  compute_edge(&coords[VINDEX(nelem,0,0)], nelem, scale1, numv);
  compute_edge(&coords[VINDEX(0,0,nelem)], nelem, scale1, numv);
  compute_edge(&coords[VINDEX(nelem,0,nelem)], nelem, scale1, numv);
    // k (stride=numv^2)
  compute_edge(&coords[VINDEX(0,0,0)], nelem, scale1, numv_sq);
  compute_edge(&coords[VINDEX(nelem,0,0)], nelem, scale1, numv_sq);
  compute_edge(&coords[VINDEX(0,nelem,0)], nelem, scale1, numv_sq);
  compute_edge(&coords[VINDEX(nelem,nelem,0)], nelem, scale1, numv_sq);

    // compute faces
    // i=0, nelem
  compute_face(&coords[VINDEX(0,0,0)], nelem, scale1, numv, numv_sq);
  compute_face(&coords[VINDEX(nelem,0,0)], nelem, scale1, numv, numv_sq);
    // j=0, nelem
  compute_face(&coords[VINDEX(0,0,0)], nelem, scale1, 1, numv_sq);
  compute_face(&coords[VINDEX(0,nelem,0)], nelem, scale1, 1, numv_sq);
    // k=0, nelem
  compute_face(&coords[VINDEX(0,0,0)], nelem, scale1, 1, numv);
  compute_face(&coords[VINDEX(0,0,nelem)], nelem, scale1, 1, numv);

    // initialize corner indices
  int i000 = VINDEX(0,0,0);
  int ia00 = VINDEX(nelem,0,0);
  int i0t0 = VINDEX(0,nelem,0);
  int iat0 = VINDEX(nelem,nelem,0);
  int i00g = VINDEX(0,0,nelem);
  int ia0g = VINDEX(nelem,0,nelem);
  int i0tg = VINDEX(0,nelem,nelem);
  int iatg = VINDEX(nelem,nelem,nelem);
  double cX, cY, cZ;
  int adaInts = nelem;
  int tseInts = nelem;
  int gammaInts = nelem;
  
  
  for (int i=1; i < nelem; i++) {
    for (int j=1; j < nelem; j++) {
      for (int k=1; k < nelem; k++) {
        idx = VINDEX(i,j,k);
        double tse = i*scale1;
        double ada = j*scale2;
        double gamma = k*scale3;
        double tm1 = 1.0 - tse;
        double am1 = 1.0 - ada;
        double gm1 = 1.0 - gamma;

        cX = gm1 *   (am1*(tm1*coords[i000] + tse*coords[i0t0])  +
                      ada*(tm1*coords[ia00] + tse*coords[iat0])) +
          gamma * (am1*(tm1*coords[i00g] + tse*coords[i0tg])  +
                   ada*(tm1*coords[ia0g] + tse*coords[iatg]));

        cY = gm1 *   (am1*(tm1*coords[i000] + tse*coords[i0t0])  +
                      ada*(tm1*coords[ia00] + tse*coords[iat0])) +
          gamma * (am1*(tm1*coords[i00g] + tse*coords[i0tg])  +
                   ada*(tm1*coords[ia0g] + tse*coords[iatg]));

        cZ = gm1 *   (am1*(tm1*coords[i000] + tse*coords[i0t0])  +
                      ada*(tm1*coords[ia00] + tse*coords[iat0])) +
          gamma * (am1*(tm1*coords[i00g] + tse*coords[i0tg])  +
                   ada*(tm1*coords[ia0g] + tse*coords[iatg]));

        double *ai0k = &coords[VINDEX(k,0,i)];
        double *aiak = &coords[VINDEX(k,adaInts,i)];
        double *a0jk = &coords[VINDEX(k,j,0)];
        double *atjk = &coords[VINDEX(k,j,tseInts)];
        double *aij0 = &coords[VINDEX(0,j,i)];
        double *aijg = &coords[VINDEX(gammaInts,j,i)];
  
        coords[VINDEX(i,j,k)] = (   am1*ai0k[0] 
                                    + ada*aiak[0] 
                                    + tm1*a0jk[0] 
                                    + tse*atjk[0]
                                    + gm1*aij0[0] 
                                    + gamma*aijg[0] )/2.0 - cX/2.0;

        coords[nelem+1+VINDEX(i,j,k)] = (   am1*ai0k[nelem+1] 
                                            + ada*aiak[nelem+1] 
                                            + tm1*a0jk[nelem+1] 
                                            + tse*atjk[nelem+1]
                                            + gm1*aij0[nelem+1] 
                                            + gamma*aijg[nelem+1] )/2.0 - cY/2.0;

        coords[2*(nelem+1)+VINDEX(i,j,k)] = (   am1*ai0k[2*(nelem+1)] 
                                                + ada*aiak[2*(nelem+1)] 
                                                + tm1*a0jk[2*(nelem+1)] 
                                                + tse*atjk[2*(nelem+1)]
                                                + gm1*aij0[2*(nelem+1)] 
                                                + gamma*aijg[2*(nelem+1)] )/2.0 - cZ/2.0;
      }
    }
  }
  

#else
  for (int i=0; i < numv; i++) {
    for (int j=0; j < numv; j++) {
      for (int k=0; k < numv; k++) {
        idx = VINDEX(i,j,k);
          // blocked coordinate ordering
        coords[idx] = i*scale1;
        coords[tot_numv+idx] = j*scale2;
        coords[2*tot_numv+idx] = k*scale3;
      }
    }
  }
#endif
}

void build_connect(const int nelem, const MBEntityHandle vstart, MBEntityHandle *&connect) 
{
    // allocate the memory
  int nume_tot = nelem*nelem*nelem;
  connect = new MBEntityHandle[8*nume_tot];

  MBEntityHandle vijk;
  int numv = nelem + 1;
  int numv_sq = numv*numv;
  int idx = 0;
  for (int i=0; i < nelem; i++) {
    for (int j=0; j < nelem; j++) {
      for (int k=0; k < nelem; k++) {
        vijk = vstart+VINDEX(i,j,k);
        connect[idx++] = vijk;
        connect[idx++] = vijk+1;
        connect[idx++] = vijk+1+numv;
        connect[idx++] = vijk+numv;
        connect[idx++] = vijk+numv*numv;
        connect[idx++] = vijk+1+numv*numv;
        connect[idx++] = vijk+1+numv+numv*numv;
        connect[idx++] = vijk+numv+numv*numv;
        assert(i <= numv*numv*numv);
      }
    }
  }
}

