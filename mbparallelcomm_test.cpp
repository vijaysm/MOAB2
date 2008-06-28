/** test of MBParallelComm functionality
 *
 * To run:
 *
 * mpirun -np <#procs> mbparallelcomm_test
 *
 */

#include "MBParallelComm.hpp"
#include "MBParallelConventions.h"
#include "ReadParallel.hpp"
#include "FileOptions.hpp"
#include "MBTagConventions.hpp"
#include "MBCore.hpp"
#include "ScdVertexData.hpp"
#include "StructuredElementSeq.hpp"
#include "SequenceManager.hpp"
#include "MBError.hpp"
#include "mpi.h"
#include <iostream>
#include <sstream>
#include <assert.h>

#define REALTFI 1

const bool debug = false;

#define ERROR(a, b) {std::cerr << a << std::endl; return b;}

#define PRINT_LAST_ERROR {\
        std::string last_error;\
        result = mbImpl->get_last_error(last_error);\
        if (last_error.empty()) std::cerr << "(none)" << std::endl;\
        else std::cerr << last_error << std::endl;\
        }
#define RRA(a) if (MB_SUCCESS != result) {\
      std::string tmp_str; mbImpl->get_last_error(tmp_str);\
      tmp_str.append("\n"); tmp_str.append(a);\
      dynamic_cast<MBCore*>(mbImpl)->get_error_handler()->set_last_error(tmp_str.c_str()); \
      return result;}

MBErrorCode create_linear_mesh(MBInterface *mbImpl,
                               int N, int M, int &nshared);

MBErrorCode create_scd_mesh(MBInterface *mbImpl,
                            int IJK, int &nshared);

MBErrorCode read_file(MBInterface *mbImpl, std::vector<std::string> &filenames,
                      const char *tag_name, int tag_val, int distrib,
                      int parallel_option, int resolve_shared, int with_ghosts);

MBErrorCode test_packing(MBInterface *mbImpl, const char *filename);

MBErrorCode report_nsets(MBInterface *mbImpl);

MBErrorCode report_iface_ents(MBInterface *mbImpl,
                              std::vector<MBParallelComm *> &pcs);

int main(int argc, char **argv) 
{
    // need to init MPI first, to tell how many procs and rank
  int err = MPI_Init(&argc, &argv);

  int nprocs, rank;
  err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // start time
  double stime, rtime, setime, dtime, ltime;
  if (0 == rank) stime = MPI_Wtime();

    // create MOAB instance based on that
  MBInterface *mbImpl = new MBCore(rank, nprocs);
  if (NULL == mbImpl) return 1;
  
  MBErrorCode result = MB_SUCCESS;

    // each interior proc has a vector of N+M vertices, sharing
    // M vertices each with lower- and upper-rank processors, except
    // procs on the end

    // get N, M from command line
  int N, M;
  if (argc < 3) {
    if (0 == rank)
      std::cerr 
        << "Usage: " << argv[0] 
        << " [readpar_option] <opt> <input> [...] where:" << std::endl
        << " readpar_option = 0 (BCAST_DELETE) (default), -1 (READ_DELETE), " << std::endl
        << "                 -2 (READ_PARALLEL), -3 (BCAST)" << std::endl
        << "opt   input" << std::endl
        << "===   =====" << std::endl
        << " 1     <linear_ints> <shared_verts> " << std::endl
        << " 2     <n_ints> " << std::endl
        << " 3*    <file_name> [<tag_name>=\"MATERIAL_SET\" [tag_val] [distribute=1] [resolve_shared=1] [with_ghosts=1]" << std::endl
        << " 4    <file_name> " << std::endl
        << "*Note: if opt 3 is used, it must be the last one." << std::endl;
    
    err = MPI_Finalize();
    return 1;
  }

  int npos = 1, tag_val, distrib, with_ghosts = 1, resolve_shared = 1;
  const char *tag_name;
  std::vector<std::string> filenames;
  int parallel_option = 0;

  while (npos < argc) {
    MBErrorCode tmp_result;
    int nshared = -1;
    int this_opt = strtol(argv[npos++], NULL, 0);
    switch (this_opt) {
      case 0:
      case -1:
      case -2:
      case -3:
        parallel_option = this_opt;
        continue;
        
      case 1:
        N = atoi(argv[npos++]);
        M = atoi(argv[npos++]);
        tmp_result = create_linear_mesh(mbImpl, N, M, nshared);
        if (MB_SUCCESS != tmp_result) {
          result = tmp_result;
          std::cerr << "Couldn't create linear mesh; error message:." 
                    << std::endl;
          PRINT_LAST_ERROR
        }
        break;
      case 2:
        N = atoi(argv[npos++]);
        tmp_result = create_scd_mesh(mbImpl, N, nshared);
        if (MB_SUCCESS != tmp_result) {
          result = tmp_result;
          std::cerr << "Couldn't create structured mesh; error message:" 
                    << std::endl;
          PRINT_LAST_ERROR
        }
        break;
          
      case 3:
          // read a file in parallel from the filename on the command line
        tag_name = "MATERIAL_SET";
        tag_val = -1;
        filenames.push_back(std::string(argv[npos++]));
        if (npos < argc) tag_name = argv[npos++];
        if (npos < argc) tag_val = strtol(argv[npos++], NULL, 0);
        if (npos < argc) distrib = strtol(argv[npos++], NULL, 0);
        else distrib = 1;
        if (npos < argc) resolve_shared = strtol(argv[npos++], NULL, 0);
        if (npos < argc) with_ghosts = strtol(argv[npos++], NULL, 0);

        tmp_result = read_file(mbImpl, filenames, tag_name, tag_val,
                               distrib, parallel_option, 
                               resolve_shared, with_ghosts);
        if (MB_SUCCESS != tmp_result) {
          result = tmp_result;
          std::cerr << "Couldn't read mesh; error message:" << std::endl;
          PRINT_LAST_ERROR
        }
        nshared = -1;
        break;

      case 4:
        filenames.push_back(argv[npos++]);
        tmp_result = test_packing(mbImpl, filenames[0].c_str());
        if (MB_SUCCESS != tmp_result) {
          result = tmp_result;
          std::cerr << "Packing test failed; error message:" << std::endl;
          PRINT_LAST_ERROR
        }
        break;

      case 5:
          // read a file in parallel from the filename on the command line
        tag_name = "MATERIAL_SET";
        distrib = 1;
        tag_val = -1;
        with_ghosts = 0;
        resolve_shared = 1;
        while (npos < argc)
          filenames.push_back(std::string(argv[npos++]));
        tmp_result = read_file(mbImpl, filenames, tag_name, tag_val,
                               distrib, parallel_option, resolve_shared,
                               with_ghosts);
        if (MB_SUCCESS != tmp_result) {
          result = tmp_result;
          std::cerr << "Couldn't read mesh; error message:" << std::endl;
          PRINT_LAST_ERROR
        }
        nshared = -1;
        break;

      default:
        std::cerr << "Unrecognized option \"" << this_opt
                  << "\"; skipping." << std::endl;
        tmp_result = MB_FAILURE;
    }
    

    if (0 == rank) rtime = MPI_Wtime();
    if (MB_SUCCESS == tmp_result && 4 != this_opt && false) {
        // now figure out which vertices are shared
      MBParallelComm *pcomm = MBParallelComm::get_pcomm(mbImpl, 0);
      assert(pcomm);

      MBRange iface_ents[7];
      for (int i = 0; i < 4; i++) {
        tmp_result = pcomm->get_iface_entities(-1, i, iface_ents[i]);
      
        if (MB_SUCCESS != tmp_result) {
          std::cerr << "get_iface_entities returned error on proc " 
                    << rank << "; message: " << std::endl;
          PRINT_LAST_ERROR;
          result = tmp_result;
        }
        if (0 != i) iface_ents[4].merge(iface_ents[i]);
      }
      result = pcomm->get_part_entities(iface_ents[6], -1);
      PRINT_LAST_ERROR;

      std::cerr << "Proc " << rank << " partition entities:" << std::endl;
      iface_ents[6].print("   ");
      
      if (0 == rank) setime = MPI_Wtime();

        // check # iface entities
      if (0 <= nshared && nshared != (int) iface_ents[0].size()) {
        std::cerr << "Didn't get correct number of iface vertices on "
                  << "processor " << rank << std::endl;
        result = MB_FAILURE;
      }

      else
        std::cerr << "Proc " << rank << " option " << this_opt
                << " succeeded." << std::endl;

      if (-1 == nshared) {
        result = mbImpl->get_adjacencies(iface_ents[4], 0, false, iface_ents[5], 
                                         MBInterface::UNION);
        
        std::cerr << "Proc " << rank << " iface entities: " << std::endl;
        for (int i = 0; i < 4; i++)
          std::cerr << "    " << iface_ents[i].size() << " "
                    << i << "d iface entities." << std::endl;
        std::cerr << "    (" << iface_ents[5].size() 
                  << " verts adj to other iface ents)" << std::endl;
      }
      
      if (debug && false) {
//      if (debug && 2 == nprocs) {
          // if I'm root, get and print handles on other procs
        std::vector<MBEntityHandle> sharedh_tags(iface_ents[0].size());
        std::fill(sharedh_tags.begin(), sharedh_tags.end(), 0);
        MBTag dumt, sharedh_tag;
        result = pcomm->get_shared_proc_tags(dumt, dumt, sharedh_tag, dumt, dumt);
        result = mbImpl->tag_get_data(sharedh_tag, iface_ents[0], &sharedh_tags[0]);
        if (MB_SUCCESS != result) {
          std::cerr << "Couldn't get shared handle tag." << std::endl;
        }
        else {
          MBRange dum_range;
          std::copy(sharedh_tags.begin(), sharedh_tags.end(), mb_range_inserter(dum_range));
          std::cerr << "Shared handles: " << std::endl;
          dum_range.print();
        }
        
      result = report_nsets(mbImpl);
      }

      if (0 == rank) ltime = MPI_Wtime();
  
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
  
  if (0 == rank) dtime = MPI_Wtime();

  err = MPI_Finalize();

  if (MB_SUCCESS == result)
    std::cerr << "Proc " << rank << ": Success." << std::endl;
    
  if (0 == rank) std::cout << "Times: " 
                           << dtime-stime << " "
                           << rtime-stime << " "
                           << setime-rtime << " "
                           << ltime-setime << " "
                           << dtime - ltime
                           << " (total/read/shared/report/delete)"
                           << std::endl;
   
  return (MB_SUCCESS == result ? 0 : 1);
}

MBErrorCode report_nsets(MBInterface *mbImpl) 
{
    // get and report various numbers...
  int rank = mbImpl->proc_rank();
  
  MBRange matsets, geomsets, parsets;
  int nsets;
  MBTag mtag = 0, gtag = 0, ptag = 0, gidtag;
  MBErrorCode result = mbImpl->tag_get_handle("MATERIAL_SET", mtag);
  result = mbImpl->tag_get_handle("GEOM_DIMENSION", gtag);
  result = mbImpl->tag_get_handle("PARALLEL_PARTITION", ptag);
  result = mbImpl->tag_get_handle("GLOBAL_ID", gidtag);

  result = mbImpl->get_number_entities_by_type(0, MBENTITYSET, nsets);
  std::cout << "Proc " << rank << ": Total of " << nsets
            << " entity sets." << std::endl;
  
#define PRINTSETS(a, b, c, p) \
  if (a) {\
    result = mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, & a,\
                                                  p, 1, b); \
    if (! b .empty()) {\
      std::vector<int> ids( b .size());\
      result = mbImpl->tag_get_data(gidtag, b, &ids[0]); \
      if (MB_SUCCESS == result) {\
        std::cout << "Proc " << rank << ": " << c \
          << " (total " << b.size() << "): " \
           << ids[0]; \
        for (unsigned int i = 1; i < b .size(); i++) \
          std::cout << ", " << ids[i]; \
        std::cout << std::endl; \
      } } }
  
  PRINTSETS(mtag, matsets, "material sets", NULL);
  
  int tval = 3;
  void *pval = &tval;
  
  PRINTSETS(gtag, geomsets, "geom sets (vols)", &pval);
  tval = 2;
  geomsets.clear();
  PRINTSETS(gtag, geomsets, "geom sets (surfs)", &pval);
  tval = 1;
  geomsets.clear();
  PRINTSETS(gtag, geomsets, "geom sets (curves)", &pval);
  tval = 0;
  geomsets.clear();
  PRINTSETS(gtag, geomsets, "geom sets (verts)", &pval);
  
  PRINTSETS(ptag, parsets, "partition sets", NULL);

  if (debug) {
      // list info on all ent sets, reuse parsets
    parsets.clear();
    result = mbImpl->get_entities_by_type(0, MBENTITYSET, parsets);
    if (MB_SUCCESS == result) {
      std::cout << "Total sets (by range): " << parsets.size() << "; sets: " << std::endl;
      parsets.print("  ");
      mbImpl->list_entities(parsets);
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode read_file(MBInterface *mbImpl, 
                      std::vector<std::string> &filenames,
                      const char *tag_name, int tag_val,
                      int distrib, int parallel_option, int resolve_shared,
                      int with_ghosts) 
{
  std::ostringstream options;
  switch (parallel_option) {
    case 0:
      options << "PARALLEL=BCAST_DELETE;PARTITION=" << tag_name;
      break;
    case -1:
      options << "PARALLEL=READ_DELETE;PARTITION=" << tag_name;
      break;
    case -2:
      options << "PARALLEL=READ_PARALLEL;PARTITION=" << tag_name;
      break;
    case -3:
      options << "PARALLEL=BCAST;PARTITION=" << tag_name;
      break;
    default:
      return MB_FAILURE;
  }
  
  if (-1 != tag_val)
    options << ";PARTITION_VAL=" << tag_val;

  if (1 == distrib)
    options << ";PARTITION_DISTRIBUTE";

  if (1 == resolve_shared)
    options << ";PARALLEL_RESOLVE_SHARED_ENTS";

  if (1 == with_ghosts)
    options << ";PARALLEL_GHOSTS=3.0.1";

  options << ";CPUTIME";

  std::vector<MBEntityHandle> filesets(filenames.size());
  std::vector<MBParallelComm*> pcs(filenames.size());
  std::vector<ReadParallel*> rps(filenames.size());
  MBErrorCode result;
  
  for (unsigned int i = 0; i < filenames.size(); i++) {
    pcs[i] = new MBParallelComm(mbImpl);
    rps[i] = new ReadParallel(mbImpl, pcs[i]);
    
    result = rps[i]->load_file(filenames[i].c_str(), filesets[i], 
                               FileOptions(options.str().c_str()), NULL, 0);
    PRINT_LAST_ERROR;
  }

  report_iface_ents(mbImpl, pcs);
  
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
  EntitySequence *dum_seq = NULL;
  ScdVertexData *vseq = NULL;
  StructuredElementSeq *eseq = NULL;
  SequenceManager *seq_mgr = 
    dynamic_cast<MBCore*>(mbImpl)->sequence_manager();
  HomCoord vseq_minmax[2] = {HomCoord(0,0,0), 
                             HomCoord(IJK, IJK, IJK)};
  MBEntityHandle vstart, estart;
  
  MBErrorCode result = 
    seq_mgr->create_scd_sequence(vseq_minmax[0], vseq_minmax[1],
                                 MBVERTEX, 1, -1, vstart, dum_seq);
  if (NULL != dum_seq) vseq = dynamic_cast<ScdVertexData*>(dum_seq->data());
  assert (MB_FAILURE != result && vstart != 0 && dum_seq != NULL && 
          vseq != NULL);
    // now the element sequence
  result = seq_mgr->create_scd_sequence(vseq_minmax[0], vseq_minmax[1], 
                                        MBHEX, 1, -1, estart, dum_seq);
  if (NULL != dum_seq) eseq = dynamic_cast<StructuredElementSeq*>(dum_seq);
  assert (MB_FAILURE != result && estart != 0 && dum_seq != NULL && 
          eseq != NULL);
  
    // only need to add one vseq to this, unity transform
    // trick: if I know it's going to be unity, just input 3 sets of equivalent points
  result = eseq->sdata()->add_vsequence(vseq, vseq_minmax[0], vseq_minmax[0], 
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

  if (my_rank != 0 && my_rank != mbImpl->proc_size()-1) nshared *= 2;
  
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

MBErrorCode test_packing(MBInterface *mbImpl, const char *filename) 
{
    // read the mesh
  MBEntityHandle file_set;
  MBErrorCode result = mbImpl->load_file(filename, file_set, NULL);
  if (MB_SUCCESS != result) {
    std::cerr << "Reading file failed; message:" << std::endl;
    PRINT_LAST_ERROR;
    return result;
  }
  
    // get 3d entities and pack a buffer with them
  MBRange ents, new_ents, whole_range;
  result = mbImpl->get_entities_by_handle(file_set, ents);
  RRA("Getting 3d ents failed.");
  
  ents.insert(file_set);
  
  MBParallelComm *pcomm = new MBParallelComm(mbImpl);
  std::vector<unsigned char> buff(1024);
  int buff_size;
  result = pcomm->pack_buffer(ents, false, true, false, false, -1,
                              whole_range, buff, buff_size);
  RRA("Packing buffer count (non-stored handles) failed.");

  result = pcomm->unpack_buffer(&buff[0], false, -1, new_ents);
  RRA("Unacking buffer (non-stored handles) failed.");

  return MB_SUCCESS;
}

MBErrorCode report_iface_ents(MBInterface *mbImpl,
                              std::vector<MBParallelComm *> &pcs) 
{
  MBRange iface_ents[6];
  MBErrorCode result = MB_SUCCESS, tmp_result;
  
    // now figure out which vertices are shared
  for (unsigned int p = 0; p < pcs.size(); p++) {
    for (int i = 0; i < 4; i++) {
      tmp_result = pcs[p]->get_iface_entities(-1, i, iface_ents[i]);
      
      if (MB_SUCCESS != tmp_result) {
        std::cerr << "get_iface_entities returned error on proc " 
                  << pcs[p]->proc_config().proc_rank() << "; message: " << std::endl;
        std::string last_error;
        result = mbImpl->get_last_error(last_error);
        if (last_error.empty()) std::cerr << "(none)" << std::endl;
        else std::cerr << last_error << std::endl;
        result = tmp_result;
      }
      if (0 != i) iface_ents[4].merge(iface_ents[i]);
    }
  }

    // report # iface entities
  result = mbImpl->get_adjacencies(iface_ents[4], 0, false, iface_ents[5], 
                                   MBInterface::UNION);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  std::cerr << "Proc " << rank << " iface entities: " << std::endl;
  for (int i = 0; i < 4; i++)
    std::cerr << "    " << iface_ents[i].size() << " "
              << i << "d iface entities." << std::endl;
  std::cerr << "    (" << iface_ents[5].size() 
            << " verts adj to other iface ents)" << std::endl;

  return result;
}

