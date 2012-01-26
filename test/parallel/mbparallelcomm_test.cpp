/** test of ParallelComm functionality
 *
 * To run:
 *
 * mpirun -np <#procs> mbparallelcomm_test
 *
 */

#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "ReadParallel.hpp"
#include "FileOptions.hpp"
#include "MBTagConventions.hpp"
#include "moab/Core.hpp"
#include "ScdVertexData.hpp"
#include "StructuredElementSeq.hpp"
#include "SequenceManager.hpp"
#include "moab/Error.hpp"
#include "moab_mpi.h"
#include <iostream>
#include <sstream>
#include <assert.h>

#define REALTFI 1

const bool debug = false;

using namespace moab;

#define ERROR(a, b) {std::cerr << a << std::endl; return b;}

#define PRINT_LAST_ERROR {\
        std::string last_error;\
        result = mbImpl->get_last_error(last_error);\
        if (last_error.empty()) std::cerr << "(none)" << std::endl;\
        else std::cerr << last_error << std::endl;\
        }
#define RRA(a) if (MB_SUCCESS != result) {\
    std::cerr << a; return result;}
    

ErrorCode create_linear_mesh(Interface *mbImpl,
                               int N, int M, int &nshared);

ErrorCode create_scd_mesh(Interface *mbImpl,
                            int IJK, int &nshared);

ErrorCode read_file(Interface *mbImpl, std::vector<std::string> &filenames,
                      const char *tag_name, int tag_val, int distrib,
                      int parallel_option, int resolve_shared, int with_ghosts, 
                      int use_mpio, bool print_parallel);

ErrorCode test_packing(Interface *mbImpl, const char *filename);

ErrorCode report_nsets(Interface *mbImpl);

ErrorCode report_iface_ents(Interface *mbImpl,
                              std::vector<ParallelComm *> &pcs);

void print_usage(const char *);

int main(int argc, char **argv) 
{
    // need to init MPI first, to tell how many procs and rank
  int err = MPI_Init(&argc, &argv);

  int nprocs, rank;
  err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // start time
  double stime = 0, rtime = 0, dtime = 0, ltime = 0;
  if (0 == rank) stime = MPI_Wtime();

    // create MOAB instance based on that
  Interface *mbImpl = new Core;
  if (NULL == mbImpl) return 1;
  
  ErrorCode result = MB_SUCCESS;

    // each interior proc has a vector of N+M vertices, sharing
    // M vertices each with lower- and upper-rank processors, except
    // procs on the end

    // get N, M from command line
  if (argc < 3) {
    if (0 == rank) print_usage(argv[0]);
    err = MPI_Finalize();
    return 1;
  }

  int npos = 1, tag_val, distrib, with_ghosts = 1, resolve_shared = 1, 
      use_mpio = 0;
  bool print_parallel = false;
  const char *tag_name;
  std::vector<std::string> filenames;
  int parallel_option = 0;
  int num_files;

  if (!strcmp(argv[npos], "-p")) print_parallel = true;

  while (npos != argc) {    
    ErrorCode tmp_result;
    int nshared = -1;
    int this_opt = strtol(argv[npos++], NULL, 0);
    switch (this_opt) {
      case 0:
      case -1:
      case -2:
      case -3:
          parallel_option = this_opt;
          continue;
        
      case 3:
            // read a file in parallel from the filename on the command line
          tag_name = "MATERIAL_SET";
          tag_val = -1;
          num_files = strtol(argv[npos++], NULL, 0);
          if (0 == num_files) {
            if (0 == rank) print_usage(argv[0]);
            err = MPI_Finalize();
            return 1;
          }
          while (num_files-- && npos < argc)
            filenames.push_back(std::string(argv[npos++]));
          if (npos < argc) tag_name = argv[npos++];
          if (npos < argc) tag_val = strtol(argv[npos++], NULL, 0);
          if (npos < argc) distrib = strtol(argv[npos++], NULL, 0);
          else distrib = 1;
          if (npos < argc) resolve_shared = strtol(argv[npos++], NULL, 0);
          if (npos < argc) with_ghosts = strtol(argv[npos++], NULL, 0);
          if (npos < argc) use_mpio = strtol(argv[npos++], NULL, 0);

          tmp_result = read_file(mbImpl, filenames, tag_name, tag_val,
                                 distrib, parallel_option, 
                                 resolve_shared, with_ghosts, use_mpio,
                                 print_parallel);
          if (MB_SUCCESS != tmp_result) {
            result = tmp_result;
            std::cerr << "Couldn't read mesh; error message:" << std::endl;
            PRINT_LAST_ERROR;
            MPI_Abort(MPI_COMM_WORLD, result);
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
                                 with_ghosts, use_mpio, print_parallel);
          if (MB_SUCCESS != tmp_result) {
            result = tmp_result;
            std::cerr << "Couldn't read mesh; error message:" << std::endl;
            PRINT_LAST_ERROR;
            MPI_Abort(MPI_COMM_WORLD, result);
          }
          nshared = -1;
          break;

      default:
          std::cerr << "Unrecognized option \"" << this_opt
                    << "\"; skipping." << std::endl;
          tmp_result = MB_FAILURE;
    }
    

    if (0 == rank) rtime = MPI_Wtime();
  }
  
  if (0 == rank) dtime = MPI_Wtime();

  result = mbImpl->delete_mesh();
  if (MB_SUCCESS != result) {
    std::cerr << "Couldn't delete mesh on rank " << rank
              << "; error message: " << std::endl;
    PRINT_LAST_ERROR;
  }
  if (0 == rank) ltime = MPI_Wtime();

  if (MB_SUCCESS == result)
    std::cerr << "Proc " << rank << ": Success." << std::endl;
    
  if (0 == rank) std::cout << "Times: " 
                           << dtime-stime << " "
                           << rtime-stime << " "
                           << ltime - dtime
                           << " (total/read/delete)"
                           << std::endl;

  err = MPI_Finalize();

  delete mbImpl;
  
  return (MB_SUCCESS == result ? 0 : 1);
}

void print_usage(const char *command) 
{
  std::cerr 
      << "Usage: " << command
      << " [readpar_option] <opt> <input> [...] where:" << std::endl
      << " readpar_option = 0 (BCAST_DELETE) (default), -1 (READ_DELETE), " << std::endl
      << "                 -2 (READ_PARALLEL), -3 (BCAST)" << std::endl
      << "opt   input" << std::endl
      << "===   =====" << std::endl
      << " 1     <linear_ints> <shared_verts> " << std::endl
      << " 2     <n_ints> " << std::endl
      << " 3*    <# files> <file_names...> [<tag_name>=\"MATERIAL_SET\" [tag_val] [distribute=1] [resolve_shared=1] [with_ghosts=1] [use_mpio=0]" << std::endl
      << " 4    <file_name> " << std::endl
      << "*Note: if opt 3 is used, it must be the last one." << std::endl;
}

ErrorCode report_nsets(Interface *mbImpl) 
{
    // get and report various numbers...
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  Range matsets, geomsets, parsets;
  int nsets;
  Tag mtag = 0, gtag = 0, ptag = 0, gidtag;
  ErrorCode result = mbImpl->tag_get_handle("MATERIAL_SET", 1, MB_TYPE_INTEGER, mtag);
  result = mbImpl->tag_get_handle("GEOM_DIMENSION", 1, MB_TYPE_INTEGER, gtag);
  result = mbImpl->tag_get_handle("PARALLEL_PARTITION", 1, MB_TYPE_INTEGER, ptag);
  result = mbImpl->tag_get_handle("GLOBAL_ID", 1, MB_TYPE_INTEGER, gidtag);

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

ErrorCode read_file(Interface *mbImpl, 
                      std::vector<std::string> &filenames,
                      const char *tag_name, int tag_val,
                      int distrib, int parallel_option, int resolve_shared,
                      int with_ghosts, int use_mpio, bool print_parallel) 
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
      options << "PARALLEL=READ_PART;PARTITION=" << tag_name;
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

  if (1 == use_mpio)
    options << ";USE_MPIO";

  options << ";CPUTIME";

  if (print_parallel) 
    options << ";PRINT_PARALLEL";

  std::vector<ParallelComm*> pcs(filenames.size());
  std::vector<ReadParallel*> rps(filenames.size());
  ErrorCode result;

  if (1 < filenames.size()) {
    for (unsigned int i = 0; i < filenames.size(); i++) {
      pcs[i] = new ParallelComm(mbImpl);
      rps[i] = new ReadParallel(mbImpl, pcs[i]);
    
      result = rps[i]->load_file(filenames[i].c_str(), 0, 
                                 FileOptions(options.str().c_str()), 0, 0);
      if (MB_SUCCESS != result) 
        PRINT_LAST_ERROR;

      if (MB_SUCCESS != result) {
        MPI_Abort(MPI_COMM_WORLD, result);
        break;
      }

        // exchange tag
      Range tmp_range;
      result = pcs[i]->exchange_tags("GLOBAL_ID", tmp_range);
      if (MB_SUCCESS != result) {
        std::cerr << "Tag exchange didn't work." << std::endl;
        break;
      }

    }
  }
  else {
    result = mbImpl->load_file(filenames[0].c_str(), 0, 
                               options.str().c_str());
    RRA("Failed to load file.");
    pcs[0] = ParallelComm::get_pcomm(mbImpl, 0);
    assert(pcs[0]);
  }
    
  if (MB_SUCCESS == result) report_iface_ents(mbImpl, pcs);
  
  return result;
}

ErrorCode test_packing(Interface *mbImpl, const char *filename) 
{
    // read the mesh
  EntityHandle file_set;
  ErrorCode result = mbImpl->create_meshset( MESHSET_SET, file_set );
  RRA("create_meshset failed.");

  result = mbImpl->load_file(filename, &file_set, NULL);
  if (MB_SUCCESS != result) {
    std::cerr << "Reading file failed; message:" << std::endl;
    PRINT_LAST_ERROR;
    return result;
  }
  
    // get 3d entities and pack a buffer with them
  Range ents, whole_range;
  std::vector<EntityHandle> new_ents;
  result = mbImpl->get_entities_by_handle(file_set, ents);
  RRA("Getting 3d ents failed.");
  
  ents.insert(file_set);
  
  ParallelComm *pcomm = new ParallelComm(mbImpl);

  ParallelComm::Buffer buff;
  result = pcomm->pack_buffer(ents, false, true, false, -1, &buff);
  RRA("Packing buffer count (non-stored handles) failed.");

  std::vector<std::vector<EntityHandle> > L1hloc, L1hrem;
  std::vector<std::vector<int> > L1p;
  std::vector<EntityHandle> L2hloc, L2hrem;
  std::vector<unsigned int> L2p;
  
  buff.reset_ptr();
  result = pcomm->unpack_buffer(buff.buff_ptr, false, -1, -1, L1hloc, L1hrem, L1p, L2hloc, 
                         L2hrem, L2p, new_ents);
  RRA("Unpacking buffer (non-stored handles) failed.");

  return MB_SUCCESS;
}

ErrorCode report_iface_ents(Interface *mbImpl,
                              std::vector<ParallelComm *> &pcs) 
{
  Range iface_ents[6];
  ErrorCode result = MB_SUCCESS, tmp_result;
  
    // now figure out which vertices are shared
  Range part_ents, part_verts;
  for (unsigned int p = 0; p < pcs.size(); p++) {
    // get entities owned by this partition
    for (Range::iterator rit = pcs[p]->partition_sets().begin();
	 rit != pcs[p]->partition_sets().end(); rit++) {
      tmp_result = mbImpl->get_entities_by_dimension(*rit, 3, part_ents, true);
      if (MB_SUCCESS != tmp_result) result = tmp_result;
    }

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

    // get non-owned vertices
  result = pcs[0]->get_pstatus_entities(0, PSTATUS_NOT_OWNED, part_verts);
  if (MB_SUCCESS != result) {
    std::cerr << "Couldn't get non-owned entities." << std::endl;
    return result;
  }
  int tot_verts;
  result = mbImpl->get_number_entities_by_dimension(0, 0, tot_verts);
  if (MB_SUCCESS != result) {
    std::cerr << "Couldn't get number of vertices." << std::endl;
    return result;
  }
  tot_verts -= part_verts.size();

    // report # iface entities
  result = mbImpl->get_adjacencies(iface_ents[4], 0, false, iface_ents[5], 
                                   Interface::UNION);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  std::cerr << "Proc " << rank << " iface entities: " << std::endl;
  for (int i = 0; i < 4; i++)
    std::cerr << "    " << iface_ents[i].size() << " "
              << i << "d iface entities." << std::endl;
  std::cerr << "    (" << iface_ents[5].size() 
            << " verts adj to other iface ents)" << std::endl;
  if (iface_ents[0].size() != iface_ents[5].size())
    std::cerr << "WARNING: number of interface vertices don't agree with "
	      << "vertex adjacencies on interface entities." << std::endl;

  // report # regions owned by this proc
  std::cout << "Proc " << rank << " owns " << part_ents.size() 
	    << " 3d entities." << std::endl;

    // get total # regions over all procs
  int num_local[2], num_total[2];
  num_local[0] = tot_verts;
  num_local[1] = part_ents.size();
  
  int failure = MPI_Reduce(&num_local, &num_total, 2,
                           MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (failure) result = MB_FAILURE;

  if (0 == rank) {
    std::cout << "Total # owned vertices = " << num_total[0] << std::endl;
    std::cout << "Total # owned regions = " << num_total[1] << std::endl;
  }
  
  return result;
}

