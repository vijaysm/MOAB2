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
                      int parallel_option, int resolve_shared, int with_ghosts, int use_mpio);

MBErrorCode test_packing(MBInterface *mbImpl, const char *filename);

MBErrorCode report_nsets(MBInterface *mbImpl);

MBErrorCode report_iface_ents(MBInterface *mbImpl,
                              std::vector<MBParallelComm *> &pcs);

void print_usage(const char *);

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
  MBInterface *mbImpl = new MBCore;
  if (NULL == mbImpl) return 1;
  
  MBErrorCode result = MB_SUCCESS;

    // each interior proc has a vector of N+M vertices, sharing
    // M vertices each with lower- and upper-rank processors, except
    // procs on the end

    // get N, M from command line
  if (argc < 3) {
    if (0 == rank) print_usage(argv[0]);
    err = MPI_Finalize();
    return 1;
  }

  int npos = 1, tag_val, distrib, with_ghosts = 0, resolve_shared = 1, use_mpio = 0;
  const char *tag_name;
  std::vector<std::string> filenames;
  int parallel_option = 0;
  int num_files;

  while (npos != argc) {    
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
                                 resolve_shared, with_ghosts, use_mpio);
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
                                 with_ghosts, use_mpio);
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

  err = MPI_Finalize();

  result = mbImpl->delete_mesh();
  if (MB_SUCCESS != result) {
    std::cerr << "Couldn't delete mesh on rank " << rank
              << "; error message: " << std::endl;
    PRINT_LAST_ERROR;
  }

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

MBErrorCode report_nsets(MBInterface *mbImpl) 
{
    // get and report various numbers...
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
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
                      int with_ghosts, int use_mpio) 
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

  if (1 == use_mpio)
    options << ";USE_MPIO";

  options << ";CPUTIME";

  std::vector<MBEntityHandle> filesets(filenames.size());
  std::vector<MBParallelComm*> pcs(filenames.size());
  std::vector<ReadParallel*> rps(filenames.size());
  MBErrorCode result;

  if (1 < filenames.size()) {
    for (unsigned int i = 0; i < filenames.size(); i++) {
      pcs[i] = new MBParallelComm(mbImpl);
      rps[i] = new ReadParallel(mbImpl, pcs[i]);
    
      result = rps[i]->load_file(filenames[i].c_str(), filesets[i], 
                                 FileOptions(options.str().c_str()), NULL, 0);
      if (MB_SUCCESS != result) 
        PRINT_LAST_ERROR;

      if (MB_SUCCESS != result) {
        MPI_Abort(MPI_COMM_WORLD, result);
        break;
      }

        // exchange tag
      result = pcs[i]->exchange_tags("GLOBAL_ID");
      if (MB_SUCCESS != result) {
        std::cerr << "Tag exchange didn't work." << std::endl;
        break;
      }

    }
  }
  else {
    result = mbImpl->load_file(filenames[0].c_str(), filesets[0], 
                               options.str().c_str());
    pcs[0] = MBParallelComm::get_pcomm(mbImpl, 0);
    assert(pcs[0]);
  }
    
  if (MB_SUCCESS == result) report_iface_ents(mbImpl, pcs);
  
  return result;
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
  std::set<unsigned int> my_addl_procs;
  std::vector<int> addl_procs;
  MBRange nonowned_ghosts[MAX_SHARING_PROCS];
  result = pcomm->pack_buffer(ents, addl_procs, false, true, false, -1,
                              nonowned_ghosts[0], buff, buff_size);
  RRA("Packing buffer count (non-stored handles) failed.");
  nonowned_ghosts[0].clear();

  result = pcomm->unpack_buffer(&buff[0], false, -1, -1, new_ents, nonowned_ghosts,
                                my_addl_procs);
  RRA("Unpacking buffer (non-stored handles) failed.");

  return MB_SUCCESS;
}

MBErrorCode report_iface_ents(MBInterface *mbImpl,
                              std::vector<MBParallelComm *> &pcs) 
{
  MBRange iface_ents[6];
  MBErrorCode result = MB_SUCCESS, tmp_result;
  
    // now figure out which vertices are shared
  MBRange part_ents;
  for (unsigned int p = 0; p < pcs.size(); p++) {
    // get entities owned by this partition
    for (MBRange::iterator rit = pcs[p]->partition_sets().begin();
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
  if (iface_ents[0].size() != iface_ents[5].size())
    std::cerr << "WARNING: number of interface vertices don't agree with "
	      << "vertex adjacencies on interface entities." << std::endl;

  // report # regions owned by this proc
  std::cout << "Proc " << rank << " owns " << part_ents.size() 
	    << " 3d entities." << std::endl;

    // get total # regions over all procs
  int num_local = part_ents.size(), num_total;
  
  int failure = MPI_Reduce(&num_local, &num_total, 1,
                           MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (failure) result = MB_FAILURE;

  if (0 == rank)
    std::cout << "Total # owned regions = " << num_total << std::endl;
    
  return result;
}

