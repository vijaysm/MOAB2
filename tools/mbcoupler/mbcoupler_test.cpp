#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "moab/Core.hpp"
#include "FileOptions.hpp"
#include "ReadParallel.hpp"
#include "Coupler.hpp"
#include "moab_mpi.h"
#include "ElemUtil.hpp"
#include <iostream>
#include <sstream>
#include <assert.h>

using namespace moab;

bool debug = false;

#define RRA(a) if (MB_SUCCESS != result) {\
      std::string tmp_str; mbImpl->get_last_error(tmp_str);\
      tmp_str.append("\n"); tmp_str.append(a);\
      dynamic_cast<Core*>(mbImpl)->get_error_handler()->set_last_error(tmp_str.c_str()); \
      return result;}

#define PRINT_LAST_ERROR \
    if (MB_SUCCESS != result) {\
      std::string tmp_str;\
      std::cout << "Failure; message:" << std::endl;\
      mbImpl->get_last_error(tmp_str);\
      std::cout << tmp_str << std::endl;\
      MPI_Finalize();                                     \
      return result;\
    }

ErrorCode get_file_options(int argc, char **argv, 
                             std::vector<const char *> &filenames,
                             std::string &tag_name,
                             std::string &out_fname,
                             std::string &opts);

ErrorCode report_iface_ents(Interface *mbImpl,
                              std::vector<ParallelComm *> &pcs,
                              bool print_results);

ErrorCode test_interpolation(Interface *mbImpl, 
                               std::string &interp_tag,
                               std::vector<ParallelComm *> &pcs,
                               std::vector<ReadParallel *> rps,
                               double &instant_time,
                               double &pointloc_time,
                               double &interp_time);

int main(int argc, char **argv) 
{
    // need to init MPI first, to tell how many procs and rank
  int err = MPI_Init(&argc, &argv);

  // Declare new mapping functions.
  moab::Element::LinearHex hex;
  moab::Element::LinearTet tet;

  if (argc < 3) {
    std::cerr << "Usage: ";
    std::cerr << argv[0] << " <nfiles> <fname1> ... <fnamen> [interp_tag] [tag_name] [tag_val] [distrib] [with_ghosts]" << std::endl;
    std::cerr << "nfiles        : number of mesh files" << std::endl;
    std::cerr << "fname1..fnamen: mesh files" << std::endl;
    std::cerr << "interp_tag    : name of tag interpolated to target mesh []" << std::endl;
    std::cerr << "output_fname  : name of output file" << std::endl;
    std::cerr << "tag_name      : name of tag used to define partitions [GEOM_DIMENSION]" << std::endl;
    std::cerr << "tag_val       : tag values denoting partition sets [--]" << std::endl;
    std::cerr << "distrib       : if non-zero, distribute the partition sets with tag_val round-robin" << std::endl;
    std::cerr << "with_ghosts   : if non-zero, after initializing in parallel, also exchange one layer of ghost elements" << std::endl;

    err = MPI_Finalize();
    
    return 1;
  }
  
  int nprocs, rank;
  err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // start time
  double stime, rtime, setime, dtime, ltime;
  if (0 == rank) stime = MPI_Wtime();

    // create MOAB instance based on that
  Interface *mbImpl = new Core(rank, nprocs);
  if (NULL == mbImpl) return 1;
  
  ErrorCode result = MB_SUCCESS;

  std::vector<const char *> filenames;
  std::string opts, interp_tag, out_fname;
  result = get_file_options(argc, argv, filenames, interp_tag, out_fname, opts);
  

    // read in mesh(es)
  std::vector<ParallelComm *> pcs(filenames.size()); 
  std::vector<ReadParallel *> rps(filenames.size()); 

  for (unsigned int i = 0; i < filenames.size(); i++) {
    pcs[i] = new ParallelComm(mbImpl);
    rps[i] = new ReadParallel(mbImpl, pcs[i]);
    
    result = rps[i]->load_file(filenames[i], 0, FileOptions(opts.c_str()));
    PRINT_LAST_ERROR;
  }

  result = report_iface_ents(mbImpl, pcs, debug);
  PRINT_LAST_ERROR;

  if (!interp_tag.empty()) {
    double instant_time, pointloc_time, interp_time;
      // test interpolation
    result = test_interpolation(mbImpl, interp_tag, pcs, rps,
                                instant_time, pointloc_time, interp_time);
    PRINT_LAST_ERROR;
    std::cout << "Coupler instantiation/point location/interpolation times = "
              << instant_time << ", " 
              << pointloc_time << ", "
              << interp_time << std::endl;
  }


    // output mesh
  const char *out_option =
//      "PARALLEL_FORMAT";
      NULL;

  if (pcs[1]->proc_config().proc_rank() == 0 && !out_fname.empty()) {
    result = mbImpl->write_file(out_fname.c_str(), NULL, out_option,
                                pcs[1]->partition_sets());
    PRINT_LAST_ERROR;
    std::cout << "Wrote " << out_fname << std::endl;
  }

  std::cout << "Success." << std::endl;

  for (unsigned int i = 0; i < filenames.size(); i++) {
    delete rps[i];
    delete pcs[i];
  }

  delete mbImpl;
  
  err = MPI_Finalize();

  return 0;
}

ErrorCode report_iface_ents(Interface *mbImpl,
                              std::vector<ParallelComm *> &pcs,
                              const bool print_results) 
{
  Range iface_ents[6];
  ErrorCode result = MB_SUCCESS, tmp_result;
  
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
                                   Interface::UNION);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (print_results || iface_ents[0].size() != iface_ents[5].size()) {
    std::cerr << "Proc " << rank << " iface entities: " << std::endl;
    for (int i = 0; i < 4; i++)
      std::cerr << "    " << iface_ents[i].size() << " "
                << i << "d iface entities." << std::endl;
    std::cerr << "    (" << iface_ents[5].size() 
              << " verts adj to other iface ents)" << std::endl;
  }
  
  return result;
}

ErrorCode get_file_options(int argc, char **argv, 
                             std::vector<const char *> &filenames,
                             std::string &interp_tag,
                             std::string &out_fname,
                             std::string &opts) 
{
  int npos = 1;
  int nfiles = atoi(argv[npos++]);
  
    // get mesh filenames
  filenames.resize(nfiles);
  for (int i = 0; i < nfiles; i++) filenames[i] = argv[npos++];

  if (npos < argc) interp_tag = argv[npos++];
  
  if (npos < argc) out_fname = argv[npos++];
  
    // get partition information
  const char *tag_name = "GEOM_DIMENSION";
  int tag_val = -1;
  int distrib = 1;
  int with_ghosts = 0;
  if (npos < argc) tag_name = argv[npos++];
  if (npos < argc) tag_val = strtol(argv[npos++], NULL, 0);
  if (npos < argc) distrib = strtol(argv[npos++], NULL, 0);
  if (npos < argc) with_ghosts = strtol(argv[npos++], NULL, 0);

  if (-1 == tag_val && !strcmp(tag_name, "GEOM_DIMENSION"))
    tag_val = 3;
  
  std::ostringstream options;
  options << "PARALLEL=READ_DELETE;PARTITION=" << tag_name;
  
  if (-1 != tag_val)
    options << ";PARTITION_VAL=" << tag_val;

  if (1 == distrib)
    options << ";PARTITION_DISTRIBUTE";

  options << ";PARALLEL_RESOLVE_SHARED_ENTS";

  if (1 == with_ghosts)
    options << ";PARALLEL_GHOSTS=3.0.1";

  options << ";CPUTIME";
    
  opts = options.str();

  return MB_SUCCESS;
}

ErrorCode test_interpolation(Interface *mbImpl, 
                               std::string &interp_tag,
                               std::vector<ParallelComm *> &pcs,
                               std::vector<ReadParallel *> rps,
                               double &instant_time,
                               double &pointloc_time,
                               double &interp_time) 
{
    // source is 1st mesh, target is 2nd
  Range src_elems, targ_elems;
  ErrorCode result = pcs[0]->get_part_entities(src_elems, 3);
  PRINT_LAST_ERROR;

  double start_time = MPI_Wtime();

    // instantiate a coupler, which also initializes the tree
  Coupler mbc(mbImpl, pcs[0], src_elems, 0);

  instant_time = MPI_Wtime();

    // get points from the target mesh to interpolate
  Range targ_verts, tmp_verts;

    // first get all vertices adj to partition entities in target mesh
  result = pcs[1]->get_part_entities(targ_elems, 3);
  result = mbImpl->get_adjacencies(targ_elems, 0, false, targ_verts, 
                                   Interface::UNION);
  PRINT_LAST_ERROR;

    // then get non-owned verts and subtract 
  result = pcs[1]->get_pstatus_entities(0, PSTATUS_NOT_OWNED, tmp_verts);
  PRINT_LAST_ERROR;
  targ_verts = subtract( targ_verts, tmp_verts);
  
    // get position of these entities; these are the target points
  std::vector<double> vpos(3*targ_verts.size());
  result = mbImpl->get_coords(targ_verts, &vpos[0]);
  PRINT_LAST_ERROR;

    // locate those points in the source mesh
  result = mbc.locate_points(&vpos[0], targ_verts.size());
  PRINT_LAST_ERROR;

  pointloc_time = MPI_Wtime();

    // now interpolate tag onto target points
  std::vector<double> field(targ_verts.size());

  if(interp_tag == "vertex_field"){
    result = mbc.interpolate(Coupler::LINEAR_FE, interp_tag, &field[0]);
  }else if(interp_tag == "element_field"){
    result = mbc.interpolate(Coupler::PLAIN_FE, interp_tag, &field[0]);
  }else{
    std::cout << "Using tag name to determine type of sourge field at the moment... Use either vertex_field or element_field\n";
    result = MB_FAILURE;
  }
  PRINT_LAST_ERROR;
  
  interp_time = MPI_Wtime();

  interp_time -= pointloc_time;
  pointloc_time -= instant_time;
  instant_time -= start_time;

    // set field values as tag on target vertices
  Tag tag;
  result = mbImpl->tag_get_handle(interp_tag.c_str(), tag); PRINT_LAST_ERROR;
  result = mbImpl->tag_set_data(tag, targ_verts, &field[0]); PRINT_LAST_ERROR;

    // done
  return MB_SUCCESS;
}
