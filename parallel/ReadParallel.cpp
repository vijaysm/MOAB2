#include "ReadParallel.hpp"
#include "MBCore.hpp"
#include "MBProcConfig.hpp"
#include "FileOptions.hpp"
#include "MBError.hpp"
#include "MBReaderWriterSet.hpp"
#include "MBReadUtilIface.hpp"
#include "MBParallelComm.hpp"
#include "MBParallelConventions.h"
#include "MBCN.hpp"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <assert.h>

const bool debug = false;

#define RR(a) if (MB_SUCCESS != result) {                               \
      dynamic_cast<MBCore*>(mbImpl)->get_error_handler()->set_last_error(a); \
      return result;}

enum ParallelActions {PA_READ=0, PA_BROADCAST, PA_DELETE_NONLOCAL,
                      PA_CHECK_GIDS_SERIAL, PA_GET_FILESET_ENTS, 
                      PA_RESOLVE_SHARED_ENTS,
                      PA_EXCHANGE_GHOSTS};
const char *ParallelActionsNames[] = {
    "PARALLEL READ",
    "PARALLEL BROADCAST", 
    "PARALLEL DELETE NONLOCAL",
    "PARALLEL CHECK_GIDS_SERIAL",
    "PARALLEL GET_FILESET_ENTS",
    "PARALLEL RESOLVE_SHARED_ENTS",
    "PARALLEL EXCHANGE_GHOSTS"
};

const char* ReadParallel::parallelOptsNames[] = { "NONE", "BCAST", "BCAST_DELETE", 
                                                  "READ_DELETE", "READ_PARALLEL", 
                                                  "FORMAT", "", 0 };
      
ReadParallel::ReadParallel(MBInterface* impl, 
                           MBParallelComm *pc) 
        : mbImpl(impl), myPcomm(pc) 
{
  if (!myPcomm) {
    myPcomm = MBParallelComm::get_pcomm(mbImpl, 0);
    if (NULL == myPcomm) myPcomm = new MBParallelComm(mbImpl);
  }
}

MBErrorCode ReadParallel::load_file(const char **file_names,
                                    const int num_files,
                                    MBEntityHandle& file_set,
                                    const FileOptions &opts,
                                    const char* set_tag_name,
                                    const int* set_tag_values,
                                    const int num_tag_values ) 
{
  MBError *merror = ((MBCore*)mbImpl)->get_error_handler();

    // Get parallel settings
  int parallel_mode;
  MBErrorCode result = opts.match_option( "PARALLEL", parallelOptsNames, 
                                          parallel_mode );
  if (MB_FAILURE == result) {
    merror->set_last_error( "Unexpected value for 'PARALLEL' option\n" );
    return MB_FAILURE;
  }
  else if (MB_ENTITY_NOT_FOUND == result) {
    parallel_mode = 0;
  }
    // Get partition setting
  std::string partition_tag_name;
  result = opts.get_option("PARTITION", partition_tag_name);
  if (MB_ENTITY_NOT_FOUND == result || partition_tag_name.empty())
    partition_tag_name = PARALLEL_PARTITION_TAG_NAME;

    // Get partition tag value(s), if any, and whether they're to be
    // distributed or assigned
  std::vector<int> partition_tag_vals;
  result = opts.get_ints_option("PARTITION_VAL", partition_tag_vals);
  bool distrib = false;
  result = opts.get_null_option("PARTITION_DISTRIBUTE");
  if (MB_SUCCESS == result) distrib = true;

    // see if we need to report times
  bool cputime = false;
  result = opts.get_null_option("CPUTIME");
  if (MB_SUCCESS == result) cputime = true;

    // get ghosting options
  std::string ghost_str;
  int bridge_dim, ghost_dim = -1, num_layers;
  result = opts.get_str_option("PARALLEL_GHOSTS", ghost_str);
  if (MB_TYPE_OUT_OF_RANGE == result) {
    ghost_dim = 3;
    bridge_dim = 0;
    num_layers = 1;
  }
  else if (MB_SUCCESS == result) {
    int num_fields = 
        sscanf(ghost_str.c_str(), "%d.%d.%d", &ghost_dim, &bridge_dim, &num_layers);
    if (3 != num_fields) {
      merror->set_last_error( "Didn't read 3 fields from PARALLEL_GHOSTS string\n" );
      return MB_FAILURE;
    }
  }

    // get resolve_shared_ents option
  std::string shared_str;
  int resolve_dim = -2, shared_dim = -1;
  result = opts.get_str_option("PARALLEL_RESOLVE_SHARED_ENTS", shared_str);
  if (MB_TYPE_OUT_OF_RANGE == result) {
    resolve_dim = -1;
    shared_dim = -1;
  }
  else if (MB_SUCCESS == result) {
    int num_fields = 
        sscanf(shared_str.c_str(), "%d.%d", &resolve_dim, &shared_dim);
    if (2 != num_fields) {
      merror->set_last_error( "Didn't read 2 fields from PARALLEL_RESOLVE_SHARED_ENTS string\n" );
      return MB_FAILURE;
    }
  }

    // get MPI IO processor rank
  int reader_rank;
  result = opts.get_int_option( "MPI_IO_RANK", reader_rank );
  if (MB_ENTITY_NOT_FOUND == result)
    reader_rank = 0;
  else if (MB_SUCCESS != result) {
    merror->set_last_error( "Unexpected value for 'MPI_IO_RANK' option\n" );
    return MB_FAILURE;
  }

    // now that we've parsed all the parallel options, make an instruction
    // queue
  std::vector<int> pa_vec;
  bool is_reader = (reader_rank == (int) myPcomm->proc_config().proc_rank());
  
  switch (parallel_mode) {
    case POPT_BCAST:
        if (is_reader) {
          pa_vec.push_back(PA_READ);
          pa_vec.push_back(PA_CHECK_GIDS_SERIAL);
          pa_vec.push_back(PA_GET_FILESET_ENTS);
        }
        pa_vec.push_back(PA_BROADCAST);
        if (!is_reader) pa_vec.push_back(PA_GET_FILESET_ENTS);

        break;
    
    case POPT_BCAST_DELETE:
        if (is_reader) {
          pa_vec.push_back(PA_READ);
          pa_vec.push_back(PA_CHECK_GIDS_SERIAL);
          pa_vec.push_back(PA_GET_FILESET_ENTS);
        }
        pa_vec.push_back(PA_BROADCAST);
        if (!is_reader) pa_vec.push_back(PA_GET_FILESET_ENTS);
        pa_vec.push_back(PA_DELETE_NONLOCAL);
        break;

    case POPT_DEFAULT:
    case POPT_READ_DELETE:
        pa_vec.push_back(PA_READ);
        pa_vec.push_back(PA_CHECK_GIDS_SERIAL);
        pa_vec.push_back(PA_GET_FILESET_ENTS);
        pa_vec.push_back(PA_DELETE_NONLOCAL);
        break;

    case POPT_FORMAT:
        merror->set_last_error( "Access to format-specific parallel read not implemented.\n");
        return MB_NOT_IMPLEMENTED;
    case POPT_READ_PARALLEL:
        merror->set_last_error( "Partitioning for PARALLEL=READ_PARALLEL not supported yet.\n");
        return MB_NOT_IMPLEMENTED;
    default:
        return MB_FAILURE;
  }

  if (-2 != resolve_dim) pa_vec.push_back(PA_RESOLVE_SHARED_ENTS);

  if (-1 != ghost_dim) pa_vec.push_back(PA_EXCHANGE_GHOSTS);
  
  
  return load_file(file_names, num_files, file_set, parallel_mode, 
                   partition_tag_name,
                   partition_tag_vals, distrib, pa_vec, opts,
                   set_tag_name, set_tag_values, num_tag_values,
                   reader_rank, cputime, 
                   resolve_dim, shared_dim,
                   ghost_dim, bridge_dim, num_layers);
}
    
MBErrorCode ReadParallel::load_file(const char **file_names,
                                    const int num_files,
                                    MBEntityHandle& file_set,
                                    int parallel_mode, 
                                    std::string &partition_tag_name, 
                                    std::vector<int> &partition_tag_vals, 
                                    bool distrib,
                                    std::vector<int> &pa_vec,
                                    const FileOptions &opts,
                                    const char* set_tag_name,
                                    const int* set_tag_values,
                                    const int num_tag_values,
                                    const int reader_rank,
                                    const bool cputime,
                                    const int resolve_dim,
                                    const int shared_dim,
                                    const int ghost_dim,
                                    const int bridge_dim,
                                    const int num_layers) 
{
  MBErrorCode result = MB_SUCCESS;
  if (myPcomm == NULL)
    myPcomm = new MBParallelComm(mbImpl);

  MBRange entities; 
  MBTag file_set_tag = 0;
  int other_sets = 0;
  MBReaderWriterSet::iterator iter;
  MBRange other_file_sets, file_sets;
  MBCore *impl = dynamic_cast<MBCore*>(mbImpl);

  std::vector<double> act_times(pa_vec.size()+1);
  double stime = 0.0;
  if (cputime) stime = MPI_Wtime();
  std::vector<int>::iterator vit;
  int i, j;
  act_times[0] = MPI_Wtime();
  
    // make a new set for the parallel read
  result = mbImpl->create_meshset(MESHSET_SET, file_set);
  if (MB_SUCCESS != result) return result;
  bool i_read = false;

  for (i = 1, vit = pa_vec.begin();
       vit != pa_vec.end(); vit++, i++) {

    MBErrorCode tmp_result = MB_SUCCESS;
    switch (*vit) {
//==================
      case PA_READ:
          i_read = true;
          
          for (j = 0; j < num_files; j++) {
            if (debug)
              std::cout << "Reading file " << file_names[j] << std::endl;

            MBEntityHandle new_file_set = 0;
            tmp_result = impl->serial_load_file( file_names[j], 
                                                 new_file_set, 
                                                 opts,
                                                 set_tag_name,
                                                 set_tag_values,
                                                 num_tag_values );
            if (MB_SUCCESS != tmp_result) break;

              // put the contents of each file set for the reader into the 
              // file set for the parallel read
            assert(0 != new_file_set);
            MBRange all_ents;
            tmp_result = mbImpl->get_entities_by_handle(new_file_set, all_ents);
            if (MB_SUCCESS != tmp_result) break;
            all_ents.insert(new_file_set);
            tmp_result = mbImpl->add_entities(file_set, all_ents);
            if (MB_SUCCESS != tmp_result) break;
          }
          if (MB_SUCCESS != tmp_result) break;
        
            // mark the file set for this parallel reader
          tmp_result = mbImpl->tag_create("__file_set", sizeof(int), 
                                          MB_TAG_SPARSE, 
                                          MB_TYPE_INTEGER, file_set_tag, 
                                          0, true);
          if (MB_SUCCESS != tmp_result) break;
        
          tmp_result = mbImpl->tag_set_data(file_set_tag, &file_set, 1, 
                                            &other_sets);
          break;

//==================
      case PA_GET_FILESET_ENTS:
          if (debug)
            std::cout << "Getting fileset entities." << std::endl;

            // get entities in the file set, and add actual file set to it;
            // mark the file set to make sure any receiving procs know which it
            // is
          tmp_result = mbImpl->get_entities_by_handle( file_set, entities );
          if (MB_SUCCESS != tmp_result)
            entities.clear();

            // add actual file set to entities too
          entities.insert(file_set);
          break;

//==================
      case PA_BROADCAST:
            // do the actual broadcast; if single-processor, ignore error
          if (debug)
            std::cout << "Broadcasting mesh." << std::endl;

          if (myPcomm->proc_config().proc_size() > 1)
            tmp_result = myPcomm->broadcast_entities( reader_rank, entities );

          if (debug) {
            std::cerr << "Bcast done; entities:" << std::endl;
            mbImpl->list_entities(0, 0);
          }

            // add the received entities to this fileset if I wasn't the reader
          if (!i_read && MB_SUCCESS == tmp_result)
            tmp_result = mbImpl->add_entities(file_set, entities);
          
          break;

//==================
      case PA_DELETE_NONLOCAL:
          if (debug)
            std::cout << "Deleting nonlocal entities." << std::endl;

          tmp_result = delete_nonlocal_entities(partition_tag_name, 
                                                partition_tag_vals, 
                                                distrib,
                                                file_set);
          if (debug) {
            std::cerr << "Delete nonlocal done; entities:" << std::endl;
            mbImpl->list_entities(0, 0);
          }
          break;

//==================
      case PA_CHECK_GIDS_SERIAL:
          if (debug)
            std::cout << "Checking global ids." << std::endl;

          tmp_result = myPcomm->check_global_ids(file_set, 0, 1, true, false);
          break;
        
//==================
      case PA_RESOLVE_SHARED_ENTS:
          if (debug)
            std::cout << "Resolving shared entities." << std::endl;

          tmp_result = myPcomm->resolve_shared_ents(file_set, resolve_dim, shared_dim);
          break;
        
//==================
      case PA_EXCHANGE_GHOSTS:
          if (debug)
            std::cout << "Exchanging ghost entities." << std::endl;

          tmp_result = myPcomm->exchange_ghost_cells(ghost_dim, bridge_dim, 
                                                     num_layers, true);
          break;
        
//==================
      default:
          return MB_FAILURE;
    }

    if (MB_SUCCESS != tmp_result) {
      result = tmp_result;
      if (myPcomm->proc_config().proc_size() != 1) {
        std::ostringstream ostr;
        ostr << "Failed in step " << ParallelActionsNames[*vit] << std::endl;
        std::string tmp_str;
        if (MB_SUCCESS == mbImpl->get_last_error(tmp_str)) ostr << tmp_str << std::endl;
        RR(ostr.str().c_str());
      }
      break;
    }

    if (cputime) act_times[i] = MPI_Wtime();
  }

  if (cputime && 0 == myPcomm->proc_config().proc_rank()) {
    std::cout << "Read times: ";
    for (i = 1, vit = pa_vec.begin();
         vit != pa_vec.end(); vit++, i++) 
      std::cout << act_times[i] - act_times[i-1] << " ";
    std::cout << "(";
    for (vit = pa_vec.begin();
         vit != pa_vec.end(); vit++) 
      std::cout << ParallelActionsNames[*vit] << "/";
    std::cout << ")" << std::endl;
  }
  
  return result;
}

MBErrorCode ReadParallel::delete_nonlocal_entities(std::string &ptag_name,
                                                   std::vector<int> &ptag_vals,
                                                   bool distribute,
                                                   MBEntityHandle file_set) 
{
  MBRange partition_sets;
  MBErrorCode result;

  MBTag ptag;
  result = mbImpl->tag_get_handle(ptag_name.c_str(), ptag); 
  RR("Failed getting tag handle in delete_nonlocal_entities.");

  result = mbImpl->get_entities_by_type_and_tag(file_set, MBENTITYSET,
                                                &ptag, NULL, 1,
                                                myPcomm->partition_sets());
  RR("Failed to get sets with partition-type tag.");

  int proc_sz = myPcomm->proc_config().proc_size();
  int proc_rk = myPcomm->proc_config().proc_rank();

  if (!ptag_vals.empty()) {
      // values input, get sets with those values
    MBRange tmp_sets;
    std::vector<int> tag_vals(myPcomm->partition_sets().size());
    result = mbImpl->tag_get_data(ptag, myPcomm->partition_sets(), &tag_vals[0]);
    RR("Failed to get tag data for partition vals tag.");
    for (std::vector<int>::iterator pit = tag_vals.begin(); 
         pit != tag_vals.end(); pit++) {
      std::vector<int>::iterator pit2 = std::find(ptag_vals.begin(),
                                                  ptag_vals.end(), *pit);
      if (pit2 != ptag_vals.end()) 
        tmp_sets.insert(myPcomm->partition_sets()[pit - tag_vals.begin()]);
    }

    myPcomm->partition_sets().swap(tmp_sets);
  }

  if (distribute) {
      // for now, require that number of partition sets be greater
      // than number of procs
    if (myPcomm->partition_sets().size() < (unsigned int) proc_sz) {
      result = MB_FAILURE;
      RR("Number of procs greater than number of partitions.");
    }
    
    MBRange tmp_sets;
      // distribute the partition sets
    unsigned int num_sets = myPcomm->partition_sets().size() / proc_sz;
    unsigned int num_leftover = myPcomm->partition_sets().size() % proc_sz;
    int begin_set = 0;
    if (proc_rk < (int) num_leftover) {
      num_sets++;
      begin_set = num_sets * proc_rk;
    }
    else
      begin_set = proc_rk * num_sets + num_leftover;
      

    for (unsigned int i = 0; i < num_sets; i++)
      tmp_sets.insert(myPcomm->partition_sets()[begin_set+i]);
    
    myPcomm->partition_sets().swap(tmp_sets);
  }

  if (debug) {
    std::cerr << "My partition sets: ";
    myPcomm->partition_sets().print();
  }
  
  result = delete_nonlocal_entities(file_set); RR(" ");
  
  if (ptag_name != PARALLEL_PARTITION_TAG_NAME) {
      // tag the partition sets with a standard tag name
    result = mbImpl->tag_create(PARALLEL_PARTITION_TAG_NAME, sizeof(int), 
                                MB_TAG_SPARSE, 
                                MB_TYPE_INTEGER, ptag, 
                                0, true);
    if (MB_ALREADY_ALLOCATED == result) {
        // this tag already exists; better check to see that tagged sets
        // agree with this partition
      MBRange tagged_sets;
      int *proc_rk_ptr = &proc_rk;
      result = mbImpl->get_entities_by_type_and_tag(file_set, MBENTITYSET, &ptag, 
                                                    (const void* const*)&proc_rk_ptr, 1,
                                                    tagged_sets); RR(" ");
      if (!tagged_sets.empty() && tagged_sets != myPcomm->partition_sets()) {
        result = mbImpl->tag_delete_data(ptag, tagged_sets); RR(" ");
      }
      else if (tagged_sets == myPcomm->partition_sets()) return MB_SUCCESS;
    }

      // if we get here, we need to assign the tag
    std::vector<int> values(myPcomm->partition_sets().size());
    for (unsigned int i = 0; i < myPcomm->partition_sets().size(); i++)
      values[i] = proc_rk;
    result = mbImpl->tag_set_data(ptag, myPcomm->partition_sets(), &values[0]); RR(" ");
  }

  return result;
}

MBErrorCode ReadParallel::delete_nonlocal_entities(MBEntityHandle file_set) 
{

  MBErrorCode result;

    // get partition entities and ents related to/used by those
    // get ents in the partition
  std::string iface_name = "MBReadUtilIface";
  MBReadUtilIface *read_iface;
  mbImpl->query_interface(iface_name, reinterpret_cast<void**>(&read_iface));
  MBRange partition_ents, all_sets;

  if (debug) std::cout << "Gathering related entities." << std::endl;
  
  result = read_iface->gather_related_ents(myPcomm->partition_sets(), partition_ents,
                                           &all_sets);
  RR("Failure gathering related entities.");

    // get pre-existing entities
  MBRange file_ents;
  result = mbImpl->get_entities_by_handle(file_set, file_ents); 
  RR("Couldn't get pre-existing entities.");

  if (debug && myPcomm->proc_config().proc_rank() == 0) {
    std::cout << "File entities: " << std::endl;
    file_ents.print("ff  ");
  }
  
    // get deletable entities by subtracting partition ents from file ents
  MBRange deletable_ents = file_ents.subtract(partition_ents);

    // cache deletable vs. keepable sets
  MBRange deletable_sets = all_sets.intersect(deletable_ents);
  MBRange keepable_sets = all_sets.subtract(deletable_sets);
  
  if (debug) std::cout << "Removing deletable entities from keepable sets." << std::endl;

    // remove deletable ents from all keepable sets
  for (MBRange::iterator rit = keepable_sets.begin();
       rit != keepable_sets.end(); rit++) {
    result = mbImpl->remove_entities(*rit, deletable_ents);
    RR("Failure removing deletable entities.");
  }

  if (debug) {
    std::cout << "Deleting deletable entities." << std::endl;

    if (myPcomm->proc_config().proc_rank() == 0) {
      std::cout << "Deletable sets: " << std::endl;
      deletable_sets.print("ff  ");
    }
  }
  
    // delete sets, then ents
  if (!deletable_sets.empty())
    result = mbImpl->delete_entities(deletable_sets);
  RR("Failure deleting sets in delete_nonlocal_entities.");

  deletable_ents = deletable_ents.subtract(deletable_sets);

  if (debug && myPcomm->proc_config().proc_rank() == 0) {
    std::cout << "Deletable entities: " << std::endl;
    deletable_ents.print("ff  ");
  }
  
  if (!deletable_ents.empty())
    result = mbImpl->delete_entities(deletable_ents);
  RR("Failure deleting entities in delete_nonlocal_entities.");

  return result;
}
