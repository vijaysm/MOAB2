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

const bool debug = false;

#define RR(a) if (MB_SUCCESS != result) {\
          dynamic_cast<MBCore*>(mbImpl)->get_error_handler()->set_last_error(a);\
          return result;}

enum ParallelActions {PA_READ=0, PA_BROADCAST, PA_DELETE_NONLOCAL,
                      PA_CHECK_GIDS_SERIAL, PA_GET_FILESET_ENTS};
const char *ParallelActionsNames[] = {
  "PARALLEL READ",
  "PARALLEL BROADCAST", 
  "PARALLEL DELETE NONLOCAL",
  "PARALLEL CHECK_GIDS_SERIAL",
  "PARALLEL GET_FILESET_ENTS"
};

MBErrorCode ReadParallel::load_file(const char *file_name,
                                    MBEntityHandle& file_set,
                                    const FileOptions &opts,
                                    const int* material_set_list,
                                    const int num_material_sets ) 
{
  MBError *merror = ((MBCore*)mbImpl)->get_error_handler();

    // Get parallel settings
  int parallel_mode;
  const char* parallel_opts[] = { "NONE", "BCAST", "BCAST_DELETE", 
                                  "READ_DELETE", "READ_PARALLEL", 
                                  "FORMAT", 0 };
  enum ParallelOpts {POPT_NONE=0, POPT_BCAST, POPT_BCAST_DELETE, 
                     POPT_READ_DELETE, POPT_READ_PARALLEL,
                     POPT_FORMAT, POPT_LAST};
      
  MBErrorCode result = opts.match_option( "PARALLEL", parallel_opts, 
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
    partition_tag_name += "PARTITION";

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
  bool is_reader = (mbImpl->proc_rank() == reader_rank);
  
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
  
  return load_file(file_name, file_set, parallel_mode, partition_tag_name,
                   partition_tag_vals, distrib, pa_vec, material_set_list,
                   num_material_sets, opts, reader_rank, cputime);
}
    
MBErrorCode ReadParallel::load_file(const char *file_name,
                                    MBEntityHandle& file_set,
                                    int parallel_mode, 
                                    std::string &partition_tag_name, 
                                    std::vector<int> &partition_tag_vals, 
                                    bool distrib,
                                    std::vector<int> &pa_vec,
                                    const int* material_set_list,
                                    const int num_material_sets,
                                    const FileOptions &opts,
                                    int reader_rank,
                                    bool cputime) 
{
  MBErrorCode result = MB_SUCCESS;
  MBParallelComm pcom( mbImpl);
  MBRange entities; 
  MBTag file_set_tag = 0;
  int other_sets;
  MBReaderIface* reader;
  MBReaderWriterSet::iterator iter;
  MBRange other_file_sets, file_sets;
  int tag_val, *tag_val_ptr = &tag_val;
  MBCore *impl = dynamic_cast<MBCore*>(mbImpl);

  double act_times[10] = {0.0};
  double stime = 0.0;
  if (cputime) stime = MPI_Wtime();

  for (std::vector<int>::iterator vit = pa_vec.begin();
       vit != pa_vec.end(); vit++) {

    MBErrorCode tmp_result = MB_SUCCESS;
    switch (*vit) {
//==================
      case PA_READ:
        if (debug)
          std::cout << "Reading file " << file_name << std::endl;
            
        reader = impl->reader_writer_set()->
          get_file_extension_reader( file_name );
        if (reader)
        { 
          tmp_result = reader->load_file( file_name, file_set, opts, 
                                          material_set_list, num_material_sets );
          delete reader;
        }
        else
        {  
            // Try all the readers
          for (iter = impl->reader_writer_set()->begin(); 
               iter != impl->reader_writer_set()->end(); ++iter)
          {
            reader = iter->make_reader( mbImpl );
            if (NULL != reader)
            {
              tmp_result = reader->load_file( file_name, file_set, opts, 
                                              material_set_list, num_material_sets );
              delete reader;
              if (MB_SUCCESS == tmp_result)
                break;
            }
          }
        }
        if (MB_SUCCESS != tmp_result) break;

          // mark the file set
        other_sets = 0;
        tmp_result = mbImpl->tag_create("__file_set", sizeof(int), 
                                        MB_TAG_SPARSE, 
                                        MB_TYPE_INTEGER, file_set_tag, 
                                        0, true);
        if (MB_ALREADY_ALLOCATED == tmp_result) {
          tmp_result = mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET,
                                                            &file_set_tag, NULL, 1,
                                                            other_file_sets);
          if (MB_SUCCESS == tmp_result) other_sets = other_file_sets.size();
        }
        if (MB_SUCCESS == tmp_result)
          tmp_result = mbImpl->tag_set_data(file_set_tag, &file_set, 1, 
                                            &other_sets);
        break;

//==================
      case PA_GET_FILESET_ENTS:
        if (debug)
          std::cout << "Getting fileset entities." << std::endl;

        if (0 == file_set_tag) {
          tmp_result = mbImpl->tag_get_handle("FILE_SET", file_set_tag);
          if (MB_SUCCESS == tmp_result) {
            other_file_sets.clear();
            file_sets.clear();
            tmp_result = mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET,
                                                              &file_set_tag, 
                                                              NULL, 1,
                                                              other_file_sets);
            if (MB_SUCCESS == tmp_result && other_file_sets.size() > 1) {
              tag_val = other_file_sets.size();
              result = mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET,
                                                            &file_set_tag, 
                                                            (void*const*) 
                                                            &tag_val_ptr, 1,
                                                            file_sets);
              if (!file_sets.empty()) other_file_sets = file_sets;
            }
            if (!other_file_sets.empty()) file_set = *other_file_sets.rbegin();
          }
        }
        
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

        if (mbImpl->proc_size() > 1)
          tmp_result = pcom.broadcast_entities( reader_rank, entities );

        if (debug) {
          std::cerr << "Bcast done; entities:" << std::endl;
          mbImpl->list_entities(0, 0);
        }
        break;

//==================
      case PA_DELETE_NONLOCAL:
        if (debug)
          std::cout << "Deleting nonlocal entities." << std::endl;

        tmp_result = delete_nonlocal_entities(partition_tag_name, 
                                              partition_tag_vals, 
                                              distrib,
                                              file_set);
        break;

//==================
      case PA_CHECK_GIDS_SERIAL:
        if (debug)
          std::cout << "Checking global ids." << std::endl;

        tmp_result = pcom.check_global_ids(file_set, 0, 1, true, false);
        break;
        
//==================
      default:
        return MB_FAILURE;
    }

    if (MB_SUCCESS != tmp_result &&
        (*vit != PA_BROADCAST || mbImpl->proc_size() != 1)) {
      result = tmp_result;
      std::ostringstream ostr;
      ostr << "Failed in step " << ParallelActionsNames[*vit] << std::endl;
      std::string tmp_str;
      if (MB_SUCCESS == mbImpl->get_last_error(tmp_str)) ostr << tmp_str << std::endl;
      RR(ostr.str().c_str());
    }

    if (cputime) act_times[*vit] = MPI_Wtime() - 
                     (*vit ? act_times[*vit - 1] : stime);
  }

  if (cputime && 0 == mbImpl->proc_rank()) {
    std::cout << "Read times: ";
    for (std::vector<int>::iterator vit = pa_vec.begin();
         vit != pa_vec.end(); vit++) 
      std::cout << act_times[*vit] << " ";
    std::cout << "(";
    for (std::vector<int>::iterator vit = pa_vec.begin();
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
                                                partition_sets);
  RR("Failed to get sets with partition-type tag.");

  int proc_sz = mbImpl->proc_size();
  int proc_rk = mbImpl->proc_rank();
  unsigned int num_partsets = partition_sets.size();

  if (!ptag_vals.empty()) {
      // values input, get sets with those values
    MBRange tmp_sets;
    std::vector<int> tag_vals(partition_sets.size());
    result = mbImpl->tag_get_data(ptag, partition_sets, &tag_vals[0]);
    RR("Failed to get tag data for partition vals tag.");
    for (std::vector<int>::iterator pit = tag_vals.begin(); 
         pit != tag_vals.end(); pit++) {
      std::vector<int>::iterator pit2 = std::find(ptag_vals.begin(),
                                                  ptag_vals.end(), *pit);
      if (pit2 != ptag_vals.end()) 
        tmp_sets.insert(partition_sets[pit - tag_vals.begin()]);
    }

    partition_sets.swap(tmp_sets);
  }

  if (distribute) {
    MBRange tmp_sets;
      // distribute the partition sets
    unsigned int num_sets = partition_sets.size() / proc_sz;
    if (proc_rk < (int) (partition_sets.size() % proc_sz)) num_sets++;

      // cut them in half if we're on one proc
    if (proc_sz == 1 && num_partsets == num_sets) num_sets /= 2;
    
    for (unsigned int i = 0; i < num_sets; i++) 
      tmp_sets.insert(partition_sets[i*proc_sz + proc_rk]);

    partition_sets.swap(tmp_sets);
  }

  if (debug) {
    std::cerr << "My partition sets: ";
    partition_sets.print();
  }
  
  return delete_nonlocal_entities(partition_sets, file_set);
}

MBErrorCode ReadParallel::delete_nonlocal_entities(MBRange &partition_sets,
                                                   MBEntityHandle file_set) 
{

  MBErrorCode result;

    // get partition entities and ents related to/used by those
    // get ents in the partition
  std::string iface_name = "MBReadUtilIface";
  MBReadUtilIface *read_iface;
  mbImpl->query_interface(iface_name, reinterpret_cast<void**>(&read_iface));
  MBRange partition_ents, all_sets;

  if (debug) std::cout << "Gathering related entities." << std::endl;
  
  result = read_iface->gather_related_ents(partition_sets, partition_ents,
                                           &all_sets);
  RR("Failure gathering related entities.");

    // get pre-existing entities
  MBRange file_ents;
  result = mbImpl->get_entities_by_handle(file_set, file_ents); 
  RR("Couldn't get pre-existing entities.");

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

  if (debug) std::cout << "Deleting deletable entities." << std::endl;
  
    // delete sets, then ents
  if (!deletable_sets.empty())
    result = mbImpl->delete_entities(deletable_sets);
  RR("Failure deleting sets in delete_nonlocal_entities.");

  deletable_ents = deletable_ents.subtract(deletable_sets);
  if (!deletable_ents.empty())
    result = mbImpl->delete_entities(deletable_ents);
  RR("Failure deleting entities in delete_nonlocal_entities.");

    // mark partition sets with partition tag, needed later for
    // establishing interface sets
  MBTag partition_set_tag;
  result = mbImpl->tag_create(PARALLEL_PARTITION_TAG_NAME, sizeof(int),
                              MB_TAG_SPARSE, MB_TYPE_INTEGER, 
                              partition_set_tag, NULL, true);
  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) {
    RR("Couldn't create/get partition set tag.");
  }
  std::vector<int> pset_vals(partition_sets.size());
  std::fill(pset_vals.begin(), pset_vals.end(), mbImpl->proc_rank());
  result = mbImpl->tag_set_data(partition_set_tag, partition_sets, 
                                &pset_vals[0]);
  RR("Couldn't set partition set tag value.");

  return result;
}
