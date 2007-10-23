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
                                  "READ_DELETE", "SCATTER", 
                                  "FORMAT", 0 };
  enum ParallelOpts {POPT_NONE=0, POPT_BCAST, POPT_BCAST_DELETE, 
                     POPT_READ_DELETE, POPT_SCATTER,
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
    case POPT_SCATTER:
      merror->set_last_error( "Partitioning for PARALLEL=SCATTER not supported yet.\n");
      return MB_NOT_IMPLEMENTED;
    default:
      return MB_FAILURE;
  }
  
  return load_file(file_name, file_set, parallel_mode, partition_tag_name,
                   partition_tag_vals, distrib, pa_vec, material_set_list,
                   num_material_sets, opts, reader_rank);
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
                                    int reader_rank) 
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
  

  for (std::vector<int>::iterator vit = pa_vec.begin();
       vit != pa_vec.end(); vit++) {

    MBErrorCode tmp_result;
    switch (*vit) {
//==================
      case PA_READ:
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
        tmp_result = pcom.broadcast_entities( reader_rank, entities );
        if (debug) {
          std::cerr << "Bcast done; entities:" << std::endl;
          mbImpl->list_entities(0, 0);
        }
        break;

//==================
      case PA_DELETE_NONLOCAL:
        tmp_result = delete_nonlocal_entities(partition_tag_name, 
                                              partition_tag_vals, 
                                              distrib,
                                              file_set);
        break;

//==================
      case PA_CHECK_GIDS_SERIAL:
        tmp_result = pcom.check_global_ids(file_set, 0, 1, true, false);
        break;
        
//==================
      default:
        return MB_FAILURE;
    }

    if (MB_SUCCESS != tmp_result) {
      result = tmp_result;
      std::ostringstream ostr;
      ostr << "Failed in step " << ParallelActionsNames[*vit] << std::endl;
      RR(ostr.str().c_str());
    }
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
    int num_sets = partition_sets.size() / proc_sz, orig_numsets = num_sets;
    if (proc_rk < partition_sets.size() % proc_sz) num_sets++;

    for (int i = 0; i < num_sets; i++) 
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

  if (debug) std::cerr << "Deleting non-local entities." << std::endl;
  
  MBErrorCode result;
  MBError *merror = ((MBCore*)mbImpl)->get_error_handler();

    // get partition entities and ents related to/used by those
    // get ents in the partition
  std::string iface_name = "MBReadUtilIface";
  MBReadUtilIface *read_iface;
  mbImpl->query_interface(iface_name, reinterpret_cast<void**>(&read_iface));
  MBRange partition_ents, all_sets;
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
  
    // remove deletable ents from all keepable sets
  for (MBRange::iterator rit = keepable_sets.begin();
       rit != keepable_sets.end(); rit++) {
    result = mbImpl->remove_entities(*rit, deletable_ents);
    RR("Failure removing deletable entities.");
  }

    // delete sets, then ents
  result = mbImpl->delete_entities(deletable_sets);
  RR("Failure deleting sets in delete_nonlocal_entities.");

  deletable_ents = deletable_ents.subtract(deletable_sets);
  result = mbImpl->delete_entities(deletable_ents);
  RR("Failure deleting entities in delete_nonlocal_entities.");

//  if (debug)
//    result = ((MBCore*)mbImpl)->check_adjacencies();

  return result;

/*  


// ================================  
    // get entities in this partition
  int my_rank = (int)mbImpl->proc_config().rank();
  if (my_rank == 0 && mbImpl->proc_config().size() == 1) my_rank = 1;
  int *my_rank_ptr = &my_rank;
  MBTag partition_tag;
  
  result = mbImpl->tag_get_handle(partition_name.c_str(), partition_tag);
  if (MB_TAG_NOT_FOUND == result) {
    merror->set_last_error( "Couldn't find partition tag\n");
    return result;
  }
  else if (MB_SUCCESS != result) return result;
    
  MBRange partition_sets;
  result = mbImpl->get_entities_by_type_and_tag(file_set, MBENTITYSET,
                                                &partition_tag, 
                                                (const void* const *) &my_rank_ptr, 
                                                1, partition_sets); RR;
  if (MB_SUCCESS != result || partition_sets.empty()) return result;
  
  MBRange file_ents, partition_ents, exist_ents, all_ents;

  for (MBRange::iterator rit = partition_sets.begin(); 
       rit != partition_sets.end(); rit++) {
    result = mbImpl->get_entities_by_handle(*rit, partition_ents, 
                                            MBInterface::UNION); RR;
  }

    // get pre-existing ents, which are all entities minus file ents
  result = mbImpl->get_entities_by_handle(0, all_ents); RR;
  result = mbImpl->get_entities_by_handle(file_set, file_ents); RR;
  exist_ents = all_ents.subtract(file_ents);

    // merge partition ents into pre-existing entities
  exist_ents.merge(partition_ents);
  
    // gather adjacent ents of lower dimension and add to existing ents
  MBRange tmp_ents;
  for (int dim = 2; dim >= 0; dim--) {
    MBEntityType lower_type = MBCN::TypeDimensionMap[dim+1].first,
      upper_type = MBCN::TypeDimensionMap[3].second;
    
    MBRange::const_iterator bit = exist_ents.lower_bound(lower_type),
      eit = exist_ents.upper_bound(upper_type);
    MBRange from_ents;
    from_ents.merge(bit, eit);
    tmp_ents.clear();
    result = mbImpl->get_adjacencies(from_ents, dim, false, tmp_ents, 
                                     MBInterface::UNION); RR;
    exist_ents.merge(tmp_ents);
  }
  
    // subtract from all ents to get deletable ents
  all_ents = all_ents.subtract(exist_ents);
  
    // go through the sets to which ones we should keep
  MBRange all_sets, deletable_sets;
  result = mbImpl->get_entities_by_type(0, MBENTITYSET, all_sets);
  for (MBRange::iterator rit = all_sets.begin(); rit != all_sets.end(); rit++) {
    tmp_ents.clear();
    result = mbImpl->get_entities_by_handle(*rit, tmp_ents, true); RR;
    MBRange tmp_ents2 = tmp_ents.intersect(exist_ents);
    
      // if the intersection is empty, set is deletable
    if (tmp_ents2.empty()) deletable_sets.insert(*rit);
    
    else if (tmp_ents.size() > tmp_ents2.size()) {
        // more elements in set or contained sets than we're keeping; delete 
        // the difference from just this set, to remove entities to be deleted below
        // it's ok if entity isn't contained, doesn't generate an error
      tmp_ents = tmp_ents.subtract(tmp_ents2);
      result = mbImpl->remove_entities(*rit, tmp_ents); RR;
    }
  }

    // take the deletable sets out of other sets so we don't end up
    // with stale set handles
  for (MBRange::iterator rit = all_sets.begin(); rit != all_sets.end(); rit++) {
    if (deletable_sets.find(*rit) == deletable_sets.end()) {
      result = mbImpl->remove_entities(*rit, deletable_sets); RR;
    }
  }

    // remove sets from all_ents, since they're dealt with separately
  all_ents = all_ents.subtract(all_sets);
  
    // now delete sets first, then ents
  result = mbImpl->delete_entities(deletable_sets); RR;
  result = mbImpl->delete_entities(all_ents); RR;
  
  result = ((MBCore*)mbImpl)->check_adjacencies();
  
  return result;

*/
}
