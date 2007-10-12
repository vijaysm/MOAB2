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

const bool debug = true;

#define RR(a) if (MB_SUCCESS != result) {\
          dynamic_cast<MBCore*>(mbImpl)->get_error_handler()->set_last_error(a);\
          return result;}

MBErrorCode ReadParallel::load_file(const char *file_name,
                                    MBEntityHandle& file_set,
                                    const FileOptions &opts,
                                    const int* material_set_list,
                                    const int num_material_sets ) 
{
  MBError *merror = ((MBCore*)mbImpl)->get_error_handler();

  MBCore *impl = dynamic_cast<MBCore*>(mbImpl);
  
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

    // Get partition tag value(s), if any
  std::vector<int> partition_tag_vals;
  result = opts.get_ints_option("PARTITION_VAL", partition_tag_vals);

    // get MPI IO processor rank
  int reader_rank;
  result = opts.get_int_option( "MPI_IO_RANK", reader_rank );
  if (MB_ENTITY_NOT_FOUND == result)
    reader_rank = 0;
  else if (MB_SUCCESS != result) {
    merror->set_last_error( "Unexpected value for 'MPI_IO_RANK' option\n" );
    return MB_FAILURE;
  }
  
    // now that we've parsed all the parallel options, return
    // failure for most of them because we haven't implemented 
    // most of them yet.
  if (parallel_mode == POPT_FORMAT) {
    merror->set_last_error( "Access to format-specific parallel read not implemented.\n");
    return MB_NOT_IMPLEMENTED;
  }

  if (parallel_mode == POPT_SCATTER) {
    merror->set_last_error( "Partitioning for PARALLEL=SCATTER not supported yet.\n");
    return MB_NOT_IMPLEMENTED;
  }

  if ((parallel_mode != POPT_SCATTER && 
       parallel_mode != POPT_BCAST_DELETE) || 
      reader_rank == (int)(mbImpl->proc_rank())) {
      // Try using the file extension to select a reader
    const MBReaderWriterSet* set = impl->reader_writer_set();
    MBReaderIface* reader = set->get_file_extension_reader( file_name );
    if (reader)
    { 
      result = reader->load_file( file_name, file_set, opts, 
                                material_set_list, num_material_sets );
      delete reader;
    }
    else
    {  
        // Try all the readers
      MBReaderWriterSet::iterator iter;
      for (iter = set->begin(); iter != set->end(); ++iter)
      {
        MBReaderIface* reader = iter->make_reader( mbImpl );
        if (NULL != reader)
        {
          result = reader->load_file( file_name, file_set, opts, 
                                    material_set_list, num_material_sets );
          delete reader;
          if (MB_SUCCESS == result)
            break;
        }
      }
    }
  }
  else {
    result = MB_SUCCESS;
  }
  
  if (MB_SUCCESS != result)
    RR("Failed initial file load.");
    
  if (parallel_mode == POPT_BCAST ||
      parallel_mode == POPT_BCAST_DELETE) {
    MBRange entities; 
    MBParallelComm pcom( mbImpl);

      // get which entities need to be broadcast, only if I'm the reader
    if (reader_rank == (int)(mbImpl->proc_rank())) {

        // if I'm root, check to make sure we have global ids (at least for
        // vertices, anyway) & generate (in serial) if not
      result = pcom.check_global_ids(file_set, 0, 1, true, false);
      RR("Failed to generate/find global ids for parallel read.");
      
      result = mbImpl->get_entities_by_handle( file_set, entities );
      if (MB_SUCCESS != result)
        entities.clear();

        // add actual file set to entities too
      entities.insert(file_set);

        // mark the file set so the receivers know which one it is
      MBTag file_set_tag;
      int other_sets = 0;
      result = mbImpl->tag_create("FILE_SET", sizeof(int), MB_TAG_SPARSE, 
                                  MB_TYPE_INTEGER, file_set_tag, 0, true);
      if (MB_ALREADY_ALLOCATED == result) {
        MBRange other_file_sets;
        result = mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET,
                                                      &file_set_tag, NULL, 1,
                                                      other_file_sets);
        if (MB_SUCCESS == result) other_sets = other_file_sets.size();
      }
      if (MB_SUCCESS == result)
        result = mbImpl->tag_set_data(file_set_tag, &file_set, 1, &other_sets);
    }

      // do the actual broadcast; if single-processor, ignore error
    result = pcom.broadcast_entities( reader_rank, entities );
    if (mbImpl->proc_size() == 1 && MB_SUCCESS != result) 
      result = MB_SUCCESS;
      
      // go get the file set if I'm not the reader
    if (MB_SUCCESS == result && 
        reader_rank != (int)(mbImpl->proc_rank())) {
      MBTag file_set_tag;
      result = mbImpl->tag_get_handle("FILE_SET", file_set_tag);
      if (MB_SUCCESS == result) {
        MBRange other_file_sets;
        result = mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET,
                                                      &file_set_tag, NULL, 1,
                                                      other_file_sets);
        if (MB_SUCCESS == result && other_file_sets.size() > 1) {
          int tag_val = other_file_sets.size(), *tag_val_ptr = &tag_val;
          MBRange file_sets;
          result = mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET,
                                                        &file_set_tag, 
                                                        (void*const*) &tag_val_ptr, 1,
                                                        file_sets);
          if (!file_sets.empty()) other_file_sets = file_sets;
        }
        if (!other_file_sets.empty()) file_set = *other_file_sets.rbegin();
      }
    }
    
    RR("Failed to broadcast mesh.");
  }

  if (parallel_mode == POPT_BCAST_DELETE ||
      parallel_mode == POPT_READ_DELETE) {
    result = delete_nonlocal_entities(partition_tag_name, 
                                      partition_tag_vals, 
                                      file_set);
    if (MB_SUCCESS != result) return result;
  }
  
  return result;
}

MBErrorCode ReadParallel::delete_nonlocal_entities(std::string &ptag_name,
                                                   std::vector<int> &ptag_vals,
                                                   MBEntityHandle file_set) 
{
  MBRange partition_sets, my_sets;
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

  if (ptag_vals.empty()) {
      // no values input, just distribute sets
    int num_sets = partition_sets.size() / proc_sz, orig_numsets = num_sets;
    if (partition_sets.size() % proc_sz != 0) {
      num_sets++;
      if (proc_rk == proc_sz-1) 
        num_sets = partition_sets.size() % num_sets;
    }

    int istart = orig_numsets * proc_rk;
    for (int i = 0; i < num_sets; i++) 
      my_sets.insert(partition_sets[istart+i]);
  }
  else {
      // values input, get sets with those values
    std::vector<int> tag_vals(partition_sets.size());
    result = mbImpl->tag_get_data(ptag, partition_sets, &tag_vals[0]);
    RR("Failed to get tag data for partition vals tag.");
    for (std::vector<int>::iterator pit = ptag_vals.begin(); 
         pit != ptag_vals.end(); pit++) {
      std::vector<int>::iterator pit2 = std::find(tag_vals.begin(),
                                                  tag_vals.end(), *pit);
      if (pit2 == tag_vals.end()) RR("Couldn't find partition tag value.");
      my_sets.insert(partition_sets[pit2 - tag_vals.begin()]);
    }
  }
  
  return delete_nonlocal_entities(my_sets, file_set);
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
