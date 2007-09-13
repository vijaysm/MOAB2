#include "ReadParallel.hpp"
#include "MBCore.hpp"
#include "MBProcConfig.hpp"
#include "FileOptions.hpp"
#include "MBError.hpp"
#include "MBReaderWriterSet.hpp"
#include "MBParallelComm.hpp"

#define RR if (MB_SUCCESS != result) return result

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
  const char* parallel_opts[] = { "NONE", "BCAST", "BCAST_DELETE", "SCATTER", 
                                  "FORMAT", 0 };
  enum ParallelOpts {POPT_NONE=0, POPT_BCAST, POPT_BCAST_DELETE, POPT_SCATTER,
                     POPT_FORMAT, POPT_LAST};
      
  MBErrorCode rval = opts.match_option( "PARALLEL", parallel_opts, 
                                        parallel_mode );
  if (MB_FAILURE == rval) {
    merror->set_last_error( "Unexpected value for 'PARALLEL' option\n" );
    return MB_FAILURE;
  }
  else if (MB_ENTITY_NOT_FOUND == rval) {
    parallel_mode = 0;
  }
    // Get partition setting
  std::string partition_tag;
  rval = opts.get_option("PARTITION", partition_tag);
  if (MB_ENTITY_NOT_FOUND == rval || partition_tag.empty())
    partition_tag += "PARTITION";

    // get MPI IO processor rank
  int reader_rank;
  rval = opts.get_int_option( "MPI_IO_RANK", reader_rank );
  if (MB_ENTITY_NOT_FOUND == rval)
    reader_rank = 0;
  else if (MB_SUCCESS != rval) {
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
  
  if (parallel_mode != POPT_SCATTER || 
      reader_rank == (int)(mbImpl->proc_config().rank())) {
      // Try using the file extension to select a reader
    const MBReaderWriterSet* set = impl->reader_writer_set();
    MBReaderIface* reader = set->get_file_extension_reader( file_name );
    if (reader)
    { 
      rval = reader->load_file( file_name, file_set, opts, 
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
          rval = reader->load_file( file_name, file_set, opts, 
                                    material_set_list, num_material_sets );
          delete reader;
          if (MB_SUCCESS == rval)
            break;
        }
      }
    }
  }
  else {
    rval = MB_SUCCESS;
  }
  
  if (parallel_mode == POPT_BCAST ||
      parallel_mode == POPT_BCAST_DELETE) {
    MBRange entities; 
    if (MB_SUCCESS == rval && 
        reader_rank == (int)(mbImpl->proc_config().rank())) {
      rval = mbImpl->get_entities_by_handle( file_set, entities );
      if (MB_SUCCESS != rval)
        entities.clear();
    }
    
    MBParallelComm tool( mbImpl, impl->tag_server(), impl->sequence_manager());
    MBErrorCode tmp_rval = tool.broadcast_entities( reader_rank, entities );
    if (MB_SUCCESS != rval && mbImpl->proc_config().size() != 1)
      tmp_rval = rval;
    else if (MB_SUCCESS != rval) rval = MB_SUCCESS;
      
    if (MB_SUCCESS == rval && 
        reader_rank != (int)(mbImpl->proc_config().rank())) {
      rval = mbImpl->create_meshset( MESHSET_SET, file_set );
      if (MB_SUCCESS == rval) {
        rval = mbImpl->add_entities( file_set, entities );
        if (MB_SUCCESS != rval) {
          mbImpl->delete_entities( &file_set, 1 );
          file_set = 0;
        }
      }
    }

    if (parallel_mode == POPT_BCAST_DELETE)
      rval = delete_nonlocal_entities(partition_tag, file_set);
    
  }
  
  return rval;
}

MBErrorCode ReadParallel::delete_nonlocal_entities(std::string &partition_name,
                                                   MBEntityHandle file_set) 
{
  MBErrorCode result;
  MBError *merror = ((MBCore*)mbImpl)->get_error_handler();
  
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

    // get ents in the partition
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
  
    // gather adjacent verts and add to existing ents
  MBRange tmp_ents;
  MBRange::iterator bit = exist_ents.lower_bound(MBEDGE),
    eit = exist_ents.upper_bound(MBENTITYSET);
  MBRange from_ents(*bit, *eit);
  from_ents = from_ents.intersect(exist_ents);
  result = mbImpl->get_adjacencies(from_ents, 0, false, tmp_ents, 
                                   MBInterface::UNION); RR;
  exist_ents.merge(tmp_ents);
  
    // subtract from all ents to get deletable ents
  all_ents = all_ents.subtract(exist_ents);
  
    // now go through the sets to see if we should keep any
  MBRange all_sets, deletable_sets;
  result = mbImpl->get_entities_by_type(0, MBENTITYSET, all_sets);
  for (MBRange::iterator rit = all_sets.begin(); rit != all_sets.end(); rit++) {
    tmp_ents.clear();
    result = mbImpl->get_entities_by_handle(*rit, tmp_ents); RR;
    tmp_ents = tmp_ents.intersect(exist_ents);
    
      // if the intersection is empty, set is deletable
    if (tmp_ents.empty()) deletable_sets.insert(*rit);
  }
  
    // now delete sets first, then ents
  result = mbImpl->delete_entities(deletable_sets); RR;
  result = mbImpl->delete_entities(all_ents); RR;
  
    // finally, look for sparse tags which have no entities, and delete
    // those too
  std::vector<MBTag> all_tags;
  result = mbImpl->tag_get_tags(all_tags);
  MBTag *tag_vec = &all_tags[0];
  for (unsigned int i = 0; i < all_tags.size(); i++) {
      // get type first, and continue if not sparse
    MBTagType this_type;
    result = mbImpl->tag_get_type(tag_vec[i], this_type); RR;
    if (MB_TAG_SPARSE != this_type) continue;
    
      // get ents with this tag; should be efficient for sparse tags
    tmp_ents.clear();
    result = mbImpl->get_entities_by_type_and_tag(0, MBMAXTYPE, 
                                                  tag_vec+i, NULL,
                                                  1, tmp_ents); RR;
    if (tmp_ents.empty()) {
        // no entities with this tag - delete the tag
      result = mbImpl->tag_delete(tag_vec[i]); RR;
    }
  }
  
  return MB_SUCCESS;
}
