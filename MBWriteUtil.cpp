#include "MBWriteUtil.hpp"
#include "MBCore.hpp"
#include "MBError.hpp"
#include "EntitySequenceManager.hpp"
#include "TagServer.hpp"

#ifdef WIN32
#ifdef _DEBUG
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif
#endif

MBWriteUtil::MBWriteUtil(MBCore* mdb, MBError* error_handler) 
    : MBWriteUtilIface(), mMB(mdb), mError(error_handler)
{
}


MBErrorCode MBWriteUtil::get_node_arrays(
    const int num_arrays,
    const int num_nodes, 
    const MBRange& entities, 
    MBTag node_id_tag,
    std::vector<double*>& arrays)
{
  MBErrorCode error = MB_SUCCESS;

  TagServer* tag_server = mMB->tag_server();

  // check the data coming into the function
  // dimension should be proper
  if(num_arrays < 1 || num_arrays > 3)
    return MB_FAILURE;

  // number of nodes should be greater than zero
  if(num_nodes < 1)
    return MB_FAILURE;
  
  // there should be some entities
  if(entities.empty())
    return MB_FAILURE;

  // memory should already be allocated for us
  for(int check_array=0; check_array<num_arrays; check_array++)
  {
    if(arrays[check_array] == NULL)
      return MB_FAILURE;
  }

  int node_index = 0;

  // lets get ready to iterate the entity sequences and copy data
  // we'll iterate the map in the sequence manager and
  // the entities range at the same time

  MBRange::const_iterator range_iter = entities.begin();
  MBRange::const_iterator range_iter_end = entities.end();

  std::map<MBEntityHandle, MBEntitySequence*>::const_iterator seq_iter
    = mMB->sequence_manager()->entity_map(MBVERTEX)->begin();
  
  std::map<MBEntityHandle, MBEntitySequence*>::const_iterator seq_iter_end
    = mMB->sequence_manager()->entity_map(MBVERTEX)->end();

  // lets find the entity sequence which holds the first entity
  std::map<MBEntityHandle, MBEntitySequence*>::const_iterator seq_iter_lookahead = seq_iter;
  seq_iter_lookahead++;
  for( ; seq_iter_lookahead != seq_iter_end && 
      seq_iter_lookahead->second->get_start_handle() < *range_iter; )
  {
    ++seq_iter;
    ++seq_iter_lookahead;
  }

  // a look ahead iterator
  MBRange::const_iterator range_iter_lookahead = range_iter;

  // our main loop
  for(; range_iter != range_iter_end && seq_iter != seq_iter_end; /* ++ is handled in loop*/ )
  {
    // find a range that fits in the current entity sequence
    for(; range_iter_lookahead != range_iter_end && 
        *range_iter_lookahead <= seq_iter->second->get_end_handle(); 
        ++range_iter_lookahead)
    {}

    int start_node_index = node_index;

    // get the coordinate array
    node_index = start_node_index;
    double* coord_array[3];
    static_cast<VertexEntitySequence*>(seq_iter->second)->get_coordinate_arrays(
        coord_array[0], coord_array[1], coord_array[2]);
    MBEntityHandle start_ent = seq_iter->second->get_start_handle();

    // for each of the entities in this entity sequence, copy data
    for(MBRange::const_iterator tmp_iter = range_iter; 
        tmp_iter != range_iter_lookahead;
        ++tmp_iter)
    {
      arrays[0][node_index] = coord_array[0][*tmp_iter - start_ent];
      arrays[1][node_index] = coord_array[1][*tmp_iter - start_ent];

      if( num_arrays == 3 )
        arrays[2][node_index] = coord_array[2][*tmp_iter - start_ent];

      ++node_index;
      tag_server->set_data(node_id_tag, *tmp_iter, &node_index);
    }

    // go to the next entity sequence
    ++seq_iter;
    // start with the next entities
    range_iter = range_iter_lookahead;
  }


  // we need to make sure we found all the nodes we were supposed to find
  // if not, we screwed up in this function
  assert(node_index == num_nodes);
  // if we hit this assert, then something wrong happened in MB.  The user only specifies meshsets to write out.
  // Therefore, we should always have valid entity handles and we should always be able to get the nodes we want.

  return error;

}

MBErrorCode MBWriteUtil::get_element_array(
    const int num_elements, 
    const int verts_per_element,
    MBTag node_id_tag,
    const MBRange& elements, 
    MBTag element_id_tag,
    int start_element_id,
    int* element_array)
{

  // check the data we got
  if(num_elements < 1)
    return MB_FAILURE;
  if(verts_per_element < 1)
    return MB_FAILURE;
  if(elements.empty())
    return MB_FAILURE;
  if(!element_array)
    return MB_FAILURE;

  TagServer* tag_server = mMB->tag_server();

  MBRange::const_iterator range_iter = elements.begin();
  MBRange::const_iterator range_iter_end = elements.end();

  std::map<MBEntityHandle, MBEntitySequence*>::const_iterator seq_iter, seq_iter_end;
  MBEntityType current_type = TYPE_FROM_HANDLE(*range_iter);
 
  seq_iter = mMB->sequence_manager()->entity_map(current_type)->begin();
  seq_iter_end = mMB->sequence_manager()->entity_map(current_type)->end();

  // lets find the entity sequence which holds the first entity
  std::map<MBEntityHandle, MBEntitySequence*>::const_iterator seq_iter_lookahead = seq_iter;
  seq_iter_lookahead++;
  for( ; seq_iter_lookahead != seq_iter_end && 
      seq_iter_lookahead->second->get_start_handle() < *range_iter; )
  {
    ++seq_iter;
    ++seq_iter_lookahead;
  }

  // a look ahead iterator
  MBRange::const_iterator range_iter_lookahead = range_iter;

  // our main loop
  for(; range_iter != range_iter_end && seq_iter != seq_iter_end; /* ++ is handled in loop*/ )
  {
    // find a range that fits in the current entity sequence
    for(; range_iter_lookahead != range_iter_end && 
        *range_iter_lookahead <= seq_iter->second->get_end_handle(); 
        ++range_iter_lookahead)
    {}
  
    if(current_type != TYPE_FROM_HANDLE(*range_iter))
    {
      current_type = TYPE_FROM_HANDLE(*range_iter);
      seq_iter = mMB->sequence_manager()->entity_map(current_type)->begin();
      seq_iter_end = mMB->sequence_manager()->entity_map(current_type)->end();

      // lets find the entity sequence which holds the first entity of this type
      std::map<MBEntityHandle, MBEntitySequence*>::const_iterator seq_iter_lookahead = seq_iter;
      seq_iter_lookahead++;
      for( ; seq_iter_lookahead != seq_iter_end && 
          seq_iter_lookahead->second->get_start_handle() < *range_iter; )
      {
        ++seq_iter;
        ++seq_iter_lookahead;
      }
    }

    int i = static_cast<ElementEntitySequence*>(seq_iter->second)->nodes_per_element();

    // get the connectivity array
    MBEntityHandle* conn_array = NULL;
    static_cast<ElementEntitySequence*>(seq_iter->second)->get_connectivity_array(conn_array);
 
    MBEntityHandle start_handle = seq_iter->second->get_start_handle();

    for(MBRange::const_iterator tmp_iter = range_iter; 
        tmp_iter != range_iter_lookahead;
        ++tmp_iter)
    {
      // set the element id tag
      tag_server->set_data(element_id_tag, *tmp_iter, &start_element_id);
      ++start_element_id;

      // for each node
      for(int j=0; j<i; j++)
      {
        MBEntityHandle node = *(conn_array + j + i*(*tmp_iter - start_handle));
        tag_server->get_data(node_id_tag, node, element_array);
        element_array++;
      }
    }

    // go to the next entity sequence
    ++seq_iter;
    // start with the next entities
    range_iter = range_iter_lookahead;
  }

  return MB_SUCCESS;
}

MBErrorCode MBWriteUtil::gather_nodes_from_elements(
      const MBRange& elements,
      const MBTag node_bit_mark_tag,
      MBRange& nodes
      )
{

  if(elements.empty())
    return MB_SUCCESS;

  TagServer* tag_server = mMB->tag_server();

  // see if we need to use our own marking tag
  MBTag exporting_nodes_tag = 0;
  if(node_bit_mark_tag)
    exporting_nodes_tag = node_bit_mark_tag;
  else
  {
    mMB->tag_create("__MBWriteUtil::exporting_nodes", 1, MB_TAG_BIT, 
                     exporting_nodes_tag, NULL);
  }
  

  MBRange::const_iterator range_iter = elements.begin();
  MBRange::const_iterator range_iter_end = elements.end();

  std::map<MBEntityHandle, MBEntitySequence*>::const_iterator seq_iter, seq_iter_end;
  MBEntityType current_type = TYPE_FROM_HANDLE(*range_iter);
 
  seq_iter = mMB->sequence_manager()->entity_map(current_type)->begin();
  seq_iter_end = mMB->sequence_manager()->entity_map(current_type)->end();

  // lets find the entity sequence which holds the first entity
  std::map<MBEntityHandle, MBEntitySequence*>::const_iterator seq_iter_lookahead = seq_iter;
  seq_iter_lookahead++;
  for( ; seq_iter_lookahead != seq_iter_end && 
      seq_iter_lookahead->second->get_start_handle() < *range_iter; )
  {
    ++seq_iter;
    ++seq_iter_lookahead;
  }

  // a look ahead iterator
  MBRange::const_iterator range_iter_lookahead = range_iter;

  // the x,y,z tag handles we need
  MBEntityHandle lower_bound = ~0, upper_bound = 0;
  
  // our main loop
  for(; range_iter != range_iter_end && seq_iter != seq_iter_end; /* ++ is handled in loop*/ )
  {
    // find a range that fits in the current entity sequence
    for(; range_iter_lookahead != range_iter_end && 
        *range_iter_lookahead <= seq_iter->second->get_end_handle(); 
        ++range_iter_lookahead)
    {}
  
    if(current_type != TYPE_FROM_HANDLE(*range_iter))
    {
      current_type = TYPE_FROM_HANDLE(*range_iter);
      seq_iter = mMB->sequence_manager()->entity_map(current_type)->begin();
      seq_iter_end = mMB->sequence_manager()->entity_map(current_type)->end();

      // lets find the entity sequence which holds the first entity of this type
      std::map<MBEntityHandle, MBEntitySequence*>::const_iterator seq_iter_lookahead = seq_iter;
      seq_iter_lookahead++;
      for( ; seq_iter_lookahead != seq_iter_end && 
          seq_iter_lookahead->second->get_start_handle() < *range_iter; )
      {
        ++seq_iter;
        ++seq_iter_lookahead;
      }
    }

    int i = static_cast<ElementEntitySequence*>(seq_iter->second)->nodes_per_element();

    // get the connectivity array
    MBEntityHandle* conn_array = NULL;
    static_cast<ElementEntitySequence*>(seq_iter->second)->get_connectivity_array(conn_array);
 
    MBEntityHandle start_handle = seq_iter->second->get_start_handle();

    for(MBRange::const_iterator tmp_iter = range_iter; 
        tmp_iter != range_iter_lookahead;
        ++tmp_iter)
    {
      // for each node
      for(int j=0; j<i; j++)
      {
        MBEntityHandle node = *(conn_array + j + i*(*tmp_iter - start_handle));
        if(node < lower_bound)
          lower_bound = node;
        if(node > upper_bound)
          upper_bound = node;
        unsigned char bit = 0x1;
        tag_server->set_data(exporting_nodes_tag, &node, 1, &bit);
      }
    }

    // go to the next entity sequence
    ++seq_iter;
    // start with the next entities
    range_iter = range_iter_lookahead;
  }

  // we can get a REALLY long loop if lower_bound is zero
  assert(lower_bound != 0);
  // gather up all the nodes
  for(; upper_bound >= lower_bound; --upper_bound)
  {
    unsigned char node_marked=0;
    tag_server->get_data(exporting_nodes_tag, &upper_bound, 1, &node_marked);
    if(node_marked == 0x1)
      nodes.insert(upper_bound);
  }

  // clean up our own marking tag
  if(node_bit_mark_tag == 0)
    mMB->tag_delete(exporting_nodes_tag);

  return MB_SUCCESS;

}


MBErrorCode MBWriteUtil::report_error( const std::string& error )
{
  if(mError)
    return mError->set_last_error(error);
  else
    return MB_FAILURE;
}


MBErrorCode MBWriteUtil::report_error( const char* error, ... )
{
  va_list args;
  va_start(args, error);
  MBErrorCode result = mError->set_last_error(error, args);
  va_end(args);
  return result;
}
