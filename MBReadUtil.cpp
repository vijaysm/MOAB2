#ifdef WIN32
#ifdef _DEBUG
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif
#endif

#include "MBReadUtil.hpp"
#include "MBCore.hpp"
#include "AEntityFactory.hpp"
#include "MBError.hpp"
#include "EntitySequenceManager.hpp"


MBReadUtil::MBReadUtil(MBCore* mdb, MBError* error_handler) 
    : MBReadUtilIface(), mMB(mdb), mError(error_handler)
{
}


MBErrorCode MBReadUtil::get_node_arrays(
    const int /*num_arrays*/,
    const int num_nodes, 
    const int preferred_start_id,
    MBEntityHandle& actual_start_handle, 
    std::vector<double*>& arrays)
{

  MBErrorCode error;
  MBEntitySequence* seq = 0;

  MBEntityHandle preferred_start_handle;
  static int err;
  preferred_start_handle = CREATE_HANDLE(MBVERTEX, preferred_start_id, err);
 
  // create an entity sequence for these nodes 
  error = mMB->sequence_manager()->create_entity_sequence(
      MBVERTEX, num_nodes, 0, preferred_start_handle, actual_start_handle,
      seq);

  if(error != MB_SUCCESS)
    return error;

  arrays.resize(3);

  error = static_cast<VertexEntitySequence*>(seq)->get_coordinate_arrays(arrays[0], arrays[1], arrays[2]);
  
  return error;
}

MBErrorCode MBReadUtil::get_element_array(
    const int num_elements, 
    const int verts_per_element,
    const MBEntityType mdb_type,
    int preferred_start_id, 
    MBEntityHandle& actual_start_handle, 
    MBEntityHandle*& array)
{

  MBErrorCode error;
  MBEntitySequence* seq;

  // make an entity sequence to hold these elements
  error = mMB->sequence_manager()->create_entity_sequence(
      mdb_type, num_elements, verts_per_element, preferred_start_id, actual_start_handle, seq);

  // get an array for the connectivity
  error = static_cast<ElementEntitySequence*>(seq)->get_connectivity_array(array);

  return error;
  
}



MBErrorCode MBReadUtil::update_adjacencies(
      const MBEntityHandle start_handle,
      const int number_elements,
      const int number_vertices_per_element,
      const MBEntityHandle* conn_array)
{

  MBEntityHandle tmp_hndl = start_handle;
  AEntityFactory* adj_fact = mMB->a_entity_factory();

  // iterator over the elements and update adjacency information
  if(adj_fact != NULL && adj_fact->vert_elem_adjacencies())
  {
    int j=0;
    for(int i=0; i<number_elements; i++)
    {
      adj_fact->notify_create_entity( tmp_hndl, (conn_array+j), number_vertices_per_element);
      tmp_hndl++;
      j+=number_vertices_per_element;
    }
  }
  return MB_SUCCESS;
}



MBErrorCode MBReadUtil::report_error( const std::string& error )
{
  if(mError)
    return mError->set_last_error(error);
  else
    return MB_FAILURE;
}


MBErrorCode MBReadUtil::report_error( const char* error, ... )
{
  va_list args;
  va_start(args, error);
  MBErrorCode result = mError->set_last_error(error, args);
  va_end(args);
  return result;
}


