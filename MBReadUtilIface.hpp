
#ifndef MB_READ_UTIL_IFACE_HPP
#define MB_READ_UTIL_IFACE_HPP


#include <vector>
#include "MBInterface.hpp"
#include <stdarg.h>

// Interface implemented in MB which provides memory for mesh reading utilities
class MB_DLL_EXPORT MBReadUtilIface
{
public:
  MBReadUtilIface(){}
  virtual ~MBReadUtilIface(){}

  //! get node arrays
  //! given number of arrays requested (1 for each dimension)
  //! given number of nodes for space requested
  //! given a prefferrable start handle
  //! returns the actual start handle
  //! returns the arrays
  //!   arrays can be indexed as ---
  //!     arrays[0] is an x-coordinate array
  //!     arrays[1] is a y-coordinate array
  //!     arrays[0][0] is the x coordinate of the first node
  virtual MBErrorCode get_node_arrays(
      const int num_arrays,
      const int num_nodes, 
      const int preferred_start_id,
      MBEntityHandle& actual_start_handle, 
      std::vector<double*>& arrays
      ) = 0;

  //! get element array
  //! given number of elements for space requested
  //! given vertices per element
  //! given a preferrable start handle
  //! returns the actual start handle
  //! returns an array
  //!  array can be indexed as ---
  //!    array[0] is the first node of the first element
  virtual MBErrorCode get_element_array(
      const int num_elements, 
      const int verts_per_element,
      const MBEntityType mdb_type,
      int preferred_start_id, 
      MBEntityHandle& actual_start_handle, 
      MBEntityHandle*& array
      ) = 0;

  //! update adjacencies
  //! given a start handle, number of elements, number of vertices per element,
  //! and a connectivity array for entities, adjacency information will be updated
  //! in MB.  Also think of it as a way of Readers telling MB what elements are 
  //! new because we aren't using the MBInterface to create elements.
  virtual MBErrorCode update_adjacencies(
      const MBEntityHandle start_handle,
      const int number_elements,
      const int number_vertices_per_element,
      const MBEntityHandle* conn_array
      ) = 0;


  //! if an error occured when reading the mesh, report it to MB
  //! it makes sense to have this as long as MBInterface has a load_mesh function
  virtual MBErrorCode report_error( const std::string& error ) = 0;

  //! overloaded report_error behaves like the above
  virtual MBErrorCode report_error( const char* error, ... ) = 0;

};

#endif 


