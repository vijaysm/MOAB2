
#ifndef MB_READ_UTIL_IFACE_HPP
#define MB_READ_UTIL_IFACE_HPP


#include <vector>
#include "MBInterface.hpp"
#include <stdarg.h>

//! Interface implemented in MOAB which provides memory for mesh reading utilities
class MB_DLL_EXPORT MBReadUtilIface
{
public:

    //! constructor 
  MBReadUtilIface(){}

    //! destructor
  virtual ~MBReadUtilIface(){}

    //! Given a requested number of vertices and number of coordinates, returns
    //! memory space which will be used to store vertex coordinates and information
    //! about what handles those new vertices are assigned; allows direct read of 
    //! coordinate data into memory
    //! \param num_arrays Number of node position arrays requested
    //! \param num_nodes Number of nodes
    //! \param preferred_start_id Preferred integer id starting value
    //! \param actual_start_handle Actual starting id value
    //! \param arrays STL vector of double*'s, point to memory storage to be used for 
    //!     these vertices
    //! \return status Success/failure of this call
  virtual MBErrorCode get_node_arrays(
    const int num_arrays,
    const int num_nodes, 
    const int preferred_start_id,
    MBEntityHandle& actual_start_handle, 
    std::vector<double*>& arrays
    ) = 0;

    //! Given requested number of elements, element type, and number of
    //! elements, returns pointer to memory space allocated to store connectivity
    //! of those elements; allows direct read of connectivity data into memory
    //! \param num_elements Number of elements being requested
    //! \param verts_per_element Number of vertices per element (incl. higher-order nodes)
    //! \param mdb_type Element type
    //! \param preferred_start_id Preferred integer id for first element
    //! \param actual_start_handle Actual integer id for first element (returned)
    //! \param array Pointer to memory allocated for storing connectivity for these elements
    //! \return status Success/failure of this call
  virtual MBErrorCode get_element_array(
    const int num_elements, 
    const int verts_per_element,
    const MBEntityType mdb_type,
    int preferred_start_id, 
    MBEntityHandle& actual_start_handle, 
    MBEntityHandle*& array
    ) = 0;

    //! update adjacencies
    //! given information about new elements, adjacency information will be updated
    //! in MOAB.  Think of this function as a way of Readers telling MOAB what elements are 
    //! new because we aren't using the MBInterface to create elements.
    //! \param start_handle Handle of first new element
    //! \param number_elements Number of new elements
    //! \param number_vertices_per_element Number of vertices in each new element
    //! \param conn_array Connectivity of new elements
    //! \return status Success/failure of this call
  virtual MBErrorCode update_adjacencies(
    const MBEntityHandle start_handle,
    const int number_elements,
    const int number_vertices_per_element,
    const MBEntityHandle* conn_array
    ) = 0;


    //! if an error occured when reading the mesh, report it to MOAB
    //! it makes sense to have this as long as MBInterface has a load_mesh function
  virtual MBErrorCode report_error( const std::string& error ) = 0;

    //! overloaded report_error behaves like the above
  virtual MBErrorCode report_error( const char* error, ... ) = 0;

};

#endif 


