
#ifndef MB_WRITE_UTIL_IFACE_HPP
#define MB_WRITE_UTIL_IFACE_HPP


#include <vector>
#include "MBInterface.hpp"
#include <stdarg.h>

// Interface implemented in MB which provides memory for mesh reading utilities
class MB_DLL_EXPORT MBWriteUtilIface
{
public:
  MBWriteUtilIface(){}
  virtual ~MBWriteUtilIface(){}

  //! get node arrays
  //! given number of arrays requested (1 for each dimension)
  //! given number of nodes for space requested
  //! given entities
  //! given node id tag handle.  MB will assign int values to the nodes 
  //!   it will write out.  These values start at 1.
  //! MB fills in the coordinates data
  //!   arrays can be indexed as ---
  //!     arrays[0] is an x-coordinate array
  //!     arrays[1] is a y-coordinate array
  //!     arrays[0][0] is the x coordinate of the first node
  virtual MBErrorCode get_node_arrays(
      const int num_arrays,
      const int num_nodes, 
      const MBRange& entities, 
      MBTag node_id_tag,
      std::vector<double*>& arrays
      ) = 0;

  //! get element array
  //! given number of elements for space requested
  //! given vertices per element
  //! tag to get integer node ids from entity handles, these integers end up in array.
  //! given the entities
  //! tag for MB to set element ids. Can be zero in which case ids aren't set
  //!   the element start id for which element ids are set.  Ids start at 1.
  //! MB fills in connectivity data
  //!  array can be indexed as ---
  //!    array[0] is the first node of the first element
  virtual MBErrorCode get_element_array(
      const int num_elements, 
      const int verts_per_element,
      MBTag node_id_tag,
      const MBRange& entities, 
      MBTag element_id_tag,
      int start_element_id,
      int* array
      ) = 0;

  //! get a set of nodes that represent a set of elements
  //! given elements
  //! given node_bit_mark_tag.  Will marks nodes with 0x1.  If zero, nodes aren't marked.
  //! returns nodes for elements, append to nodes list if list has things in it
  virtual MBErrorCode gather_nodes_from_elements(
      const MBRange& elements,
      const MBTag node_bit_mark_tag,
      MBRange& nodes
      ) = 0;

    //! assign ids to input elements starting with start_id, written to id_tag
    //! if id_tag is zero, assigns to GLOBAL_ID_TAG_NAME and passes back
  virtual MBErrorCode assign_ids(MBRange &elements,
                                 MBTag id_tag,
                                 const int start_id) = 0;

  //! if an error occured when reading the mesh, report it to MB
  //! it makes sense to have this as long as MBInterface has a write_mesh function
  virtual MBErrorCode report_error( const std::string& error ) = 0;
  
  //! overloaded report_error behaves like the above
  virtual MBErrorCode report_error( const char* error, ... ) = 0;

};

#endif 


