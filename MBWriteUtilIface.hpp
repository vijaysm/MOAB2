
#ifndef MB_WRITE_UTIL_IFACE_HPP
#define MB_WRITE_UTIL_IFACE_HPP


#include <vector>
#include "MBInterface.hpp"
#include <stdarg.h>

//! Interface implemented in MOAB which provides memory for mesh reading utilities
class MB_DLL_EXPORT MBWriteUtilIface
{
public:

    //! constructor
  MBWriteUtilIface(){}

    //! destructor
  virtual ~MBWriteUtilIface(){}

    //! Given information about the nodes to be written, and pointers to memory
    //! to which coordinates will be written, writes coordinate data there, and
    //! also assigns global ids to nodes & writes to a tag
    //! \param num_arrays Number of coordinate arrays requested
    //! \param num_nodes Number of nodes to be written
    //! \param entities Range of nodes to be written
    //! \param node_id_tag Tag used to write ids to nodes
    //! \param start_node_id Starting value for node ids
    //! \param arrays Pointers to memory where coordinate data will be written
    //! \return status Return status
  virtual MBErrorCode get_node_arrays(
    const int num_arrays,
    const int num_nodes, 
    const MBRange& entities, 
    MBTag node_id_tag,
    const int start_node_id,
    std::vector<double*>& arrays
    ) = 0;

    //! Given information about elements to be written and a pointer to memory
    //! where connectivity for those elements should be written, writes connectivity
    //! to that memory; uses node ids stored in a tag during call to <em>get_node_arrays</em>
    //! function
    //! \param num_elements Number of elements to be written
    //! \param verts_per_element Number of vertices per element
    //! \param node_id_tag Tag used to store node ids
    //! \param entities Range of elements to be written
    //! \param element_id_tag Tag which should be used to store element ids
    //! \param start_element_id Starting value for element ids
    //! \param array Pointer to memory where connectivity data will be written
    //! \return status Return status
  virtual MBErrorCode get_element_array(
    const int num_elements, 
    const int verts_per_element,
    MBTag node_id_tag,
    const MBRange& entities, 
    MBTag element_id_tag,
    int start_element_id,
    int* array
    ) = 0;

    //! given elements to be written, gather all the nodes which define those elements
    //! \param elements Range of elements to be written
    //! \param node_bit_mark_tag Bit tag to use to identify nodes
    //! \param nodes Range of nodes gathered from elements (returned)
    //! \return status Return status
  virtual MBErrorCode gather_nodes_from_elements(
    const MBRange& elements,
    const MBTag node_bit_mark_tag,
    MBRange& nodes
    ) = 0;

    //! assign ids to input entities starting with start_id, written to id_tag
    //! if id_tag is zero, assigns to GLOBAL_ID_TAG_NAME
    //! \param elements Entities to be written
    //! \param id_tag Tag used to store entity id
    //! \param start_id Starting value for entity ids
    //! \return status Return status
  virtual MBErrorCode assign_ids(MBRange &elements,
                                 MBTag id_tag,
                                 const int start_id) = 0;

    //! if an error occured when reading the mesh, report it to MB
    //! it makes sense to have this as long as MBInterface has a write_mesh function
    //! \return status Return status
  virtual MBErrorCode report_error( const std::string& error ) = 0;
  
    //! overloaded report_error behaves like the above
    //! \return status Return status
  virtual MBErrorCode report_error( const char* error, ... ) = 0;

};

#endif 


