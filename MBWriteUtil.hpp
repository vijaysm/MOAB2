
#ifndef MB_WRITE_UTIL_HPP
#define MB_WRITE_UTIL_HPP

#ifndef IS_BUILDING_MB
#error "MBWriteUtil.hpp isn't supposed to be included into an application"
#endif

#include "MBWriteUtilIface.hpp"

class MBCore;
class MBError;

class MBWriteUtil : public MBWriteUtilIface
{
private:
  //! pointer to the MBCore
  MBCore* mMB;
  MBError* mError;
public:

  //! constructor takes MBCore pointer
  MBWriteUtil(MBCore* mdb, MBError* error_handler);

  //! destructor
  ~MBWriteUtil(){}

  //! gets arrays for coordinate data from the MB
  MBErrorCode get_node_arrays(
      const int num_arrays,
      const int num_nodes, 
      const MBRange& entities,
      MBTag node_id_tag,
      std::vector<double*>& arrays
      );

  //! get array for connectivity data from the MB
  MBErrorCode get_element_array(
      const int num_elements, 
      const int verts_per_element,
      MBTag node_id_tag,
      const MBRange& entities, 
      MBTag element_id_tag,
      int start_element_id,
      int* array
      );
  
  //! get a set of nodes that represent a set of elements
  MBErrorCode gather_nodes_from_elements(
      const MBRange& elements,
      const MBTag node_bit_mark_tag,
      MBRange& nodes
      );
  
 
  //! tell MB there was an error when writing the mesh
  //! it makes sense to have this as long as MBInterface has a write_mesh function
  MBErrorCode report_error( const std::string& error );

  MBErrorCode report_error( const char* error, ... );

};

#endif


