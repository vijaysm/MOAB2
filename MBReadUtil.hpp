
#ifndef MB_READ_UTIL_HPP
#define MB_READ_UTIL_HPP

#ifndef IS_BUILDING_MB
#error "MBReadUtil.hpp isn't supposed to be included into an application"
#endif

#include "MBReadUtilIface.hpp"

class MBCore;
class MBError;

class MBReadUtil : public MBReadUtilIface
{
private:
  //! pointer to the MBCore
  MBCore* mMB;
  MBError* mError;
public:

  //! constructor takes MBCore pointer
  MBReadUtil(MBCore* mdb, MBError* error_handler);

  //! destructor
  ~MBReadUtil(){}

  //! gets arrays for coordinate data from the MB
  MBErrorCode get_node_arrays(
      const int num_arrays,
      const int num_nodes, 
      const int preferred_start_id,
      MBEntityHandle& actual_start_handle, 
      std::vector<double*>& arrays
      );

  //! get array for connectivity data from the MB
  MBErrorCode get_element_array(
      const int num_elements, 
      const int verts_per_element,
      const MBEntityType mdb_type,
      int preferred_start_id, 
      MBEntityHandle& actual_start_handle, 
      MBEntityHandle*& array
      );
 
  //! tell MB which elements have been added to the database
  MBErrorCode update_adjacencies(
      const MBEntityHandle start_handle,
      const int number_elements,
      const int number_vertices_per_element,
      const MBEntityHandle* conn_array);

  //! tell MB there was an error when reading the mesh
  //! it makes sense to have this as long as MBInterface has a load_mesh function
  MBErrorCode report_error( const std::string& error );

  MBErrorCode report_error( const char* error, ... );

};

#endif


