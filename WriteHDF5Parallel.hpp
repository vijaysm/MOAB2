/** 
 * \class WriteHDF5Parallel
 * \brief Write MOAB HDF5 file in parallel.
 * \author Jason Kraftcheck
 * \data   22 July 2004
 */

#ifndef WRITE_HDF5_PARALLEL_HPP
#define WRITE_HDF5_PARALLEL_HPP

#include "WriteHDF5.hpp"
#include <mpi.h>

class MB_DLL_EXPORT WriteHDF5Parallel : public WriteHDF5
{
  public:
    
    WriteHDF5Parallel( MBInterface* iface ) : WriteHDF5( iface ) {}
    
  
  protected:
  
    virtual MBErrorCode create_file( const char* filename,
                                     bool overwrite,
                                     std::vector<std::string>& qa_records,
                                     int dimension = 3 );
                                     
    struct Interface {
      MBRange range;
      MBEntityHandle set;
      int rank;
      bool remote;
      bool operator<( Interface& other )
        { return rank < other.rank; }
      
      std::vector<MBEntityHandle> buffer;
      MPI_Request req;
    };
    
    std::list<Interface> interfaceList;
    
    MBErrorCode gather_interface_meshes();
    MBErrorCode communicate_remote_ids();
    void sort_tags_by_name();

    MBErrorCode mesh_from_geom_sets( const MBRange& sets, MBRange& mesh );
    MBErrorCode start_send( Interface& );
    MBErrorCode start_recv( Interface& );
    MBErrorCode finish_recv( Interface& );
};

#endif
