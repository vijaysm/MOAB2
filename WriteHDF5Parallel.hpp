/** 
 * \class WriteHDF5Parallel
 * \brief Write MOAB HDF5 file in parallel.
 * \author Jason Kraftcheck
 * \data   22 July 2004
 */

#ifndef WRITE_HDF5_PARALLEL_HPP
#define WRITE_HDF5_PARALLEL_HPP

#include "WriteHDF5.hpp"
#include "MBTagConventions.hpp"
#include <mpi.h>

class MB_DLL_EXPORT WriteHDF5Parallel : public WriteHDF5
{
  public:
    
    WriteHDF5Parallel( MBInterface* iface,
                       const char*[] multiproc_set_tags =
                       { MATERIAL_SET_TAG_NAME,
                         DIRICHLET_SET_TAG_NAME,
                         NEUMANN_SET_TAG_NAME,
                         0 } ) ;
    
  
  protected:
  
      //! Called by normal (non-parallel) writer.  Sets up
      //! necessary data for parallel write.
    virtual MBErrorCode create_file( const char* filename,
                                     bool overwrite,
                                     std::vector<std::string>& qa_records,
                                     int dimension = 3 );
                                     
      //! An array of interface mesh which is to be written by
      //! remote processors.  Indexed by MPI rank (processor number).
    std::vector<MBRange> remoteMesh;
    
      //! The meshsets which a) exist on this processor and
      //! b) span multiple processors.
    MBRange sharedSets;
    std::vector<long> sharedContentsOffests, sharedChildrenOffsets;
   
    MBErrorCode gather_interface_meshes();
    MBErrorCode communicate_remote_ids(MBEntityType type);
    void sort_tags_by_name();
    
    MBErrorCode create_node_table( int dimension );
    MBErrorCode negotiate_type_list();
    MBErrorCode create_element_tables();
    MBErrorCode create_adjacency_tables();
    MBErrorCode negotiate_shared_meshsets( long* offests );
    MBErrorCode create_meshset_tables();
    MBErrorCode create_tag_tables();
    
    
  private:
    
    int numProc;
    int myRank;
    
    std::vector<std::string> multiProcSetTags;
    
};

#endif
