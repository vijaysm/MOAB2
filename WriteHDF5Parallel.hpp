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

struct RemoteSetData;

class MB_DLL_EXPORT WriteHDF5Parallel : public WriteHDF5
{
  public:
    
      /** Constructor
       *\param multiproc_set_tags Null-terminated list strings.
       *
       *multiproc_set_tags is a null-terminated list of tag names.
       *Each tag specified must have an native integer (int) data 
       *type.  The tag data is used to identify meshsets that span
       *multiple processors such that they are written as a single
       *meshset in the resulting file.  The default behavior if the
       *argument is null is to use MATERIAL_SET_TAG_NAME, 
       *DIRICHLET_SET_TAG_NAME, and NEUMANN_SET_TAG_NAME.  To 
       *disable this functionality entirely, pass in an array containing
       *only the terminated NULL.
       */
    WriteHDF5Parallel( MBInterface* iface,
                       const char** multiproc_set_tags = 0 );
    
  
  protected:
  
      //! Called by normal (non-parallel) writer.  Sets up
      //! necessary data for parallel write.
    virtual MBErrorCode create_file( const char* filename,
                                     bool overwrite,
                                     std::vector<std::string>& qa_records,
                                     int dimension = 3 );
    
      //! Figure out which mesh local mesh is duplicated on
      //! remote processors and which processor will write
      //! that mesh.
    MBErrorCode gather_interface_meshes();
    
      //! For entities that will be written by another 
      //! processor, get the file Ids that will be assigned
      //! to those so they can be referenced by
      //! entities to be written on this processor.
    MBErrorCode communicate_remote_ids(MBEntityType type);
    
      //! Sort the list of tag information in the parent
      //! class by name so all procs have them in the same
      //! order.
    void sort_tags_by_name();
    
      //! Create the node table in the file.
    MBErrorCode create_node_table( int dimension );
    
      //! Communicate with other processors to negotiate 
      //! the types of elements that will be written
      //! (the union of the types defined on each proc.)
    MBErrorCode negotiate_type_list();
    
      //! Create tables to hold element connectivity
    MBErrorCode create_element_tables();
    
      //! Create tables to hold element adjacencies.
    MBErrorCode create_adjacency_tables();
    
      //! Identify and set up meshsets that span multiple
      //! processors.
      //!\param offsets Output array of three values.
    MBErrorCode negotiate_shared_meshsets( long* offsets );
    
      //! Setup meshsets spanning multiple processors
    MBErrorCode get_remote_set_data( const char* tagname,
                                     RemoteSetData& data,
                                     long& offset );
    
      //! Determine offsets in contents and children tables for 
      //! meshsets shared between processors.
    MBErrorCode negotiate_remote_set_contents( RemoteSetData& data,
                                               long* offsets );
    
      //! Create tables for mesh sets
    MBErrorCode create_meshset_tables();
    
      //! Write tag descriptions and create tables to hold tag data.
    MBErrorCode create_tag_tables();
    
      //! Mark multiple-processor meshsets with correct file Id
      //! from the set description offset stored in that tag by
      //! negotiate_shared_meshsets(..).
    MBErrorCode fix_remote_set_ids( RemoteSetData& data, long first_id );
      
      //! Write set descriptions for multi-processor meshsets.
      //! Virtual function called by non-parallel code after
      //! the normal (single-processor) meshset descriptions have
      //! been written.
    MBErrorCode write_shared_set_descriptions( hid_t table );
       
      //! Write set contents for multi-processor meshsets.
      //! Virtual function called by non-parallel code after
      //! the normal (single-processor) meshset contents have
      //! been written.
    MBErrorCode write_shared_set_contents( hid_t table );
       
      //! Write set children for multi-processor meshsets.
      //! Virtual function called by non-parallel code after
      //! the normal (single-processor) meshset children have
      //! been written.
    MBErrorCode write_shared_set_children( hid_t table );
       
      //! Write set children for multi-processor meshsets.
      //! Virtual function called by non-parallel code after
      //! the normal (single-processor) meshset children have
      //! been written.
    MBErrorCode write_shared_set_parents( hid_t table );
  
      //! Virtual function overridden from WriteHDF5.  
      //! Release memory by clearing member lists.
    MBErrorCode write_finished();
    
      //! Remove any remote mesh entities from the passed range.
    void remove_remote_entities( MBRange& range );
    void remove_remote_entities( std::vector<MBEntityHandle>& vect );
    
  private:
    
      //! MPI environment
    int numProc, myRank;
                                     
      //! An array of interface mesh which is to be written by
      //! remote processors.  Indexed by MPI rank (processor number).
    std::vector<MBRange> remoteMesh;
    
      //! Tag names for identifying multi-processor meshsets
    std::vector<std::string> multiProcSetTags;
    
      //! Struct describing a multi-processor meshset
    struct ParallelSet {
      MBEntityHandle handle;
      long contentsOffset;
      long childrenOffset;
      long parentsOffset;
      long contentsCount;
      long childrenCount;
      long parentsCount;
      bool description;
    };
    
      //! List of multi-processor meshsets
    std::list<ParallelSet> parallelSets;
    
    void printrange( MBRange& );
};

#endif
