
#ifndef WRITE_HDF5_PARALLEL_HPP
#define WRITE_HDF5_PARALLEL_HPP

#include "WriteHDF5.hpp"
#include <map>

namespace moab {

struct RemoteSetData;
class ParallelComm;
class IODebugTrack;

/** 
 * \brief Write MOAB HDF5 file in parallel.
 * \author Jason Kraftcheck
 * \data   22 July 2004
 */
class WriteHDF5Parallel : public WriteHDF5
{
  public:

    static WriterIface* factory( Interface* );
    
      /** Consturctor
       *
       * This constructor will automatically register the tags for
       * material set (block), dirichlet set (nodeset), neumann set
       * (sideset), and geometry grouping sets for use in identifying
       * sets that are shared across multiple processors.  To explicitly
       * disable this functionality, call one of the other construtors
       * with an empty list of tags.
       */
    WriteHDF5Parallel( Interface* iface );
     
    
      /** Constructor
       *\param multiproc_set_tags Null-terminated list strings.
       *
       * multiproc_set_tags is a null-terminated list of tag names.
       * Each tag specified must have an native integer (int) data 
       * type.  The tag data is used to identify meshsets that span
       * multiple processors such that they are written as a single
       * meshset in the resulting file.  
       *
       * NOTE: This list must be identical on all processors, including
       *       the order!
       */
    WriteHDF5Parallel( Interface* iface,
                       const std::vector<std::string>& multiproc_set_tags );

  virtual ~WriteHDF5Parallel();
  
    /**\brief Define tags used to identify sets spanning multiple procesors */
    class MultiProcSetTags {
      friend class WriteHDF5Parallel;
      public:

        /**Specify the name of a tag used to identify parallel entity sets.
         * The tag must have an native integer (int) data type.  The value
         * of the tag will be used to match sets on different processors.
         */
      void add( const std::string& name );
 
        /**Specify separate tags for identifying parallel entity sets and
         * matching them across processors.
         *\param filter_name The name of a tag used to identify parallel entity sets
         *\param value_name  The name of a tag having a native integer (int) data
         *                   type.  The value of this tag is used as an ID to match
         *                   entity sets on different processors.
         */
      void add( const std::string& filter_name, const std::string& value_name );
 
        /**Specify separate tags for identifying parallel entity sets and
         * matching them across processors.
         *\param filter_name The name of a tag used to identify parallel entity sets.
         *                   The data type of this tag must be a native integer (int).
         *\param filter_value The value of the filter_name tag to use to identify
         *                   parallel entity sets.
         *\param value_name  The name of a tag having a native integer (int) data
         *                   type.  The value of this tag is used as an ID to match
         *                   entity sets on different processors.
         */
      void add( const std::string& filter_name, int filter_value, const std::string& value_name );
      
      private:
      class Data;
      std::vector<Data> list;
    };
     
      /** Constructor
       *\param multiproc_set_tags Data used to identify sets spanning multiple processors.
       *                          NOTE:  This must be identical on all processors, including
       *                          the order in which tags were added to the object!
       */
    WriteHDF5Parallel( Interface* iface, const MultiProcSetTags& multiproc_set_tags );
      
    
  
  protected:
    
    virtual void debug_barrier_line(int lineno);

  
    virtual void print_times( const double* times ) const;
  
      //! Called by normal (non-parallel) writer.  Sets up
      //! necessary data for parallel write.
    virtual ErrorCode parallel_create_file( const char* filename,
                                     bool overwrite,
                                     const std::vector<std::string>& qa_records,
                                     const Tag* user_tag_list = 0,
                                     int user_tag_count = 0,
                                     int dimension = 3,
                                     int pcomm_no = 0,
                                     double* times = 0);
    
      //! Figure out which mesh local mesh is duplicated on
      //! remote processors and which processor will write
      //! that mesh.
      //!\param non_local_ents Output list of entities that are not to
      //!                  be written by this processor but are
      //!                  referenced by other entities that are
      //!                  to be written.
    ErrorCode gather_interface_meshes( Range& non_local_ents );
    
      //! For entities that will be written by another 
      //! processor but are referenced by entities on this
      //! processor, get the file Ids that will be assigned
      //! to those so they can be referenced by
      //! entities to be written on this processor.
      //!\param non_local_ents List of entities that are not to
      //!                  be written by this processor but are
      //!                  referenced by other entities that are
      //!                  to be written.
    ErrorCode exchange_file_ids( const Range& non_local_ents );
    
      //! Create the node table in the file.
    ErrorCode create_node_table( int dimension );
    
      //! Communicate with other processors to negotiate 
      //! the types of elements that will be written
      //! (the union of the types defined on each proc.)
    ErrorCode negotiate_type_list();
    
      //! Create tables to hold element connectivity
    ErrorCode create_element_tables();
    
      //! Create tables to hold element adjacencies.
    ErrorCode create_adjacency_tables();
    
      //! Setup meshsets spanning multiple processors
    ErrorCode get_remote_set_data( const MultiProcSetTags::Data& tag,
                                     RemoteSetData& data,
                                     long& offset );
    
      //! Determine offsets in contents and children tables for 
      //! meshsets shared between processors.
    ErrorCode negotiate_remote_set_contents( RemoteSetData& data,
                                               long* offsets );
    
      //! Create tables for mesh sets
    ErrorCode create_meshset_tables(double* times);
    
      //! Write tag descriptions and create tables to hold tag data.
    ErrorCode create_tag_tables();
    
      //! Mark multiple-processor meshsets with correct file Id
      //! from the set description offset stored in that tag by
      //! negotiate_shared_meshsets(..).
    ErrorCode set_shared_set_ids( RemoteSetData& data, long& start_id );
      
      //! Write set descriptions for multi-processor meshsets.
      //! Virtual function called by non-parallel code after
      //! the normal (single-processor) meshset descriptions have
      //! been written.
    ErrorCode write_shared_set_descriptions( hid_t table, IODebugTrack* );
       
      //! Write set contents for multi-processor meshsets.
      //! Virtual function called by non-parallel code after
      //! the normal (single-processor) meshset contents have
      //! been written.
    ErrorCode write_shared_set_contents( hid_t table, IODebugTrack* );
       
      //! Write set children for multi-processor meshsets.
      //! Virtual function called by non-parallel code after
      //! the normal (single-processor) meshset children have
      //! been written.
    ErrorCode write_shared_set_children( hid_t table, IODebugTrack* );
       
      //! Write set children for multi-processor meshsets.
      //! Virtual function called by non-parallel code after
      //! the normal (single-processor) meshset children have
      //! been written.
    ErrorCode write_shared_set_parents( hid_t table, IODebugTrack* );
  
      //! Virtual function overridden from WriteHDF5.  
      //! Release memory by clearing member lists.
    ErrorCode write_finished();
   
      //! Remove any remote mesh entities from the passed range.
    void remove_remote_entities( EntityHandle relative, Range& range );
    void remove_remote_entities( EntityHandle relative, std::vector<EntityHandle>& vect );
    void remove_remote_sets( EntityHandle relative, Range& range );
    void remove_remote_sets( EntityHandle relative, std::vector<EntityHandle>& vect );
  
      //! get any existing tags which aren't excluded and add to shared set tags
    ErrorCode get_sharedset_tags();

    ErrorCode append_serial_tag_data( std::vector<unsigned char>& buffer,
                                      const WriteHDF5::TagDesc& tag );
  
    ErrorCode write_shared_set_data( hid_t table,
                                     WriteUtilIface::EntityListType which_data,
                                     IODebugTrack* dbg_track );

      //! helper function for create_tag_tables
    ErrorCode check_serial_tag_data( 
                               const std::vector<unsigned char>& buffer,
                               std::vector<TagDesc*>* missing = 0,
                               std::vector<TagDesc*>* newlist = 0 );
  
  private:
    
      //! Tag names for identifying multi-processor meshsets
    MultiProcSetTags multiProcSetTags;
    
      //! Struct describing a multi-processor meshset
    struct ParallelSet {
      EntityHandle handle;// set handle on this processor
      long contentsOffset;  // offset in table at which to write set contents
      long childrenOffset;  // offset in table at which to write set children
      long parentsOffset;   // offset in table at which to write set parents
      long contentsCount;   // total size of set contents (all processors)
      long childrenCount;   // total number of set children (all processors)
      long parentsCount;    // total numoer of set parents (all processors)
      bool description;     // true if this processor 'ownes' the set
    };
    
      //! List of multi-processor meshsets
    std::list<ParallelSet> parallelSets;
    
      //! Vector indexed by MPI rank, containing the list
      //! of parallel sets that each processor knows about.
    std::map<unsigned,Range> cpuParallelSets;
    typedef std::map<unsigned,Range>::iterator proc_iter;

    //! pcomm controlling parallel nature of mesh
  ParallelComm *myPcomm;

    //! whether this instance allocated (and dtor should delete) the pcomm
  bool pcommAllocated;
};



class WriteHDF5Parallel::MultiProcSetTags::Data
{
  public:
  Data( const std::string& name ) 
   : filterTag(name), dataTag(name), useFilterValue(false) {}
  Data( const std::string& fname, const std::string& dname )
   : filterTag(fname), dataTag(dname), useFilterValue(false) {}
  Data( const std::string& fname, const std::string& dname, int fval )
   : filterTag(fname), dataTag(dname), filterValue(fval), useFilterValue(true) {}
   
  std::string filterTag;
  std::string dataTag;
  int filterValue;
  bool useFilterValue;
  
  ErrorCode get_sets( Interface* moab, Range& sets, Tag& id_tag_out ) const;
};

} // namespace moab

#endif
