#ifndef SEQUENCE_MANAGER_HPP
#define SEQUENCE_MANAGER_HPP

#include "TypeSequenceManager.hpp"
#include "MBHandleUtils.hpp"

class HomCoord;
class TagServer;

class SequenceManager 
{
  public:
  
    SequenceManager( const MBHandleUtils& handle_utils ) 
      : handleUtils(handle_utils)
      {}
      
      /** Delete all contained data */
    void clear();
    
    MBErrorCode find( MBEntityHandle handle, EntitySequence*& sequence_out )
      { 
        return typeData[TYPE_FROM_HANDLE(handle)].find( handle, sequence_out );
      }
    
    MBErrorCode find( MBEntityHandle handle, const EntitySequence*& sequence_out ) const
      { 
        return typeData[TYPE_FROM_HANDLE(handle)].find( handle, sequence_out );
      }
    
    void get_entities( MBEntityType type, MBRange& entities_out ) const
      { typeData[type].get_entities( entities_out ); }
    
    MBEntityID get_number_entities( MBEntityType type ) const
      { return typeData[type].get_number_entities(); }
    
    MBErrorCode replace_subsequence( EntitySequence* new_seq, TagServer* ts );
    
    MBErrorCode check_valid_entities( const MBRange& entities ) const;
    
    MBErrorCode delete_entity( MBEntityHandle entity );
    
    MBErrorCode delete_entities( const MBRange& entities );
    
    MBErrorCode create_vertex( unsigned processor_id,
                               const double coords[3],
                               MBEntityHandle& handle_out );
    
    MBErrorCode create_element( MBEntityType type,
                                unsigned processor_id,
                                const MBEntityHandle* conn_array,
                                unsigned num_vertices,
                                MBEntityHandle& handle_out );
    
    MBErrorCode create_mesh_set( unsigned processor_id,
                                 unsigned flags,
                                 MBEntityHandle& handle_out );
    
    MBErrorCode allocate_mesh_set( MBEntityHandle at_this_handle,
                                   unsigned flags );
    
    MBErrorCode create_entity_sequence( MBEntityType type,
                                        MBEntityID num_entities,
                                        int nodes_per_entity,
                                        MBEntityID start_id_hint,
                                        int processor_id,
                                        MBEntityHandle& first_handle_out,
                                        EntitySequence*& sequence_out );
    
    MBErrorCode create_meshset_sequence( MBEntityID num_sets,
                                         MBEntityID start_id_hint,
                                         int processor_id,
                                         const unsigned* flags,
                                         MBEntityHandle& first_handle_out,
                                         EntitySequence*& sequence_out );
    
    MBErrorCode create_meshset_sequence( MBEntityID num_sets,
                                         MBEntityID start_id_hint,
                                         int processor_id,
                                         unsigned flags,
                                         MBEntityHandle& first_handle_out,
                                         EntitySequence*& sequence_out );
    
    MBErrorCode create_scd_sequence( int imin, int jmin, int kmin,
                                     int imax, int jmax, int kmax,
                                     MBEntityType type,
                                     MBEntityID start_id_hint,
                                     int processor_id,
                                     MBEntityHandle& first_handle_out,
                                     EntitySequence*& sequence_out );
    
    MBErrorCode create_scd_sequence( const HomCoord& coord_min,
                                     const HomCoord& coord_max,
                                     MBEntityType type,
                                     MBEntityID start_id_hint,
                                     int processor_id,
                                     MBEntityHandle& first_handle_out,
                                     EntitySequence*& sequence_out );
                                     
    
    TypeSequenceManager& entity_map( MBEntityType type )
      { return typeData[type]; }
    
    const TypeSequenceManager& entity_map( MBEntityType type ) const
      { return typeData[type]; }
    
    void get_memory_use( unsigned long& total_entity_storage,
                         unsigned long& total_storage ) const;
                         
    void get_memory_use( MBEntityType type,
                         unsigned long& total_entity_storage,
                         unsigned long& total_storage ) const;
    
    void get_memory_use( const MBRange& entities,
                         unsigned long& total_entity_storage,
                         unsigned long& total_amortized_storage ) const;
    
  
  
    /* Dense Tag Functions */
    
    void reset_tag_data();
    
    MBErrorCode reserve_tag_id( unsigned tag_size, MBTagId tag_id );
    MBErrorCode release_tag( MBTagId tag_id );
    
    MBErrorCode remove_tag_data( MBTagId tag_id, 
                                 MBEntityHandle handle,
                                 const void* default_tag_value );
    MBErrorCode set_tag_data( MBTagId tag_id,
                              MBEntityHandle handle,
                              const void* value,
                              const void* default_value );
    MBErrorCode set_tag_data( MBTagId tag_id,
                              const MBRange& handles,
                              const void* values,
                              const void* default_value );
    MBErrorCode get_tag_data( MBTagId tag_id,
                              MBEntityHandle handle,
                              void* value ) const;
    MBErrorCode get_tag_data( MBTagId tag_id,
                              const MBRange& handles,
                              void* values,
                              const void* default_value ) const;
  
    MBErrorCode get_entity_tags( MBEntityHandle entity,
                                 std::vector<MBTag>& tags_out ) const;

    MBErrorCode get_tagged_entities( MBTagId tag_id, 
                                     MBEntityType type,
                                     MBRange& entities_out ) const;
    MBErrorCode count_tagged_entities( MBTagId tag, 
                                       MBEntityType type, 
                                       int& result ) const;

    MBErrorCode get_entities_with_tag_value( MBTagId id,
                                             MBEntityType type,
                                             MBRange& entities_out,
                                             const void* value ) const;
    MBErrorCode get_entities_with_tag_value( const MBRange& range,
                                             MBTagId id,
                                             MBEntityType type,
                                             MBRange& entities_out,
                                             const void* value ) const;
    
    MBErrorCode get_tag_memory_use( MBTagId id, 
                                    unsigned long& total, 
                                    unsigned long& per_entity ) const;
    
    
  private:
  
    /**\brief Utility function for allocate_mesh_set (and similar)
     *
     * Given a block of available handles, determine the non-strict
     * subset at which to create a new EntitySequence.
     */
    void trim_sequence_block( MBEntityHandle start_handle,
                              MBEntityHandle& end_handle_in_out,
                              unsigned maximum_sequence_size );
  
  
      /**\brief Get range of handles in which to create an entity sequence
       *
       * Get range of handles in whcih to place a new entity sequence.
       *\param type              The MBEntityType for the contents of the sequence
       *\param entity_count      The number of entities in the range
       *\param values_per_entity Vertices per element, zero for other types
       *\param start_id_hint     Preferred id of first handle
       *\param processor_rank    MPI processor ID
       *\param data_out  Output: Either NULL or an existing SequenceData
       *                         with a sufficiently large block to accomodate 
       *                         the handle range.
       *\return zero if no available handle range, start handle otherwise.
       */
    MBEntityHandle sequence_start_handle(   MBEntityType type,
                                              MBEntityID entity_count,
                                                     int values_per_entity,
                                              MBEntityID start_id_hint,
                                                     int processor_rank,
                                          SequenceData*& data_out );
  
    const MBHandleUtils handleUtils;
    TypeSequenceManager typeData[MBMAXTYPE];
    
    std::vector<unsigned> tagSizes;
};

#endif
