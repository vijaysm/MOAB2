#ifndef SEQUENCE_MANAGER_HPP
#define SEQUENCE_MANAGER_HPP

#include "TypeSequenceManager.hpp"
#include "MBHandleUtils.hpp"
#include "TagInfo.hpp"

class HomCoord;
class TagServer;

class SequenceManager 
{
  public:
  
    SequenceManager( const MBHandleUtils& handle_utils ) 
      : handleUtils(handle_utils)
      {}
    
    ~SequenceManager();
    
      /** Delete all contained data */
    void clear();
    
      /** Find entity sequence containing specified handle.
       *\return MB_SUCCESS or MB_ENTITY_NOT_FOUND
       */
    MBErrorCode find( MBEntityHandle handle, EntitySequence*& sequence_out )
      { 
        return typeData[TYPE_FROM_HANDLE(handle)].find( handle, sequence_out );
      }
    
      /** Find entity sequence containing specified handle.
       *\return MB_SUCCESS or MB_ENTITY_NOT_FOUND
       */
    MBErrorCode find( MBEntityHandle handle, const EntitySequence*& sequence_out ) const
      { 
        return typeData[TYPE_FROM_HANDLE(handle)].find( handle, sequence_out );
      }
    
      /** Get all entities of a given MBEntityType */
   void get_entities( MBEntityType type, MBRange& entities_out ) const
      { typeData[type].get_entities( entities_out ); }
    
      /** Count entities of a given MBEntityType */
    MBEntityID get_number_entities( MBEntityType type ) const
      { return typeData[type].get_number_entities(); }
      
      /** Get most recently accessed sequence for a given type */
    const EntitySequence* get_last_accessed_sequence( MBEntityType type ) const
      { return typeData[type].get_last_accessed(); }
    
      /**\brief Replace subset of existing sequence with new 
       *        sequence (splits existing sequence)
       *
       * Used for converting number of nodes for fixed-connectivity-length
       * elements.  Input sequence must be a non-strict subset of an existing
       * sequence.  Existing sequence will be removed, modified, or split
       * into two prevent it from overlapping the new sequence.
       */
    MBErrorCode replace_subsequence( EntitySequence* new_seq, TagServer* ts );
    
      /** Check if passed entity handles are valid */
    MBErrorCode check_valid_entities( const MBRange& entities ) const;
    
      /** Delete an entity.  Deletes sequence if only contained entity. */
    MBErrorCode delete_entity( MBEntityHandle entity );
    
      /** Delete entities */
    MBErrorCode delete_entities( const MBRange& entities );
    
      /** Allocate a vertex (possibly in an existing sequence) and 
       *  assign it the passed coordinate values.
       */
    MBErrorCode create_vertex( unsigned processor_id,
                               const double coords[3],
                               MBEntityHandle& handle_out );
    
      /** Allocate a element (possibly in an existing sequence) and 
       *  assign it the passed connectivity.
       */
    MBErrorCode create_element( MBEntityType type,
                                unsigned processor_id,
                                const MBEntityHandle* conn_array,
                                unsigned num_vertices,
                                MBEntityHandle& handle_out );
    
      /** Allocate an entity set (possibly in an existing sequence) */
    MBErrorCode create_mesh_set( unsigned processor_id,
                                 unsigned flags,
                                 MBEntityHandle& handle_out );
      /** Allocate an entity set with the specified handle. 
       *\return MB_ALREADY_ALLOCATED if handle is in use, MB_SUCCESS otherwise.
       */
    MBErrorCode allocate_mesh_set( MBEntityHandle at_this_handle,
                                   unsigned flags );
    
      /**\brief Allocate a block of consecutive entity handles
       *
       * Allocate a block of consecutive entity handles.  Handles
       * may be appended or prepended to an existing entity sequence.
       *\param type The type of of entity for which to allocate handles
       *\param num_entities Number of entities to allocate
       *\param nodes_per_entity Number of nodes in connectivity for elements,
       *                    ignored MBVERTEX, MBPOLYGON, MBPOLYHEDRON, and
       *                    MBENTITYSET types.
       *\param start_id_hint Preferred ID portion for first handle.  
       *                    May be ignored if not available.
       *\param processor_id Processor ID to embed in handles
       *\param first_handle_out First allocated handle.  Allocated handles
       *                    are [first_handle_out, first_handle_out+num_entities-1].
       *\param sequence_out The sequence in which the entities were allocated.
       *                    NOTE: first_handle_out may not be first handle in
       *                    sequence.
       */
    MBErrorCode create_entity_sequence( MBEntityType type,
                                        MBEntityID num_entities,
                                        int nodes_per_entity,
                                        MBEntityID start_id_hint,
                                        int processor_id,
                                        MBEntityHandle& first_handle_out,
                                        EntitySequence*& sequence_out );
    
      /**\brief Allocate a block of consecutive mesh sets
       *
       * Allocate a block of consecutive entity handles.  Handles
       * may be appended or prepended to an existing entity sequence.
       *\param type The type of of entity for which to allocate handles
       *\param num_sets     Number of entities to allocate
       *\param start_id_hint Preferred ID portion for first handle.  
       *                    May be ignored if not available.
       *\param processor_id Processor ID to embed in handles
       *\param flags        Array of length 'num_sets' containing entity set
       *                    creating flags.
       *\param first_handle_out First allocated handle.  Allocated handles
       *                    are [first_handle_out, first_handle_out+num_entities-1].
       *\param sequence_out The sequence in which the entities were allocated.
       *                    NOTE: first_handle_out may not be first handle in
       *                    sequence.
       */
    MBErrorCode create_meshset_sequence( MBEntityID num_sets,
                                         MBEntityID start_id_hint,
                                         int processor_id,
                                         const unsigned* flags,
                                         MBEntityHandle& first_handle_out,
                                         EntitySequence*& sequence_out );
    
      /**\brief Allocate a block of consecutive mesh sets
       *
       * Alternate form that creates all mesh sets with same flags.
       */
    MBErrorCode create_meshset_sequence( MBEntityID num_sets,
                                         MBEntityID start_id_hint,
                                         int processor_id,
                                         unsigned flags,
                                         MBEntityHandle& first_handle_out,
                                         EntitySequence*& sequence_out );
    
      /** Create structured mesh */
    MBErrorCode create_scd_sequence( int imin, int jmin, int kmin,
                                     int imax, int jmax, int kmax,
                                     MBEntityType type,
                                     MBEntityID start_id_hint,
                                     int processor_id,
                                     MBEntityHandle& first_handle_out,
                                     EntitySequence*& sequence_out );
    
      /** Create structured mesh */
    MBErrorCode create_scd_sequence( const HomCoord& coord_min,
                                     const HomCoord& coord_max,
                                     MBEntityType type,
                                     MBEntityID start_id_hint,
                                     int processor_id,
                                     MBEntityHandle& first_handle_out,
                                     EntitySequence*& sequence_out );
                                     
      /** Get data for a specific MBEntityType */
    TypeSequenceManager& entity_map( MBEntityType type )
      { return typeData[type]; }
    
      /** Get data for a specific MBEntityType */
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
    
      /** Release all dense tag data storage */
    void reset_tag_data();
    
      /** Allocate a tag ID
       *\param tag_id   The ID to allocate/reserve
       *\param tag_size The size of the tag value for each entity
       */
    MBErrorCode reserve_tag_id( int tag_size, MBTagId tag_id );
    
      /** Release a reserved tag ID any any associated storage */
    MBErrorCode release_tag( MBTagId tag_id );
    
      /** If tag data is allocated for the specified entity,
       *  change it to the passed default value.  Otherwise 
       *  do nothing. 
       */
    MBErrorCode remove_tag_data( MBTagId tag_id, 
                                 MBEntityHandle handle,
                                 const void* default_tag_value,
                                 int default_value_size = 0 );
                                 
      /** Set fixed-length tag values.
       *\NOTE Default value must be given because it is often
       *      necessary to allocate storage for additional entities
       *
       *\NOTE Will fail for variable-length tag data.
       */
    MBErrorCode set_tag_data( MBTagId tag_id,
                              const MBEntityHandle* handles,
                              int num_handles,
                              const void* values,
                              const void* default_value );
                              
      /** Set tag values for array of entity handles
       *\param tag_id      The tag.
       *\param handles     Array of entity handles.
       *\param num_handles Length of 'handles' array.
       *\param values      Array of pointers to tag values, one pointer for each handle
       *\param lengths     Length of each tag value.  Ignored for fixed-length tags.
       *\param default_value Used to initialize any additional tag storage.  Ignored
       *                   for variable-length tags.
       */
    MBErrorCode set_tag_data( MBTagId tag_id,
                              const MBEntityHandle* handles,
                              int num_handles,
                              void const* const* values,
                              const int* lengths,
                              const void* default_value );
                              
      /** Set fixed-length tag value for an MBRange of entities
       *\NOTE Default value must be given because it is often
       *      necessary to allocate storage for other entities
       *\NOTE Will fail for variable-length tag data
       */
    MBErrorCode set_tag_data( MBTagId tag_id,
                              const MBRange& handles,
                              const void* values,
                              const void* default_value );
                              
      /** Set tag data for an MBRange of entities.
       *
       *\param tag_id  The tag
       *\param handles The entities
       *\param values  An array of pointers, one per entity, pointing to the
       *               tag value for the corresponding entity.
       *\param lengths An array of integers, one per entity, indicating the
       *               length of the tag value for each entity.  Ingored
       *               for fixed-length tags.
       *\param default_value The default value for the tag.  Ignored for
       *               variable-length tags.
       */
     MBErrorCode set_tag_data( MBTagId tag_id,
                               const MBRange& handles,
                               void const* const* values,
                               const int* lengths,
                               const void* default_value );

      /** Get fixed-length tag values for array of entities
       *\NOTE Will fail with MB_VARIABLE_DATA_LENGTH if called
       *      for variable-length tag.
       */
    MBErrorCode get_tag_data( MBTagId tag_id,
                              const MBEntityHandle* handles,
                              int num_handles,
                              void* values, 
                              const void* default_value ) const;
                              
      /** Get pointers to tag data for array of entities
       *\param tag_id      The Tag.
       *\param handles     Array of entity handles.
       *\param num_handles Length of 'handles' array.
       *\param tag_ptrs    Pointers to tag values, one pointer for each input handle.
       *\param lengths     Length of each tag value.  Ignored for fixed-length tags.
       *\param default_value Pointer to default value for tag, or NULL if none.
       *\param default_value_length  Length of default tag value.  Ingored for
       *                   fixed-length tags.
       */
    MBErrorCode get_tag_data( MBTagId tag_id,
                              const MBEntityHandle* handles,
                              int num_handles,
                              const void** tag_ptrs,
                              int* lengths,
                              const void* default_value,
                              int default_value_length ) const;
                              
      /** Get fixed-length tag value for an MBRange of entities
       *\NOTE Will fail with MB_VARIABLE_DATA_LENGTH if called
       *      for variable-length tag.
       */
    MBErrorCode get_tag_data( MBTagId tag_id,
                              const MBRange& handles,
                              void* values,
                              const void* default_value ) const;
                              
      /** Get pointers to tag data for an MBRange of entities.
       *
       *\param tag_id   The tag.
       *\param handles  The entities.
       *\param values   Array of pointers of type 'const void*'.  Array
       *                must be the same length as the size of 'entities'.
       *                The array will be populated with pointers to the
       *                internal storage of the tag data for each entity.
       *\param lengths  Array of integers.  Will be populated with the 
       *                length of the tag value for each entity.  Argument
       *                is optional for fixed-length tags.
       *\param default_value The default value for the tag.
       *\param default_value_length The length of the default tag value.
       */
    MBErrorCode get_tag_data( MBTagId tag_id,
                              const MBRange& handles,
                              const void** values,
                              int* lengths,
                              const void* default_value,
                              int defaut_value_length ) const;
  
      /** Get all tags for which values (possibly the default value)
       *  have been allocated for a given entity.
       *\NOTE For variable-length data, will only return tag if
       *      data length is greater than zero.
       */
    MBErrorCode get_entity_tags( MBEntityHandle entity,
                                 std::vector<MBTag>& tags_out ) const;

      /** Get all entities for which storage for a specific tag has
       *  been allocated.
       *\NOTE For variable-length data, will only return entities for
       *      which data length is greater than zero.
       */
    MBErrorCode get_tagged_entities( MBTagId tag_id, 
                                     MBEntityType type,
                                     MBRange& entities_out ) const;

      /** Count all entities for which storage for a specific tag has
       *  been allocated.
       *\NOTE For variable-length data, will only count entities for
       *      which data length is greater than zero.
       */
    MBErrorCode count_tagged_entities( MBTagId tag, 
                                       MBEntityType type, 
                                       int& result ) const;

      /** Get entities by type and tag value (intersection) */
    MBErrorCode get_entities_with_tag_value( MBTagId id,
                                             const TagInfo& tag_info,
                                             MBEntityType type,
                                             MBRange& entities_out,
                                             const void* value,
                                             int value_size ) const;
                                             
      /** Get subset of entities with a type and a tag value (intersection)
       *\param range Intersect result with this range
       *\param id    The ID of the tag to match
       *\param type  Only check entities of this type 
       *\param entities_out Result (subset of input 'range')
       *\param value The tag value
       */
    MBErrorCode get_entities_with_tag_value( const MBRange& range,
                                             MBTagId id,
                                             const TagInfo& tag_info,
                                             MBEntityType type,
                                             MBRange& entities_out,
                                             const void* value,
                                             int value_size ) const;
    
    MBErrorCode get_tag_memory_use( MBTagId id, 
                                    unsigned long& total, 
                                    unsigned long& per_entity ) const;
    
      /**\brief Get default size of POLYGON and POLYHEDRON SequenceData */
    static MBEntityID default_poly_sequence_size( int entity_connectivity_length );
    
      /**\brief Size to allocate for new SquenceData */
    MBEntityID new_sequence_size( MBEntityHandle start_handle, 
                                  MBEntityID reqested_size,
                                  MBEntityID default_size ) const;
    
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
    
    std::vector<int> tagSizes;
};

#endif
