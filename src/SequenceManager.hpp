#ifndef SEQUENCE_MANAGER_HPP
#define SEQUENCE_MANAGER_HPP

#include "TypeSequenceManager.hpp"
#include "TagInfo.hpp"
#include <vector>

namespace moab {

class HomCoord;
class TagServer;

class SequenceManager 
{
  public:
    
    ~SequenceManager();
    
      /** Delete all contained data */
    void clear();
    
      /** Find entity sequence containing specified handle.
       *\return MB_SUCCESS or MB_ENTITY_NOT_FOUND
       */
    ErrorCode find( EntityHandle handle, EntitySequence*& sequence_out )
      { 
        return typeData[TYPE_FROM_HANDLE(handle)].find( handle, sequence_out );
      }
    
      /** Find entity sequence containing specified handle.
       *\return MB_SUCCESS or MB_ENTITY_NOT_FOUND
       */
    ErrorCode find( EntityHandle handle, const EntitySequence*& sequence_out ) const
      { 
        return typeData[TYPE_FROM_HANDLE(handle)].find( handle, sequence_out );
      }
    
      /** Get all entities of a given EntityType, return all entities
       *  if type == MBMAXTYPE */
    void get_entities( EntityType type, Range& entities_out ) const
      { 
        if (type == MBMAXTYPE)
          get_entities( entities_out );
        else
          typeData[type].get_entities( entities_out ); 
      }

    void get_entities( Range& entities_out ) const;
    
      /** Get all entities of a given EntityType, return all entities
       *  if type == MBMAXTYPE */
    void get_entities( EntityType type, std::vector<EntityHandle>& entities_out ) const
      { 
        if (type == MBMAXTYPE)
          get_entities( entities_out );
        else
          typeData[type].get_entities( entities_out ); 
      }
    
    void get_entities( std::vector<EntityHandle>& entities_out ) const;
    
      /** Count entities of a given EntityType */
    EntityID get_number_entities( EntityType type ) const
      { return type == MBMAXTYPE ? get_number_entities() : typeData[type].get_number_entities(); }
    
      /** Count entities of a given EntityType */
    EntityID get_number_entities( ) const;
      
      /** Get most recently accessed sequence for a given type */
    const EntitySequence* get_last_accessed_sequence( EntityType type ) const
      { return typeData[type].get_last_accessed(); }
    
      /**\brief Replace subset of existing sequence with new 
       *        sequence (splits existing sequence)
       *
       * Used for converting number of nodes for fixed-connectivity-length
       * elements.  Input sequence must be a non-strict subset of an existing
       * sequence.  Existing sequence will be removed, modified, or split
       * into two prevent it from overlapping the new sequence.
       */
    ErrorCode replace_subsequence( EntitySequence* new_seq, TagServer* ts );
    
      /** Check if passed entity handles are valid */
    ErrorCode check_valid_entities( const Range& entities ) const;
    
      /** Check if passed entity handles are valid */
    ErrorCode check_valid_entities( const EntityHandle entities[],
                                      size_t num_entities ) const;
    
      /** Delete an entity.  Deletes sequence if only contained entity. */
    ErrorCode delete_entity( EntityHandle entity );
    
      /** Delete entities */
    ErrorCode delete_entities( const Range& entities );
    
      /** Allocate a vertex (possibly in an existing sequence) and 
       *  assign it the passed coordinate values.
       */
    ErrorCode create_vertex( const double coords[3],
                               EntityHandle& handle_out );
    
      /** Allocate a element (possibly in an existing sequence) and 
       *  assign it the passed connectivity.
       */
    ErrorCode create_element( EntityType type,
                                const EntityHandle* conn_array,
                                unsigned num_vertices,
                                EntityHandle& handle_out );
    
      /** Allocate an entity set (possibly in an existing sequence) */
    ErrorCode create_mesh_set( unsigned flags,
                                 EntityHandle& handle_out );
      /** Allocate an entity set with the specified handle. 
       *\return MB_ALREADY_ALLOCATED if handle is in use, MB_SUCCESS otherwise.
       */
    ErrorCode allocate_mesh_set( EntityHandle at_this_handle,
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
       *\param first_handle_out First allocated handle.  Allocated handles
       *                    are [first_handle_out, first_handle_out+num_entities-1].
       *\param sequence_out The sequence in which the entities were allocated.
       *                    NOTE: first_handle_out may not be first handle in
       *                    sequence.
       */
    ErrorCode create_entity_sequence( EntityType type,
                                        EntityID num_entities,
                                        int nodes_per_entity,
                                        EntityID start_id_hint,
                                        EntityHandle& first_handle_out,
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
    ErrorCode create_meshset_sequence( EntityID num_sets,
                                         EntityID start_id_hint,
                                         const unsigned* flags,
                                         EntityHandle& first_handle_out,
                                         EntitySequence*& sequence_out );
    
      /**\brief Allocate a block of consecutive mesh sets
       *
       * Alternate form that creates all mesh sets with same flags.
       */
    ErrorCode create_meshset_sequence( EntityID num_sets,
                                         EntityID start_id_hint,
                                         unsigned flags,
                                         EntityHandle& first_handle_out,
                                         EntitySequence*& sequence_out );
    
      /** Create structured mesh */
    ErrorCode create_scd_sequence( int imin, int jmin, int kmin,
                                     int imax, int jmax, int kmax,
                                     EntityType type,
                                     EntityID start_id_hint,
                                     EntityHandle& first_handle_out,
                                     EntitySequence*& sequence_out );
    
      /** Create structured mesh */
    ErrorCode create_scd_sequence( const HomCoord& coord_min,
                                     const HomCoord& coord_max,
                                     EntityType type,
                                     EntityID start_id_hint,
                                     EntityHandle& first_handle_out,
                                     EntitySequence*& sequence_out );

      /** Create swept mesh */
    ErrorCode create_sweep_sequence( int imin, int jmin, int kmin,
				       int imax, int jmax, int kmax,
				       int* Cq,
				       EntityType type,
				       EntityID start_id_hint,
				       EntityHandle& first_handle_out,
				       EntitySequence*& sequence_out );
    
      /** Create swept mesh */
    ErrorCode create_sweep_sequence( const HomCoord& coord_min,
				       const HomCoord& coord_max,
				       int* Cq,
				       EntityType type,
				       EntityID start_id_hint,
				       EntityHandle& first_handle_out,
				       EntitySequence*& sequence_out );

    /** Add a structured vertex sequence to this structured element sequence;
     * see comments in ScdElementData */
  ErrorCode add_vsequence(EntitySequence *vert_seq,
                            EntitySequence *elem_seq,
                            const HomCoord &p1, const HomCoord &q1,
                            const HomCoord &p2, const HomCoord &q2,
                            const HomCoord &p3, const HomCoord &q3,
                            bool bb_input = false,
                            const HomCoord *bb_min = NULL,
                            const HomCoord *bb_max = NULL);

      /** Get data for a specific EntityType */
    TypeSequenceManager& entity_map( EntityType type )
      { return typeData[type]; }
    
      /** Get data for a specific EntityType */
    const TypeSequenceManager& entity_map( EntityType type ) const
      { return typeData[type]; }
    
    void get_memory_use( unsigned long& total_entity_storage,
                         unsigned long& total_storage ) const;
                         
    void get_memory_use( EntityType type,
                         unsigned long& total_entity_storage,
                         unsigned long& total_storage ) const;
    
    void get_memory_use( const Range& entities,
                         unsigned long& total_entity_storage,
                         unsigned long& total_amortized_storage ) const;
    
  
  
    /* Dense Tag Functions */
    
      /** Release all dense tag data storage */
    void reset_tag_data();
    
      /** Allocate a tag ID
       *\param tag_id   The ID to allocate/reserve
       *\param tag_size The size of the tag value for each entity
       */
    ErrorCode reserve_tag_id( int tag_size, TagId tag_id );
    
      /** Release a reserved tag ID any any associated storage */
    ErrorCode release_tag( TagId tag_id );
    
      /** If tag data is allocated for the specified entity,
       *  change it to the passed default value.  Otherwise 
       *  do nothing. 
       */
    ErrorCode remove_tag_data( TagId tag_id, 
                                 EntityHandle handle,
                                 const void* default_tag_value,
                                 int default_value_size = 0 );
                                 
      /** Set fixed-length tag values.
       *\NOTE Default value must be given because it is often
       *      necessary to allocate storage for additional entities
       *
       *\NOTE Will fail for variable-length tag data.
       */
    ErrorCode set_tag_data( TagId tag_id,
                              const EntityHandle* handles,
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
       *\param one_value   If true, tag on all entities is set to the (same)
       *                   first passed tag value.
       */
    ErrorCode set_tag_data( TagId tag_id,
                              const EntityHandle* handles,
                              int num_handles,
                              void const* const* values,
                              const int* lengths,
                              const void* default_value,
                              bool one_value = false );
                              
      /** Set fixed-length tag value for an Range of entities
       *\NOTE Default value must be given because it is often
       *      necessary to allocate storage for other entities
       *\NOTE Will fail for variable-length tag data
       */
    ErrorCode set_tag_data( TagId tag_id,
                              const Range& handles,
                              const void* values,
                              const void* default_value );
                              
      /** Set tag data for an Range of entities.
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
       *\param one_value   If true, tag on all entities is set to the (same)
       *                   first passed tag value.
       */
     ErrorCode set_tag_data( TagId tag_id,
                               const Range& handles,
                               void const* const* values,
                               const int* lengths,
                               const void* default_value,
                               bool one_value = false );

  /**\brief Access tag data via direct pointer into contiguous blocks
   *
   * Iteratively obtain direct access to contiguous blocks of tag
   * storage.  This function cannot be used with bit tags because
   * of the compressed bit storage.  This function cannot be used
   * with variable length tags because it does not provide a mechanism
   * to determine the length of the value for each entity.  This
   * function may be used with sparse tags, but if it is used, it
   * will return data for a single entity at a time.  
   *
   *\param tag_id   The ID of the tag for which to access data
   *\param iter     As input, the first entity for which to return
   *                data.  As output, one past the last entity for
   *                which data was returned.
   *\param end      One past the last entity for which data is desired
   *\param data_ptr Output: pointer to tag storage.
   *  
   *\Note If this function is called for entities for which no tag value
   *      has been set, but for which a default value exists, it will 
   *      force the allocation of explicit storage for each such entity
   *      even though MOAB would normally not explicitly store tag values
   *      for such entities.
   */
    ErrorCode tag_iterate( TagId tag_id,
                           Range::iterator& iter,
                           const Range::iterator& end,
                           void*& data_ptr,
                           const void* default_value );

      /** Get fixed-length tag values for array of entities
       *\NOTE Will fail with MB_VARIABLE_DATA_LENGTH if called
       *      for variable-length tag.
       */
    ErrorCode get_tag_data( TagId tag_id,
                              const EntityHandle* handles,
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
    ErrorCode get_tag_data( TagId tag_id,
                              const EntityHandle* handles,
                              int num_handles,
                              const void** tag_ptrs,
                              int* lengths,
                              const void* default_value,
                              int default_value_length ) const;
                              
      /** Get fixed-length tag value for an Range of entities
       *\NOTE Will fail with MB_VARIABLE_DATA_LENGTH if called
       *      for variable-length tag.
       */
    ErrorCode get_tag_data( TagId tag_id,
                              const Range& handles,
                              void* values,
                              const void* default_value ) const;
                              
      /** Get pointers to tag data for an Range of entities.
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
    ErrorCode get_tag_data( TagId tag_id,
                              const Range& handles,
                              const void** values,
                              int* lengths,
                              const void* default_value,
                              int defaut_value_length ) const;
  
      /** Get all tags for which values (possibly the default value)
       *  have been allocated for a given entity.
       *\NOTE For variable-length data, will only return tag if
       *      data length is greater than zero.
       */
    ErrorCode get_entity_tags( EntityHandle entity,
                                 std::vector<Tag>& tags_out ) const;

      /** Get all entities for which storage for a specific tag has
       *  been allocated.
       *\NOTE For variable-length data, will only return entities for
       *      which data length is greater than zero.
       */
    ErrorCode get_tagged_entities( TagId tag_id, 
                                     EntityType type,
                                     Range& entities_out ) const;

      /** Count all entities for which storage for a specific tag has
       *  been allocated.
       *\NOTE For variable-length data, will only count entities for
       *      which data length is greater than zero.
       */
    ErrorCode count_tagged_entities( TagId tag, 
                                       EntityType type, 
                                       int& result ) const;

      /** Get entities by type and tag value (intersection) */
    ErrorCode get_entities_with_tag_value( TagId id,
                                             const TagInfo& tag_info,
                                             EntityType type,
                                             Range& entities_out,
                                             const void* value,
                                             int value_size ) const;
                                             
      /** Get subset of entities with a type and a tag value (intersection)
       *\param range Intersect result with this range
       *\param id    The ID of the tag to match
       *\param type  Only check entities of this type 
       *\param entities_out Result (subset of input 'range')
       *\param value The tag value
       */
    ErrorCode get_entities_with_tag_value( const Range& range,
                                             TagId id,
                                             const TagInfo& tag_info,
                                             EntityType type,
                                             Range& entities_out,
                                             const void* value,
                                             int value_size ) const;
    
    ErrorCode get_tag_memory_use( TagId id, 
                                    unsigned long& total, 
                                    unsigned long& per_entity ) const;
    
      /**\brief Get default size of POLYGON and POLYHEDRON SequenceData */
    static EntityID default_poly_sequence_size( int entity_connectivity_length );
    
      /**\brief Size to allocate for new SquenceData */
    EntityID new_sequence_size( EntityHandle start_handle, 
                                  EntityID reqested_size,
                                  EntityID default_size ) const;
    
  private:
   
    /**\brief Utility function for allocate_mesh_set (and similar)
     *
     * Given a block of available handles, determine the non-strict
     * subset at which to create a new EntitySequence.
     */
    void trim_sequence_block( EntityHandle start_handle,
                              EntityHandle& end_handle_in_out,
                              unsigned maximum_sequence_size );
  
  
      /**\brief Get range of handles in which to create an entity sequence
       *
       * Get range of handles in whcih to place a new entity sequence.
       *\param type              The EntityType for the contents of the sequence
       *\param entity_count      The number of entities in the range
       *\param values_per_entity Vertices per element, zero for other types
       *\param start_id_hint     Preferred id of first handle
       *\param processor_rank    MPI processor ID
       *\param data_out  Output: Either NULL or an existing SequenceData
       *                         with a sufficiently large block to accomodate 
       *                         the handle range.
       *\return zero if no available handle range, start handle otherwise.
       */
    EntityHandle sequence_start_handle( EntityType type,
                                          EntityID entity_count,
                                          int values_per_entity,
                                          EntityID start_id_hint,
                                          SequenceData*& data_out,
                                          EntityID &data_size );
  
    TypeSequenceManager typeData[MBMAXTYPE];
    
    std::vector<int> tagSizes;
};

} // namespace moab

#endif
