#include "SequenceManager.hpp"
#include "VertexSequence.hpp"
#include "UnstructuredElemSeq.hpp"
#include "ScdVertexData.hpp"
#include "MeshSetSequence.hpp"
#include "SweptElementSeq.hpp"
#include "StructuredElementSeq.hpp"
#include "moab/HomXform.hpp"
#include "PolyElementSeq.hpp"
#include "SysUtil.hpp"
#include "TagCompare.hpp"

#include <assert.h>
#include <new>
#include <algorithm>

#ifndef NDEBUG
#include <iostream>
#endif 

namespace moab {

const EntityID DEFAULT_VERTEX_SEQUENCE_SIZE = 4096;
const EntityID DEFAULT_ELEMENT_SEQUENCE_SIZE = DEFAULT_VERTEX_SEQUENCE_SIZE;
const EntityID DEFAULT_POLY_SEQUENCE_SIZE = 4 * DEFAULT_ELEMENT_SEQUENCE_SIZE;
const EntityID DEFAULT_MESHSET_SEQUENCE_SIZE = DEFAULT_VERTEX_SEQUENCE_SIZE;

static inline
const unsigned char* tag_array( const EntitySequence* seq, 
                                EntityHandle h,
                                int tag_id, 
                                int tag_size )
{
  const void* mem = seq->data()->get_tag_data(tag_id);
  return mem ? reinterpret_cast<const unsigned char*>(mem)
    + tag_size * (h - seq->data()->start_handle()) : 0;
}

static inline
unsigned char* tag_array( EntitySequence* seq, 
                          EntityHandle h, 
                          int tag_id, 
                          int tag_size )
{
  void* mem = seq->data()->get_tag_data(tag_id);
  return mem ? reinterpret_cast<unsigned char*>(mem)
    + tag_size * (h - seq->data()->start_handle()) : 0;
}

static inline
unsigned char* make_tag( EntitySequence* seq, 
                         EntityHandle h, 
                         int tag_id, 
                         int tag_size,
                         const void* default_value )
{
  void* mem = seq->data()->get_tag_data(tag_id);
  if (!mem)
    mem = seq->data()->create_tag_data( tag_id, tag_size, default_value );
  return reinterpret_cast<unsigned char*>(mem)
    + tag_size * (h - seq->data()->start_handle());
}

static inline
const VarLenTag* vtag_array( const EntitySequence* seq, 
                             EntityHandle h, 
                             int tag_id )
{
  const void* mem = seq->data()->get_tag_data(tag_id);
  return mem ? reinterpret_cast<const VarLenTag*>(mem) + h - seq->data()->start_handle() : 0;
}

static inline
VarLenTag* vtag_array( EntitySequence* seq, 
                       EntityHandle h, 
                       int tag_id )
{
  void* mem = seq->data()->get_tag_data(tag_id);
  return mem ? reinterpret_cast<VarLenTag*>(mem) + h - seq->data()->start_handle() : 0;
}

static inline
VarLenTag* make_vtag( EntitySequence* seq, 
                      EntityHandle h, 
                      int tag_id )
{
  void* mem = seq->data()->get_tag_data(tag_id);
  if (!mem)
    mem = seq->data()->create_tag_data( tag_id, sizeof(VarLenTag), 0 );
  return reinterpret_cast<VarLenTag*>(mem) + h - seq->data()->start_handle();
}

EntityID SequenceManager::default_poly_sequence_size( int conn_len )
  {  return std::max( DEFAULT_POLY_SEQUENCE_SIZE / conn_len, (EntityID)1 ); }

SequenceManager::~SequenceManager()
{
    // release variable-length tag data
  for (unsigned i = 0; i < tagSizes.size(); ++i)
    if (tagSizes[i] == MB_VARIABLE_LENGTH)
      release_tag( i );
}

void SequenceManager::clear()
{
    // release variable-length tag data
  for (unsigned i = 0; i < tagSizes.size(); ++i)
    if (tagSizes[i] == MB_VARIABLE_LENGTH)
      release_tag( i );

    // destroy all TypeSequenceManager instances
  for (EntityType t = MBVERTEX; t < MBMAXTYPE; ++t)
    typeData[t].~TypeSequenceManager();
    
    // now re-create TypeSequenceManager instances
  for (EntityType t = MBVERTEX; t < MBMAXTYPE; ++t)
    new (typeData+t) TypeSequenceManager();
}  

void SequenceManager::get_entities( Range& entities_out ) const
{
  for (EntityType t = MBENTITYSET; t >= MBVERTEX; --t)
    typeData[t].get_entities( entities_out );
}

void SequenceManager::get_entities( std::vector<EntityHandle>& entities_out ) const
{
  for (EntityType t = MBVERTEX; t != MBMAXTYPE; ++t)
    typeData[t].get_entities( entities_out );
}

EntityID SequenceManager::get_number_entities( ) const
{
  EntityID sum = 0;
  for (EntityType t = MBVERTEX; t != MBMAXTYPE; ++t)
    sum += typeData[t].get_number_entities();
  return sum;
}



ErrorCode SequenceManager::check_valid_entities( const Range& entities ) const
{
  ErrorCode rval;
  Range::const_pair_iterator i;
  for (i = entities.const_pair_begin(); i != entities.const_pair_end(); ++i) {
    const EntityType type1 = TYPE_FROM_HANDLE(i->first);
    const EntityType type2 = TYPE_FROM_HANDLE(i->second);
    if (type1 == type2) {
      rval = typeData[type1].check_valid_handles( i->first, i->second );
      if (MB_SUCCESS != rval)
        return rval;
    }
    else {
      int junk;
      EntityHandle split = CREATE_HANDLE( type2, 0, junk );
      rval = typeData[type1].check_valid_handles( i->first, split-1 );
      if (MB_SUCCESS != rval)
        return rval;
      rval = typeData[type2].check_valid_handles( split, i->second );
      if (MB_SUCCESS != rval)
        return rval;
    }
  }
  return MB_SUCCESS;
}

ErrorCode SequenceManager::check_valid_entities( const EntityHandle* entities,
                                                   size_t num_entities ) const
{
  ErrorCode rval = MB_SUCCESS;
  const EntitySequence* ptr = 0;
  
  const EntityHandle* const end = entities + num_entities;
  for (; entities < end; ++entities) {
    rval = find(*entities, ptr);
    if (MB_SUCCESS != rval)
      break;
  }
  
  return rval;
}

ErrorCode SequenceManager::delete_entity( EntityHandle entity )
{
  return typeData[TYPE_FROM_HANDLE(entity)].erase( entity );
}

ErrorCode SequenceManager::delete_entities( const Range& entities )
{
  ErrorCode rval = check_valid_entities( entities );
  if (MB_SUCCESS != rval)
    return rval;
  
  ErrorCode result = MB_SUCCESS;
  Range::const_pair_iterator i;
  for (i = entities.const_pair_begin(); i != entities.const_pair_end(); ++i) {
    const EntityType type1 = TYPE_FROM_HANDLE(i->first);
    const EntityType type2 = TYPE_FROM_HANDLE(i->second);
    if (type1 == type2) {
      rval = typeData[type1].erase( i->first, i->second );
      if (MB_SUCCESS != rval)
        return result = rval;
    }
    else {
      int junk;
      EntityHandle split = CREATE_HANDLE( type2, 0, junk );
      rval = typeData[type1].erase( i->first, split-1 );
      if (MB_SUCCESS != rval)
        return result = rval;
      rval = typeData[type2].erase( split, i->second );
      if (MB_SUCCESS != rval)
        return result = rval;
    }
  }
  return result;
}
  
ErrorCode SequenceManager::create_vertex( const double coords[3],
                                            EntityHandle& handle )
{
  const EntityHandle start = CREATE_HANDLE( MBVERTEX, MB_START_ID );
  const EntityHandle   end = CREATE_HANDLE( MBVERTEX,   MB_END_ID );
  bool append;
  TypeSequenceManager::iterator seq = typeData[MBVERTEX].find_free_handle( start, end, append );
  VertexSequence* vseq;
  
  if (seq == typeData[MBVERTEX].end()) {
    SequenceData* seq_data = 0;
    EntityID seq_data_size = 0;
    handle = typeData[MBVERTEX].find_free_sequence( DEFAULT_VERTEX_SEQUENCE_SIZE, start, end, seq_data, seq_data_size );
    if (!handle) 
      return MB_FAILURE;
    
    if (seq_data) 
      vseq = new VertexSequence( handle, 1, seq_data );
    else
      vseq = new VertexSequence( handle, 1, DEFAULT_VERTEX_SEQUENCE_SIZE );
      
    ErrorCode rval = typeData[MBVERTEX].insert_sequence( vseq );
    if (MB_SUCCESS != rval) {
      SequenceData* vdata = vseq->data();
      delete vseq;
      if (!seq_data)
        delete vdata;
    
      return rval;
    }
  }
  else {  
    vseq = reinterpret_cast<VertexSequence*>(*seq);
    if (append) {
      vseq->push_back( 1 );
      handle = vseq->end_handle();
      typeData[MBVERTEX].notify_appended( seq );
    }
    else {
      vseq->push_front( 1 );
      handle = vseq->start_handle();
      typeData[MBVERTEX].notify_prepended( seq );
    }
  }
  
  return vseq->set_coordinates( handle, coords );
}

  
ErrorCode SequenceManager::create_element( EntityType type,
                                             const EntityHandle* conn,
                                             unsigned conn_len,
                                             EntityHandle& handle )
{
  if (type <= MBVERTEX || type >= MBENTITYSET)
    return MB_TYPE_OUT_OF_RANGE;
  
  const EntityHandle start = CREATE_HANDLE( type, MB_START_ID );
  const EntityHandle   end = CREATE_HANDLE( type,   MB_END_ID );
  bool append;
  TypeSequenceManager::iterator seq = typeData[type].find_free_handle( start, end, append, conn_len );
  UnstructuredElemSeq* eseq;
  
  if (seq == typeData[type].end()) {
    SequenceData* seq_data = 0;
    unsigned size = DEFAULT_ELEMENT_SEQUENCE_SIZE;
    if (type == MBPOLYGON || type == MBPOLYHEDRON) {
      size = default_poly_sequence_size( conn_len );
    }
    EntityID seq_data_size = 0;
    handle = typeData[type].find_free_sequence( size, start, end, seq_data, seq_data_size, conn_len );
    if (!handle) 
      return MB_FAILURE;
    
    if (MBPOLYGON == type || MBPOLYHEDRON == type) {
      if (seq_data) 
        eseq = new PolyElementSeq( handle, 1, conn_len, seq_data );
      else
        eseq = new PolyElementSeq( handle, 1, conn_len, size );
    }
    else {
      if (seq_data) 
        eseq = new UnstructuredElemSeq( handle, 1, conn_len, seq_data );
      else
        eseq = new UnstructuredElemSeq( handle, 1, conn_len, size );
    }
    
    ErrorCode rval = typeData[type].insert_sequence( eseq );
    if (MB_SUCCESS != rval) {
      SequenceData* vdata = eseq->data();
      delete eseq;
      if (!seq_data)
        delete vdata;
    
      return rval;
    }
  }
  else {  
    eseq = reinterpret_cast<UnstructuredElemSeq*>(*seq);
    if (append) {
      eseq->push_back( 1 );
      handle = eseq->end_handle();
      typeData[type].notify_appended( seq );
    }
    else {
      eseq->push_front( 1 );
      handle = eseq->start_handle();
      typeData[type].notify_prepended( seq );
    }
  }
  
  return eseq->set_connectivity( handle, conn, conn_len );
}


  
ErrorCode SequenceManager::create_mesh_set( unsigned flags,
                                              EntityHandle& handle )
{
  const EntityHandle start = CREATE_HANDLE( MBENTITYSET, MB_START_ID );
  const EntityHandle   end = CREATE_HANDLE( MBENTITYSET,   MB_END_ID );
  bool append;
  TypeSequenceManager::iterator seq = typeData[MBENTITYSET].find_free_handle( start, end, append );
  MeshSetSequence* msseq;
  
  if (seq == typeData[MBENTITYSET].end()) {
    SequenceData* seq_data = 0;
    EntityID seq_data_size = 0;
    handle = typeData[MBENTITYSET].find_free_sequence( DEFAULT_MESHSET_SEQUENCE_SIZE, start, end, seq_data, seq_data_size );
    if (!handle) 
      return MB_FAILURE;
    
    if (seq_data) 
      msseq = new MeshSetSequence( handle, 1, flags, seq_data );
    else
      msseq = new MeshSetSequence( handle, 1, flags, DEFAULT_MESHSET_SEQUENCE_SIZE );
      
    ErrorCode rval = typeData[MBENTITYSET].insert_sequence( msseq );
    if (MB_SUCCESS != rval) {
      SequenceData* vdata = msseq->data();
      delete msseq;
      if (!seq_data)
        delete vdata;
    
      return rval;
    }
  }
  else {  
    msseq = reinterpret_cast<MeshSetSequence*>(*seq);
    if (append) {
      msseq->push_back( 1, &flags );
      handle = msseq->end_handle();
      typeData[MBENTITYSET].notify_appended( seq );
    }
    else {
      msseq->push_front( 1, &flags );
      handle = msseq->start_handle();
      typeData[MBENTITYSET].notify_prepended( seq );
    }
  }
  
  return MB_SUCCESS;
}

ErrorCode SequenceManager::allocate_mesh_set( EntityHandle handle,
                                                unsigned flags )
{
  SequenceData* data = 0;
  TypeSequenceManager::iterator seqptr; 
  EntityHandle block_start = 1, block_end = 0;
  ErrorCode rval = typeData[MBENTITYSET].is_free_handle( handle, seqptr, data, block_start, block_end );
  if (MB_SUCCESS != rval)
    return rval;
  
  MeshSetSequence* seq;
  if (seqptr != typeData[MBENTITYSET].end()) {
    seq = static_cast<MeshSetSequence*>(*seqptr);
    if (seq->start_handle() - 1 == handle) {
      rval = seq->push_front( 1, &flags );
      if (MB_SUCCESS == rval) {
        rval = typeData[MBENTITYSET].notify_prepended( seqptr );
        if (MB_SUCCESS != rval)
          seq->pop_front( 1 );
      }
      return rval;
    }
    else if (seq->end_handle() + 1 == handle) {
      rval = seq->push_back( 1, &flags );
      if (MB_SUCCESS == rval) {
        rval = typeData[MBENTITYSET].notify_appended( seqptr );
        if (MB_SUCCESS != rval)
          seq->pop_back( 1 );
      }
      return rval;
    }
    else
      return MB_FAILURE; // should be unreachable
  }
  else {
    if (data) {
      seq = new MeshSetSequence( handle, 1, flags, data );
    }
    else {
      assert( handle >= block_start && handle <= block_end );
      trim_sequence_block( handle, block_end, DEFAULT_MESHSET_SEQUENCE_SIZE );
      seq = new MeshSetSequence( handle, 1, flags, block_end - handle + 1 );
    }
    
    rval = typeData[MBENTITYSET].insert_sequence( seq );
    if (MB_SUCCESS != rval) {
      SequenceData* vdata = seq->data();
      delete seq;
      if (!data)
        delete vdata;
      return rval;
    }
  
    return MB_SUCCESS;
  }
}

void
SequenceManager::trim_sequence_block( EntityHandle start_handle,
                                      EntityHandle& end_handle,
                                      unsigned max_size )
{
  assert( end_handle >= start_handle );
  assert( (int)max_size > 0 ); // cast to int also prohibits some rediculously large values
  
    // if input range is larger than preferred size, trim it
  if (end_handle - start_handle >= max_size)
    end_handle = start_handle + max_size - 1;
}

EntityHandle 
SequenceManager::sequence_start_handle( EntityType type,
                                        EntityID count,
                                        int size,
                                        EntityID start,
                                        SequenceData*& data,
                                        EntityID &data_size)
{
  TypeSequenceManager &tsm = typeData[type];
  data = 0;
  EntityHandle handle = CREATE_HANDLE( type, start );
  if (start < MB_START_ID ||
      !tsm.is_free_sequence( handle, count, data, size )) {
    EntityHandle pstart = CREATE_HANDLE( type, MB_START_ID );
    EntityHandle pend   = CREATE_HANDLE( type,   MB_END_ID );
    handle = tsm.find_free_sequence( count, pstart, pend, data, data_size, size);
  }
  return handle;
}


EntityID SequenceManager::new_sequence_size( EntityHandle start,
                                               EntityID requested_size,
                                               EntityID default_size ) const
{
  if (requested_size >= default_size)
    return requested_size;
  
  EntityHandle last = typeData[TYPE_FROM_HANDLE(start)].last_free_handle( start );
// tjt - when start is 41427, last comes back 41685, when there's really an entity
    // at 41673, and 41467+246-1=41672
  if (!last) {
    assert( false );
    return 0;
  }
  
  EntityID available_size = last - start + 1;
  if (default_size < available_size)
    return default_size;
  else
    return available_size;
}

ErrorCode 
SequenceManager::create_entity_sequence( EntityType type,
                                         EntityID count,
                                         int size,
                                         EntityID start,
                                         EntityHandle& handle,
                                         EntitySequence*& sequence )
{
  SequenceData* data = 0;
  EntityID data_size = 0;
  handle = sequence_start_handle( type, count, size, start, data, data_size );
    
  if (!handle)
    return MB_MEMORY_ALLOCATION_FAILED;
  
  switch (type) {
  case MBENTITYSET:
  case MBMAXTYPE:
    return MB_TYPE_OUT_OF_RANGE;
  
  case MBVERTEX:
    if (size != 0)
      return MB_INDEX_OUT_OF_RANGE;

    if (data)
      sequence = new VertexSequence( handle, count, data );
    else {
      if (!data_size)
        data_size = new_sequence_size(handle, count, 
                                      DEFAULT_VERTEX_SEQUENCE_SIZE);
      sequence = new VertexSequence( handle, count, data_size );
    }
    break;
  
  case MBPOLYGON:
  case MBPOLYHEDRON:
    if (size == 0)
      return MB_INDEX_OUT_OF_RANGE;

    if (data)
      sequence = new PolyElementSeq( handle, count, size, data );
    else {
      if (!data_size)
        data_size = new_sequence_size(handle, count, 
                                      default_poly_sequence_size(size));
      sequence = new PolyElementSeq( handle, count, size, data_size );
    }
    break;
  
  default:
    if (size == 0)
      return MB_INDEX_OUT_OF_RANGE;

    if (data)
      sequence = new UnstructuredElemSeq( handle, count, size, data );
    else {
      if (!data_size)
        data_size = new_sequence_size(handle, count, 
                                      DEFAULT_ELEMENT_SEQUENCE_SIZE);
      sequence = new UnstructuredElemSeq( handle, count, size, data_size );
    }
      // tjt calling new_sequence_size 'cuz don't have a sequence data;
      // start 41467, count 246
    break;
  }
  
  ErrorCode result = typeData[type].insert_sequence( sequence );
  if (MB_SUCCESS != result) {
      // change to NULL if had an existing data or if no existing data,
      // change to the new data created
    data = data ? 0 : sequence->data();
    delete sequence;
    delete data;
    return result;
  }
  
  return MB_SUCCESS;
}


ErrorCode 
SequenceManager::create_meshset_sequence( EntityID count,
                                          EntityID start,
                                          const unsigned* flags,
                                          EntityHandle& handle,
                                          EntitySequence*& sequence )
{
  SequenceData* data = 0;
  EntityID data_size = 0;
  handle = sequence_start_handle( MBENTITYSET, count, 0, start, data, data_size );

  if (!handle)
    return MB_MEMORY_ALLOCATION_FAILED;
  
  if (data)
    sequence = new MeshSetSequence( handle, count, flags, data );
  else
    sequence = new MeshSetSequence( handle, count, flags, count );
  
  
  
  ErrorCode result = typeData[MBENTITYSET].insert_sequence( sequence );
  if (MB_SUCCESS != result) {
      // change to NULL if had an existing data or if no existing data,
      // change to the new data created
    data = data ? 0 : sequence->data();
    delete sequence;
    delete data;
    return result;
  }
  
  return MB_SUCCESS;
}


ErrorCode 
SequenceManager::create_meshset_sequence( EntityID count,
                                          EntityID start,
                                          unsigned flags,
                                          EntityHandle& handle,
                                          EntitySequence*& sequence )
{
  SequenceData* data = 0;
  EntityID data_size = 0;
  handle = sequence_start_handle( MBENTITYSET, count, 0, start, data, data_size );
  if (!handle)
    return MB_MEMORY_ALLOCATION_FAILED;
  
  if (data)
    sequence = new MeshSetSequence( handle, count, flags, data );
  else
    sequence = new MeshSetSequence( handle, count, flags, count );
  
  
  
  ErrorCode result = typeData[MBENTITYSET].insert_sequence( sequence );
  if (MB_SUCCESS != result) {
      // change to NULL if had an existing data or if no existing data,
      // change to the new data created
    data = data ? 0 : sequence->data();
    delete sequence;
    delete data;
    return result;
  }
  
  return MB_SUCCESS;
}

ErrorCode
SequenceManager::create_scd_sequence( int imin, int jmin, int kmin,
                                      int imax, int jmax, int kmax,
                                      EntityType type,
                                      EntityID start_id_hint,
                                      EntityHandle& handle,
                                      EntitySequence*& sequence )
{
  int this_dim = CN::Dimension(type);

    // use > instead of != in the following assert to also catch cases where imin > imax, etc.
  assert((this_dim < 3 || kmax > kmin) &&
         (this_dim < 2 || jmax > jmin) &&
         (this_dim < 1 || imax > imin));

    // compute # entities; not as easy as it would appear...
  EntityID num_ent;
  if (MBVERTEX == type)
    num_ent = (EntityID)(imax-imin+1)*(EntityID)(jmax-jmin+1)*(EntityID)(kmax-kmin+1);
  else {
    num_ent = (imax-imin) *
      (this_dim >= 2 ? (jmax-jmin) : 1) *
      (this_dim >= 3 ? (kmax-kmin) : 1);
  }
  
    // get a start handle
  SequenceData* data = 0;
  EntityID data_size = 0;
  handle = sequence_start_handle( type, num_ent, -1, start_id_hint, data, data_size );

  if (!handle)
    return MB_MEMORY_ALLOCATION_FAILED;
  assert(!data);
  
  switch (type) {
  case MBVERTEX:
    data = new ScdVertexData( handle, imin, jmin, kmin, imax, jmax, kmax );
    sequence = new VertexSequence( handle, data->size(), data );
    break;
  case MBEDGE:
  case MBQUAD:
  case MBHEX:
    sequence = new StructuredElementSeq( handle, imin, jmin, kmin, imax, jmax, kmax );
    break;
  default:
    return MB_TYPE_OUT_OF_RANGE;
  }
  
  ErrorCode result = typeData[type].insert_sequence( sequence );
  if (MB_SUCCESS != result) {
    data = sequence->data();
    delete sequence;
    delete data;
    return result;
  }
  
  return MB_SUCCESS;
}

ErrorCode
SequenceManager::create_scd_sequence( const HomCoord& coord_min,
                                      const HomCoord& coord_max,
                                      EntityType type,
                                      EntityID start_id_hint,
                                      EntityHandle& first_handle_out,
                                      EntitySequence*& sequence_out )
{
  return create_scd_sequence( coord_min.i(), coord_min.j(), coord_min.k(),
                              coord_max.i(), coord_max.j(), coord_max.k(),
                              type, start_id_hint,
                              first_handle_out, sequence_out );
}

ErrorCode
SequenceManager::create_sweep_sequence( int imin, int jmin, int kmin,
					int imax, int jmax, int kmax,
					int* Cq,
					EntityType type,
					EntityID start_id_hint,
					EntityHandle& handle,
					EntitySequence*& sequence )
{
  int this_dim = CN::Dimension(type);

  assert((this_dim < 3 || kmax > kmin) &&
         (this_dim < 2 || jmax > jmin) &&
         (this_dim < 1 || imax > imin));

  EntityID num_ent;
  if (MBVERTEX == type)
    num_ent = (EntityID)(imax-imin+1)*(EntityID)(jmax-jmin+1)*(EntityID)(kmax-kmin+1);
  else {
    num_ent = (imax-imin) *
      (this_dim >= 2 ? (jmax-jmin) : 1) *
      (this_dim >= 3 ? (kmax-kmin) : 1);
  }
  
    // get a start handle
  SequenceData* data = 0;
  EntityID data_size = 0;
  handle = sequence_start_handle( type, num_ent, -1, start_id_hint, data, data_size );

  if (!handle)
    return MB_MEMORY_ALLOCATION_FAILED;
  assert(!data);
  
  switch (type) {
  case MBVERTEX:
    data = new ScdVertexData( handle, imin, jmin, kmin, imax, jmax, kmax );
    sequence = new VertexSequence( handle, data->size(), data );
    break;
  case MBEDGE:
  case MBQUAD:
  case MBHEX:
    sequence = new SweptElementSeq( handle, imin, jmin, kmin, imax, jmax, kmax, Cq );
    break;
  default:
    return MB_TYPE_OUT_OF_RANGE;
  }
  
  ErrorCode result = typeData[type].insert_sequence( sequence );
  if (MB_SUCCESS != result) {
    data = sequence->data();
    delete sequence;
    delete data;
    return result;
  }
  
  return MB_SUCCESS;
}

ErrorCode
SequenceManager::create_sweep_sequence( const HomCoord& coord_min,
					const HomCoord& coord_max,
					int* Cq,
					EntityType type,
					EntityID start_id_hint,
					EntityHandle& first_handle_out,
					EntitySequence*& sequence_out )
{
  return create_sweep_sequence( coord_min.i(), coord_min.j(), coord_min.k(),
				coord_max.i(), coord_max.j(), coord_max.k(),
				Cq,
				type, start_id_hint,
				first_handle_out, sequence_out );
}

ErrorCode 
SequenceManager::add_vsequence(EntitySequence *vert_seq,
                               EntitySequence *elem_seq,
                               const HomCoord &p1, const HomCoord &q1,
                               const HomCoord &p2, const HomCoord &q2,
                               const HomCoord &p3, const HomCoord &q3,
                               bool bb_input,
                               const HomCoord *bb_min,
                               const HomCoord *bb_max) 
{
    // check first that they're structured vtx/elem sequences
  ScdVertexData *scd_vd = dynamic_cast<ScdVertexData*>(vert_seq->data());
  if (!scd_vd) return MB_FAILURE;
  
  ScdElementData *scd_ed = dynamic_cast<ScdElementData*>(elem_seq->data());
  if (!scd_ed) return MB_FAILURE;

  if (bb_min && bb_max)
    return scd_ed->add_vsequence(scd_vd, p1, q1, p2, q2, p3, q3, 
                                 bb_input, *bb_min, *bb_max);
  else
    return scd_ed->add_vsequence(scd_vd, p1, q1, p2, q2, p3, q3, 
                                 bb_input, HomCoord::unitv[0], HomCoord::unitv[0]);
}
 
ErrorCode
SequenceManager::replace_subsequence( EntitySequence* new_seq, TagServer* ts )
{
  const EntityType type = TYPE_FROM_HANDLE(new_seq->start_handle());
  return typeData[type].replace_subsequence( new_seq, ts );
}

void SequenceManager::get_memory_use( unsigned long& total_entity_storage,
                                      unsigned long& total_storage ) const

{
  total_entity_storage = 0;
  total_storage = 0;
  unsigned long temp_entity, temp_total;
  for (EntityType i = MBVERTEX; i < MBMAXTYPE; ++i) {
    temp_entity = temp_total = 0;
    get_memory_use( i, temp_entity, temp_total );
    total_entity_storage += temp_entity;
    total_storage        += temp_total;
  }
}

void SequenceManager::get_memory_use( EntityType type,
                                      unsigned long& total_entity_storage,
                                      unsigned long& total_storage ) const
{
  typeData[type].get_memory_use( total_entity_storage, total_storage );
}

void SequenceManager::get_memory_use( const Range& entities,
                                      unsigned long& total_entity_storage,
                                      unsigned long& total_amortized_storage ) const
{
  total_entity_storage = 0;
  total_amortized_storage = 0;
  unsigned long temp_entity, temp_total;
  Range::const_pair_iterator i;
  for (i = entities.const_pair_begin(); i != entities.const_pair_end(); ++i) {
    const EntityType t1 = TYPE_FROM_HANDLE(i->first);
    const EntityType t2 = TYPE_FROM_HANDLE(i->second);
    if (t1 == t2) {
      temp_entity = temp_total = 0;
      typeData[t1].get_memory_use( i->first, i->second, temp_entity, temp_total );
      total_entity_storage += temp_entity;
      total_amortized_storage += temp_total;
    }
    else {
      int junk;
 
      temp_entity = temp_total = 0;
      typeData[t1].get_memory_use( i->first, CREATE_HANDLE(t1,MB_END_ID,junk), temp_entity, temp_total );
      total_entity_storage += temp_entity;
      total_amortized_storage += temp_total;

      temp_entity = temp_total = 0;
      typeData[t2].get_memory_use( CREATE_HANDLE(t2,MB_START_ID,junk), i->second, temp_entity, temp_total );
      total_entity_storage += temp_entity;
      total_amortized_storage += temp_total;
    }
  }
}

void SequenceManager::reset_tag_data() {
  for (EntityType t = MBVERTEX; t <= MBENTITYSET; ++t) {
    TypeSequenceManager& seqs = entity_map(t);
    for (TypeSequenceManager::iterator i = seqs.begin(); i != seqs.end(); ++i)
      (*i)->data()->release_tag_data( &tagSizes[0], tagSizes.size() );
  }
}

ErrorCode SequenceManager::reserve_tag_id( int size, TagId tag_id )
{
  if (size < 1 && size != MB_VARIABLE_LENGTH)
    return MB_INVALID_SIZE;
  if (tag_id >= tagSizes.size())
    tagSizes.resize( tag_id+1, 0 );
  if (tagSizes[tag_id])
    return MB_ALREADY_ALLOCATED;
  tagSizes[tag_id] = size;
  return MB_SUCCESS;
}

ErrorCode SequenceManager::release_tag( TagId tag_id )
{
  if (tag_id >= tagSizes.size() || !tagSizes[tag_id])
    return MB_TAG_NOT_FOUND;
  tagSizes[tag_id] = 0;
  
  for (EntityType t = MBVERTEX; t <= MBENTITYSET; ++t) {
    TypeSequenceManager& seqs = entity_map(t);
    for (TypeSequenceManager::iterator i = seqs.begin(); i != seqs.end(); ++i)
      (*i)->data()->release_tag_data(tag_id, tagSizes[tag_id]);
  }
  return MB_SUCCESS;
}

ErrorCode SequenceManager::remove_tag_data( TagId tag_id, 
                                              EntityHandle handle,
                                              const void* default_tag_value,
                                              int default_value_size )
{
  if (tag_id >= tagSizes.size() || !tagSizes[tag_id])
    return MB_TAG_NOT_FOUND;

  EntitySequence* seq = 0;
  ErrorCode rval = find( handle, seq );
  if (MB_SUCCESS != rval)
    return rval;
  
  if (tagSizes[tag_id] == MB_VARIABLE_LENGTH) {
    VarLenTag* tag_data = vtag_array( seq, handle, tag_id );
    if (!tag_data)
      return MB_TAG_NOT_FOUND;
    VarLenTag* vdata = reinterpret_cast<VarLenTag*>(tag_data);
    if (default_tag_value)
      vdata->set( default_tag_value, default_value_size );
    else
      vdata->clear();
  }
  else {
    void* tag_data = tag_array( seq, handle, tag_id, tagSizes[tag_id] );
    if (!tag_data)
      return MB_TAG_NOT_FOUND;
    if (default_tag_value)  
      memcpy( tag_data, default_tag_value, tagSizes[tag_id] );
    else
      memset( tag_data, 0, tagSizes[tag_id] );
  }
  return MB_SUCCESS;
}

ErrorCode SequenceManager::set_tag_data( TagId tag_id,
                                           const EntityHandle* handles,
                                           int num_handles,
                                           const void* values,
                                           const void* default_value )
{
  if (tag_id >= tagSizes.size() || tagSizes[tag_id] < 1) {
    if (tag_id < tagSizes.size() && tagSizes[tag_id] == MB_VARIABLE_LENGTH)
      return MB_VARIABLE_DATA_LENGTH;
    else
      return MB_TAG_NOT_FOUND;
  }
  
  ErrorCode result = MB_SUCCESS;
  const unsigned char* ptr = reinterpret_cast<const unsigned char*>(values);
  const EntityHandle* const end = handles + num_handles;
  for (const EntityHandle* i = handles; i != end; ++i, ptr += tagSizes[tag_id] ) {
    EntitySequence* seq = 0;
    ErrorCode rval = find( *i, seq );
    if (MB_SUCCESS != rval) {
      result = rval;
      continue;
    }
  
    unsigned char* tag_data = make_tag( seq, *i, tag_id, tagSizes[tag_id], default_value );
    memcpy( tag_data, ptr, tagSizes[tag_id] );
  }

  return result;
}

ErrorCode SequenceManager::set_tag_data( TagId tag_id,
                                           const EntityHandle* handles,
                                           int num_handles,
                                           void const* const* values,
                                           const int* lengths,
                                           const void* default_value,
                                           bool one_value )
{
  ErrorCode result = MB_SUCCESS;
  const EntityHandle* const end = handles + num_handles;
  const bool step = !one_value;
  
  if (tagSizes[tag_id] == MB_VARIABLE_LENGTH) {
    if (!lengths)
      return MB_VARIABLE_DATA_LENGTH;
    
    for (const EntityHandle* i = handles; i != end; ++i, values += step, lengths += step ) {
        // find sequence for entity
      EntitySequence* seq = 0;
      ErrorCode rval = find( *i, seq );
      if (MB_SUCCESS != rval) {
        result = rval;
        continue;
      }
        // set value
      VarLenTag* data = make_vtag( seq, *i, tag_id );
      data->set( *values, *lengths );
    }
  }
  else {
    for (const EntityHandle* i = handles; i != end; ++i, values += step ) {
        // find sequence for entity
      EntitySequence* seq = 0;
      ErrorCode rval = find( *i, seq );
      if (MB_SUCCESS != rval) {
        result = rval;
        continue;
      }
        // set value
      unsigned char* data = make_tag( seq, *i, tag_id, tagSizes[tag_id], default_value );
      memcpy( data, *values, tagSizes[tag_id] );
    }
  }
  
  return result;
}

ErrorCode SequenceManager::set_tag_data( TagId tag_id,
                                           const Range& handles,
                                           const void* values,
                                           const void* default_value )
{
  ErrorCode rval, result = MB_SUCCESS;
    // NOTE: Comparison of size to 1 should also catch 
    //       case where tag is variable-length.  
  if (tag_id >= tagSizes.size() || tagSizes[tag_id] < 1) {
    if (tag_id < tagSizes.size() && tagSizes[tag_id] == MB_VARIABLE_LENGTH)
      return MB_VARIABLE_DATA_LENGTH;
    else
      return MB_TAG_NOT_FOUND;
  }
  
  const char* data = reinterpret_cast<const char*>(values);

  for (Range::const_pair_iterator p = handles.const_pair_begin(); 
       p != handles.const_pair_end(); ++p) {
       
    EntityHandle start = p->first;
    while (start <= p->second) {
      
      EntitySequence* seq = 0;
      rval = find( start, seq );
      if (MB_SUCCESS != rval) {
        result = rval;
        ++start;
        data += tagSizes[tag_id];
        continue;
      }
      
      const EntityHandle finish = std::min( p->second, seq->end_handle() );
      const EntityID count = finish - start + 1;
      
      void* tag_array = seq->data()->get_tag_data( tag_id );
      if (!tag_array)
        tag_array = seq->data()->create_tag_data( tag_id, tagSizes[tag_id], default_value );

      char* tag_data = reinterpret_cast<char*>(tag_array) + 
                       tagSizes[tag_id] * (start - seq->data()->start_handle());
      memcpy( tag_data, data, tagSizes[tag_id] * count );
      data += tagSizes[tag_id] * count;
    
      start = finish + 1;
    }
  }
  
  return result;
}
      

ErrorCode SequenceManager::set_tag_data( TagId tag_id,
                                           const Range& handles,
                                           void const* const* values,
                                           const int* lengths,
                                           const void* default_value,
                                           bool one_value )
{
  ErrorCode rval, result = MB_SUCCESS;
  if (tag_id >= tagSizes.size() || !tagSizes[tag_id])
    return MB_TAG_NOT_FOUND;

  int alloc_size = tagSizes[tag_id];
  if (alloc_size == MB_VARIABLE_LENGTH) {
    alloc_size = sizeof(VarLenTag);
  
    if (!lengths)
      return MB_VARIABLE_DATA_LENGTH;
    
    // Ignore default value for var-len tags.  Just zero the array
    // data, which results in VarLenTag structs with zero-length
    // per-entity values.  The query code will return the default
    // value for such cases.  
    default_value = NULL;
  }

  const bool step = !one_value;

  for (Range::const_pair_iterator p = handles.const_pair_begin(); 
       p != handles.const_pair_end(); ++p) {
       
    EntityHandle start = p->first;
    while (start <= p->second) {
      
      EntitySequence* seq = 0;
      rval = find( start, seq );
      if (MB_SUCCESS != rval) {
        result = rval;
        ++start;
        values += step;
        if (lengths)
          lengths += step;
        continue;
      }
      
      const EntityHandle finish = std::min( p->second, seq->end_handle() ) + 1;
      const EntityID count = finish - start;
      
      void* tag_array = seq->data()->get_tag_data( tag_id );
      if (!tag_array)
        tag_array = seq->data()->create_tag_data( tag_id, alloc_size, default_value );

      if (tagSizes[tag_id] == MB_VARIABLE_LENGTH) {
        VarLenTag* tag_data = reinterpret_cast<VarLenTag*>(tag_array) +
                              start - seq->data()->start_handle();
        VarLenTag* end_data = tag_data + count;
        while (tag_data != end_data) {
          tag_data->set( *values, *lengths );
          ++tag_data;
          values += step;
          lengths += step;
        }
      }
      else {
        char* tag_data = reinterpret_cast<char*>(tag_array) + 
                         alloc_size * (start - seq->data()->start_handle());
        char* end_data = tag_data + alloc_size * count;
        while (tag_data != end_data) {
          memcpy( tag_data, *values, alloc_size );
          tag_data += alloc_size;
          values += step;
        }
      }
    
      start = finish;
    }
  }
  
  return result;
}
      

ErrorCode SequenceManager::tag_iterate( TagId tag_id,
                                        Range::iterator& iter,
                                        const Range::iterator& end,
                                        void*& data_ptr_out,
                                        const void* default_value )
{
  if (tag_id >= tagSizes.size() || tagSizes[tag_id] < 1) {
    if (tag_id < tagSizes.size() && tagSizes[tag_id] == MB_VARIABLE_LENGTH)
      return MB_VARIABLE_DATA_LENGTH;
    else
      return MB_TAG_NOT_FOUND;
  }

    // If asked for nothing, successfully return nothing.
  if (iter == end)
    return MB_SUCCESS;
  
  EntitySequence* seq = 0;
  ErrorCode rval = find( *iter, seq );
  if (MB_SUCCESS != rval)
    return rval;
  
  size_t count1 = *(iter.end_of_block()) - *iter + 1;
  size_t count2 = seq->end_handle() - *iter + 1;
  size_t count = std::min( count1, count2 );
  size_t offset = *iter - seq->data()->start_handle();
  
  void* tag_array = seq->data()->get_tag_data( tag_id );
  if (!tag_array)
    tag_array = seq->data()->create_tag_data( tag_id, tagSizes[tag_id], default_value );
  data_ptr_out = (unsigned char*)tag_array + offset*tagSizes[tag_id];
  iter += count;
  return MB_SUCCESS;
}

ErrorCode SequenceManager::get_tag_data( TagId tag_id,
                                           const EntityHandle* handles,
                                           int num_handles,
                                           void* values,
                                           const void* default_value ) const
{
    // NOTE: Comparison of size to 1 should also catch 
    //       case where tag is variable-length.  
  if (tag_id >= tagSizes.size() || tagSizes[tag_id] < 1) {
    if (tag_id < tagSizes.size() && tagSizes[tag_id] == MB_VARIABLE_LENGTH)
      return MB_VARIABLE_DATA_LENGTH;
    else
      return MB_TAG_NOT_FOUND;
  }

  const int len = tagSizes[tag_id];
  unsigned char* ptr = reinterpret_cast<unsigned char*>(values);
  const EntityHandle *const end = handles + num_handles;
  for (const EntityHandle* i = handles; i != end; ++i, ptr += len) {
    
    const EntitySequence* seq = 0;
    ErrorCode rval = find( *i, seq );
    // keep MOAB 3.0 behavior : return default value for invalid handles
    if (MB_SUCCESS != rval)
      return rval;
  
    const unsigned char* tag_data = tag_array( seq, *i, tag_id, len );
    if (!tag_data) {
      if (default_value) 
        memcpy( ptr, default_value, tagSizes[tag_id] );
      else
        return MB_TAG_NOT_FOUND;
    } 
    else {
      memcpy( ptr, tag_data, len );
    }
  }
  return MB_SUCCESS;
}

ErrorCode SequenceManager::get_tag_data( TagId tag_id,
                                           const EntityHandle* handles,
                                           int num_handles,
                                           const void** values,
                                           int* lengths,
                                           const void* default_value,
                                           int default_value_length ) const
{
  if (tag_id >= tagSizes.size() || !tagSizes[tag_id])
    return MB_TAG_NOT_FOUND;
  
  ErrorCode result = MB_SUCCESS;
  const EntityHandle *const end = handles + num_handles;
  const int len = tagSizes[tag_id];
  
  if (len == MB_VARIABLE_LENGTH) {
    for (const EntityHandle* i = handles; i != end; ++i) {
      
      const EntitySequence* seq = 0;
      ErrorCode rval = find( *i, seq );
      if (MB_SUCCESS != rval) {
        result = rval;
        *values = 0;
        *lengths = 0;
      }
      else {
        const VarLenTag* tag_data = vtag_array( seq, *i, tag_id );
        if (tag_data && tag_data->size()) {
          *values = tag_data->data();
          *lengths = tag_data->size();
        }
        else if (default_value) {
          *values = default_value;
          *lengths = default_value_length;
        }
        else {
          result = MB_TAG_NOT_FOUND;
          *values = 0;
          *lengths = 0;
        }
      }
      
      ++values;
      ++lengths;
    }
  }
  else {
    if (lengths) 
      SysUtil::setmem( lengths, &len, sizeof(int), num_handles );
  
    for (const EntityHandle* i = handles; i != end; ++i) {
      const EntitySequence* seq = 0;
      ErrorCode rval = find( *i, seq );
      if (MB_SUCCESS != rval) {
        result = rval;
        *values = 0;
      }
      else {
        *values = tag_array( seq, *i, tag_id, len ); 
        if (!*values) {
          if (default_value)
            *values = default_value;
          else
            result = MB_TAG_NOT_FOUND;
        }
      }
      ++values;
    }
  }
  
  return result;
}

ErrorCode SequenceManager::get_tag_data( TagId tag_id,
                                           const Range& handles,
                                           void* values,
                                           const void* default_value ) const
{
  ErrorCode rval;
    // NOTE: Comparison of size to 1 should also catch 
    //       case where tag is variable-length.  
  if (tag_id >= tagSizes.size() || tagSizes[tag_id] < 1) {
    if (tag_id < tagSizes.size() && tagSizes[tag_id] == MB_VARIABLE_LENGTH)
      return MB_VARIABLE_DATA_LENGTH;
    else
      return MB_TAG_NOT_FOUND;
  }
    
  char* data = reinterpret_cast<char*>(values);

  for (Range::const_pair_iterator p = handles.const_pair_begin(); 
       p != handles.const_pair_end(); ++p) {
       
    EntityHandle start = p->first;
    while (start <= p->second) {
      
      const EntitySequence* seq = 0;
      rval = find( start, seq );
      if (MB_SUCCESS != rval)
        return rval;
     
      const EntityHandle finish = std::min( p->second, seq->end_handle() );
      const EntityID count = finish - start + 1;
      
      const void* tag_array = seq->data()->get_tag_data( tag_id );
      if (tag_array) {
        const char* tag_data = reinterpret_cast<const char*>(tag_array) + 
                         tagSizes[tag_id] * (start - seq->data()->start_handle());
        memcpy( data, tag_data, tagSizes[tag_id] * count );
      }
      else if (default_value) {
        SysUtil::setmem( data, default_value, tagSizes[tag_id], count );
      }
      else {
        return MB_TAG_NOT_FOUND;
      }
      
      data += tagSizes[tag_id] * count;
      start = finish + 1;
    }
  }
  
  return MB_SUCCESS;
}

ErrorCode SequenceManager::get_tag_data( TagId tag_id,
                                           const Range& handles,
                                           const void** values,
                                           int* lengths,
                                           const void* default_value,
                                           int default_value_length ) const
{
  ErrorCode rval, result = MB_SUCCESS;
  if (!default_value)
    default_value_length = 0;

  if (tag_id >= tagSizes.size() || !tagSizes[tag_id]) 
    return MB_TAG_NOT_FOUND;
  
  if (tagSizes[tag_id] == MB_VARIABLE_LENGTH) {
    if (!lengths || (default_value && !default_value_length))
      return MB_VARIABLE_DATA_LENGTH;
  }
  else if (lengths) {
    int len = tagSizes[tag_id];
    SysUtil::setmem( lengths, &len, sizeof(int), handles.size() );
  }
  
  for (Range::const_pair_iterator p = handles.const_pair_begin(); 
       p != handles.const_pair_end(); ++p) {

    EntityHandle start = p->first;
    while (start <= p->second) {

      const EntitySequence* seq = 0;
      rval = find( start, seq );
      if (MB_SUCCESS != rval) {
        *values = 0; 
        ++values;
        if (lengths) {
          *lengths = 0;
          ++lengths;
        }
        result = rval;
        ++start;
        continue;
      }

      const EntityHandle finish = std::min( p->second, seq->end_handle() ) + 1;
      const EntityID count = finish - start;
      if (tagSizes[tag_id] == MB_VARIABLE_LENGTH) {
        const VarLenTag* tag_data = vtag_array( seq, start, tag_id );
        if (!tag_data) {
          SysUtil::setmem( values, &default_value, sizeof(void*), count );
          SysUtil::setmem( lengths, &default_value_length, sizeof(int), count );
          values += count;
          lengths += count;
          if (!default_value)
            result = MB_TAG_NOT_FOUND;
          continue;
        }

        const VarLenTag* end_data = tag_data + count;
        while (tag_data != end_data) {
          if (tag_data->size()) {
            *values = tag_data->data();
            *lengths = tag_data->size();
          }
          else if (default_value) {
            *values = default_value;
            *lengths = default_value_length;
          }
          else {
            *values = 0;
            *lengths = 0;
            result = MB_TAG_NOT_FOUND;
          }
          ++values;
          ++lengths;
          ++tag_data;
        }
      }
      else {
        const unsigned char* tag_data = tag_array( seq, start, tag_id, tagSizes[tag_id] );
        if (!tag_data) {
          SysUtil::setmem( values, &default_value, sizeof(void*), count );
          values += count;
          if (!default_value)
            result = MB_TAG_NOT_FOUND;
        }
        else {
          const unsigned char* end_data = tag_data + count * tagSizes[tag_id];
          while (tag_data != end_data) {
            *values = tag_data; 
            ++values;
            tag_data += tagSizes[tag_id];
          }
        }
      }
      start = finish;
    }
  }
  
  return result;
}

ErrorCode SequenceManager::get_entity_tags(  EntityHandle entity,
                                 std::vector<Tag>& tags_out ) const
{
  const EntitySequence* seq = 0;
  ErrorCode rval = find( entity, seq );
  if (MB_SUCCESS != rval)
    return rval;
  
  for (TagId i = 0; i < tagSizes.size(); ++i) {
    if (tagSizes[i] == MB_VARIABLE_LENGTH) {
      const void* data_array = seq->data()->get_tag_data(i);
      if (data_array) {
        const VarLenTag* tag_ptr = reinterpret_cast<const VarLenTag*>(data_array);
        tag_ptr += (entity - seq->data()->start_handle());
        if (tag_ptr->size())
          tags_out.push_back( TAG_HANDLE_FROM_ID( i, MB_TAG_DENSE ) );
      }
    }
    else {
      if (seq->data()->get_tag_data(i))
        tags_out.push_back( TAG_HANDLE_FROM_ID( i, MB_TAG_DENSE ) );
    }
  }
  
  return MB_SUCCESS;
}

ErrorCode SequenceManager::get_tagged_entities( TagId tag_id, 
                                                  EntityType type,
                                                  Range& entities_out ) const
{
  if (tag_id >= tagSizes.size() || !tagSizes[tag_id])
    return MB_TAG_NOT_FOUND;

  if (type == MBMAXTYPE) {
    ErrorCode rval;
    for (EntityType t = MBVERTEX; t < MBMAXTYPE; ++t) {
      rval = get_tagged_entities( tag_id, t, entities_out );
      if (MB_SUCCESS != rval)
        return rval;
    }
    return MB_SUCCESS;
  }      

  Range::iterator insert = entities_out.begin();
  const TypeSequenceManager& map = entity_map( type );
  if (tagSizes[tag_id] == MB_VARIABLE_LENGTH) {
    const VarLenTag *data, *iter, *end;
    for (TypeSequenceManager::const_iterator i = map.begin(); i != map.end(); ++i) {
      data = reinterpret_cast<const VarLenTag*>((*i)->data()->get_tag_data(tag_id));
      if (!data)
        continue;
      end = data + (*i)->end_handle() - (*i)->data()->start_handle() + 1;
      iter = data + (*i)->start_handle() - (*i)->data()->start_handle();
      EntityHandle handle = (*i)->start_handle();
      for (; iter != end; ++iter, ++handle)
        if (iter->size()) 
          insert = entities_out.insert( insert, handle, handle );
    }
  }
  else {
    for (TypeSequenceManager::const_iterator i = map.begin(); i != map.end(); ++i) 
      if ((*i)->data()->get_tag_data(tag_id))
        insert = entities_out.insert( insert, (*i)->start_handle(), (*i)->end_handle() );
  }
  
  return MB_SUCCESS;
}

ErrorCode SequenceManager::count_tagged_entities( TagId tag_id, 
                                                    EntityType type,
                                                    int& count ) const
{
  if (tag_id >= tagSizes.size() || !tagSizes[tag_id])
    return MB_TAG_NOT_FOUND;

  count = 0;

  if (type == MBMAXTYPE) {
    ErrorCode rval;
    for (EntityType t = MBVERTEX; t < MBMAXTYPE; ++t) {
      int tmp = 0;
      rval = count_tagged_entities( tag_id, t, tmp );
      if (MB_SUCCESS != rval)
        return rval;
      count += tmp;
    }
    return MB_SUCCESS;
  }      

  const TypeSequenceManager& map = entity_map( type );
  
  if (tagSizes[tag_id] == MB_VARIABLE_LENGTH) {
    const VarLenTag *data, *iter, *end;
    for (TypeSequenceManager::const_iterator i = map.begin(); i != map.end(); ++i) {
      data = reinterpret_cast<const VarLenTag*>((*i)->data()->get_tag_data(tag_id));
      if (!data)
        continue;
      end = data + (*i)->end_handle() - (*i)->data()->start_handle();
      iter = data + (*i)->start_handle() - (*i)->data()->start_handle();
      for (; iter != end; ++iter)
        if (iter->size()) 
          ++count;
    }
  }
  else {
    for (TypeSequenceManager::const_iterator i = map.begin(); i != map.end(); ++i) 
      if ((*i)->data()->get_tag_data(tag_id))
        count += (*i)->size();
  }
  
  return MB_SUCCESS;
}

ErrorCode SequenceManager::get_entities_with_tag_value( TagId id,
                                                          const TagInfo& tag_info,
                                                          EntityType type,
                                                          Range& entities_out,
                                                          const void* value,
                                                          int size ) const
{
  if (id >= tagSizes.size() || !tagSizes[id])
    return MB_TAG_NOT_FOUND;

  if (type == MBMAXTYPE) {
    ErrorCode rval;
    for (EntityType t = MBVERTEX; t < MBMAXTYPE; ++t) {
      rval = get_entities_with_tag_value( id, tag_info, t, entities_out, value, size );
      if (MB_SUCCESS != rval)
        return rval;
    }
    return MB_SUCCESS;
  }      

  const TypeSequenceManager& map = entity_map( type );
  for (TypeSequenceManager::const_iterator i = map.begin(); i != map.end(); ++i) {
    if (const void* data = (*i)->data()->get_tag_data(id)) {
      ByteArrayIterator start( (*i)->data()->start_handle(), data, tag_info );
      ByteArrayIterator end( (*i)->end_handle() + 1, 0, 0 );
      start += (*i)->start_handle() - (*i)->data()->start_handle();
      find_tag_values_equal( tag_info, value, size, start, end, entities_out );
    }
  }
  
  return MB_SUCCESS;
}

ErrorCode SequenceManager::get_entities_with_tag_value( const Range& range,
                                                          TagId id,
                                                          const TagInfo& tag_info,
                                                          EntityType type,
                                                          Range& entities_out,
                                                          const void* value,
                                                          int size ) const
{
  ErrorCode rval;
  if (id >= tagSizes.size() || !tagSizes[id])
    return MB_TAG_NOT_FOUND;
    
  Range::const_pair_iterator p = type == MBMAXTYPE ? range.begin() : range.lower_bound(type);         
  for (; 
       p != range.const_pair_end() && 
       (MBMAXTYPE == type || TYPE_FROM_HANDLE(p->first) == type); 
       ++p) {
    
    EntityHandle start = p->first;
    while (start <= p->second) {
      
      const EntitySequence* seq = 0;
      rval = find( start, seq );
      if (MB_SUCCESS != rval) 
        return rval;
     
      const EntityHandle finish = std::min( p->second, seq->end_handle() );
      const void* tag_array = seq->data()->get_tag_data( id );
      if (tag_array) {
        ByteArrayIterator istart( seq->data()->start_handle(), tag_array, tag_info );
        ByteArrayIterator iend( seq->end_handle() + 1, 0, 0 );
        istart += seq->start_handle() - seq->data()->start_handle();
        find_tag_values_equal( tag_info, value, size, istart, iend, entities_out );
      }
      start = finish + 1;
    }
  }
  
  return MB_SUCCESS;
}

ErrorCode SequenceManager::get_tag_memory_use( TagId id, 
                                       unsigned long& total, 
                                       unsigned long& per_entity ) const
{
  if (id >= tagSizes.size() || !tagSizes[id])
    return MB_TAG_NOT_FOUND;
    
  per_entity = tagSizes[id];
  total = 0;
  for (EntityType t = MBVERTEX; t <= MBENTITYSET; ++t) {
    const TypeSequenceManager& map = entity_map(t);
    const SequenceData* prev_data = 0;
    for (TypeSequenceManager::const_iterator i = map.begin(); i != map.end(); ++i) {
      if ((*i)->data() != prev_data && (*i)->data()->get_tag_data(id)) {
        prev_data = (*i)->data();
        total += tagSizes[id] * (*i)->data()->size();
      }
    }
  }
      
  return MB_SUCCESS;
}
      

// These are meant to be called from the debugger (not declared in any header)
// so leave them out of release builds (-DNDEBUG).
#ifndef NDEBUG

std::ostream& operator<<( std::ostream& s, const TypeSequenceManager& seq_man )
{
  const SequenceData* prev_data = 0;
  for (TypeSequenceManager::const_iterator i = seq_man.begin(); i != seq_man.end(); ++i) {
    const EntitySequence* seq = *i;
    if (seq->data() != prev_data) {
      prev_data = seq->data();
      s << "SequenceData [" 
        << ID_FROM_HANDLE(seq->data()->start_handle())
        << ","
        << ID_FROM_HANDLE(seq->data()->end_handle())
        << "]"
        << std::endl;
    }
    s << "  Sequence [" 
      << ID_FROM_HANDLE(seq->start_handle())
      << ","
      << ID_FROM_HANDLE(seq->end_handle())
      << "]"
      << std::endl;
  }
  return s;
}

std::ostream& operator<<( std::ostream& s, const SequenceManager& seq_man )
{
  for (EntityType t = MBVERTEX; t < MBMAXTYPE; ++t) 
    if (!seq_man.entity_map(t).empty()) 
      s << std::endl 
        << "****************** " << CN::EntityTypeName( t ) << " ******************"
        << std::endl << seq_man.entity_map(t) << std::endl;
  return s;
}

void print_sequences( const SequenceManager& seqman )
{
  std::cout << seqman << std::endl;
}

void print_sequences( const TypeSequenceManager& seqman )
{
  std::cout << seqman << std::endl;
}

#endif

} // namespace moab
