#include "SequenceManager.hpp"
#include "VertexSequence.hpp"
#include "UnstructuredElemSeq.hpp"
#include "ScdVertexData.hpp"
#include "MeshSetSequence.hpp"
#include "StructuredElementSeq.hpp"
#include "HomXform.hpp"

#include <assert.h>

const int DEFAULT_VERTEX_SEQUENCE_SIZE = 4096;
const int DEFAULT_ELEMENT_SEQUENCE_SIZE = DEFAULT_VERTEX_SEQUENCE_SIZE;
const int DEFAULT_POLY_SEQUENCE_SIZE = 4 * DEFAULT_ELEMENT_SEQUENCE_SIZE;
const int DEFAULT_MESHSET_SEQUENCE_SIZE = DEFAULT_VERTEX_SEQUENCE_SIZE;

MBErrorCode SequenceManager::check_valid_entities( const MBRange& entities ) const
{
  MBErrorCode rval;
  MBRange::const_pair_iterator i;
  for (i = entities.const_pair_begin(); i != entities.const_pair_end(); ++i) {
    const MBEntityType type1 = TYPE_FROM_HANDLE(i->first);
    const MBEntityType type2 = TYPE_FROM_HANDLE(i->second);
    if (type1 == type2) {
      rval = typeData[type1].check_valid_handles( i->first, i->second );
      if (MB_SUCCESS != rval)
        return rval;
    }
    else {
      int junk;
      MBEntityHandle split = CREATE_HANDLE( type2, 0, junk );
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

MBErrorCode SequenceManager::delete_entity( MBEntityHandle entity )
{
  return typeData[TYPE_FROM_HANDLE(entity)].erase( entity );
}

MBErrorCode SequenceManager::delete_entities( const MBRange& entities )
{
  MBErrorCode rval = check_valid_entities( entities );
  if (MB_SUCCESS != rval)
    return rval;
  
  MBErrorCode result = MB_SUCCESS;
  MBRange::const_pair_iterator i;
  for (i = entities.const_pair_begin(); i != entities.const_pair_end(); ++i) {
    const MBEntityType type1 = TYPE_FROM_HANDLE(i->first);
    const MBEntityType type2 = TYPE_FROM_HANDLE(i->second);
    if (type1 == type2) {
      rval = typeData[type1].erase( i->first, i->second );
      if (MB_SUCCESS != rval)
        return result = rval;
    }
    else {
      int junk;
      MBEntityHandle split = CREATE_HANDLE( type2, 0, junk );
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
  
MBErrorCode SequenceManager::create_vertex( unsigned proc_id,
                                            const double coords[3],
                                            MBEntityHandle& handle )
{
  if (proc_id == (unsigned)-1)
    proc_id = handleUtils.proc_rank();

  const MBEntityHandle start = handleUtils.create_handle( MBVERTEX, handleUtils.first_id(proc_id), proc_id );
  const MBEntityHandle   end = handleUtils.create_handle( MBVERTEX, handleUtils. last_id(proc_id), proc_id );
  bool append;
  TypeSequenceManager::iterator seq = typeData[MBVERTEX].find_free_handle( start, end, append );
  VertexSequence* vseq;
  
  if (seq == typeData[MBVERTEX].end()) {
    SequenceData* seq_data = 0;
    handle = typeData[MBVERTEX].find_free_sequence( DEFAULT_VERTEX_SEQUENCE_SIZE, start, end, seq_data );
    if (!handle) 
      return MB_FAILURE;
    
    if (seq_data) 
      vseq = new VertexSequence( handle, 1, seq_data );
    else
      vseq = new VertexSequence( handle, 1, DEFAULT_VERTEX_SEQUENCE_SIZE );
      
    MBErrorCode rval = typeData[MBVERTEX].insert_sequence( vseq );
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

  
MBErrorCode SequenceManager::create_element( MBEntityType type,
                                             unsigned proc_id,
                                             const MBEntityHandle* conn,
                                             unsigned conn_len,
                                             MBEntityHandle& handle )
{
  if (proc_id == (unsigned)-1)
    proc_id = handleUtils.proc_rank();

  if (type <= MBVERTEX || type >= MBENTITYSET)
    return MB_TYPE_OUT_OF_RANGE;
  
  const MBEntityHandle start = handleUtils.create_handle( type, handleUtils.first_id(proc_id), proc_id );
  const MBEntityHandle   end = handleUtils.create_handle( type, handleUtils. last_id(proc_id), proc_id );
  bool append;
  TypeSequenceManager::iterator seq = typeData[type].find_free_handle( start, end, append, conn_len );
  UnstructuredElemSeq* eseq;
  
  if (seq == typeData[type].end()) {
    SequenceData* seq_data = 0;
    unsigned size = DEFAULT_ELEMENT_SEQUENCE_SIZE;
    if (type == MBPOLYGON || type == MBPOLYHEDRON) {
      size = DEFAULT_POLY_SEQUENCE_SIZE / conn_len;
      if (!size)
        size = 1;
    }
    
    handle = typeData[type].find_free_sequence( size, start, end, seq_data, conn_len );
    if (!handle) 
      return MB_FAILURE;
    
    if (seq_data) 
      eseq = new UnstructuredElemSeq( handle, 1, conn_len, seq_data );
    else
      eseq = new UnstructuredElemSeq( handle, 1, conn_len, size );
      
    MBErrorCode rval = typeData[type].insert_sequence( eseq );
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


  
MBErrorCode SequenceManager::create_mesh_set( unsigned proc_id,
                                              unsigned flags,
                                              MBEntityHandle& handle )
{
  if (proc_id == (unsigned)-1)
    proc_id = handleUtils.proc_rank();

  const MBEntityHandle start = handleUtils.create_handle( MBENTITYSET, handleUtils.first_id(proc_id), proc_id );
  const MBEntityHandle   end = handleUtils.create_handle( MBENTITYSET, handleUtils. last_id(proc_id), proc_id );
  bool append;
  TypeSequenceManager::iterator seq = typeData[MBENTITYSET].find_free_handle( start, end, append );
  MeshSetSequence* msseq;
  
  if (seq == typeData[MBENTITYSET].end()) {
    SequenceData* seq_data = 0;
    handle = typeData[MBENTITYSET].find_free_sequence( DEFAULT_MESHSET_SEQUENCE_SIZE, start, end, seq_data );
    if (!handle) 
      return MB_FAILURE;
    
    if (seq_data) 
      msseq = new MeshSetSequence( handle, 1, flags, seq_data );
    else
      msseq = new MeshSetSequence( handle, 1, flags, DEFAULT_MESHSET_SEQUENCE_SIZE );
      
    MBErrorCode rval = typeData[MBENTITYSET].insert_sequence( msseq );
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

MBErrorCode SequenceManager::allocate_mesh_set( MBEntityHandle handle,
                                                unsigned flags )
{
  SequenceData* data = 0;
  TypeSequenceManager::iterator seqptr; 
  MBEntityHandle block_start = 1, block_end = 0;
  MBErrorCode rval = typeData[MBENTITYSET].is_free_handle( handle, seqptr, data, block_start, block_end );
  if (MB_SUCCESS != rval)
    return rval;
  
  MeshSetSequence* seq;
  if (seqptr != typeData[MBENTITYSET].end()) {
    seq = static_cast<MeshSetSequence*>(*seqptr);
    if (seq->start_handle() + 1 == handle) {
      rval = seq->push_front( 1, &flags );
      if (MB_SUCCESS == rval) {
        rval = typeData[MBENTITYSET].notify_prepended( seqptr );
        if (MB_SUCCESS != rval)
          seq->pop_front( 1 );
      }
      return rval;
    }
    else if (seq->end_handle() == handle + 1) {
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
    
    MBErrorCode rval = typeData[MBENTITYSET].insert_sequence( seq );
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
SequenceManager::trim_sequence_block( MBEntityHandle start_handle,
                                      MBEntityHandle& end_handle,
                                      unsigned max_size )
{
  assert( end_handle >= start_handle );
  assert( (int)max_size > 0 ); // cast to int also prohibits some rediculously large values
  
    // if input range is larger than preferred size, trim it
  if (end_handle - start_handle >= max_size)
    end_handle = start_handle + max_size - 1;

    // if range spans more than one proc, trim it to one
  const unsigned rank = handleUtils.rank_from_handle( start_handle );
  if (handleUtils.rank_from_handle(end_handle) != rank) 
    end_handle = handleUtils.create_handle( TYPE_FROM_HANDLE(start_handle),
                                            handleUtils.max_id(),
                                            rank );
}

MBEntityHandle 
SequenceManager::sequence_start_handle( MBEntityType type,
                                        MBEntityID count,
                                        int size,
                                        MBEntityID start,
                                        int proc,
                                        SequenceData*& data )
{
  TypeSequenceManager &tsm = typeData[type];
  data = 0;
  MBEntityHandle handle = handleUtils.create_handle( type, start, proc );
  if (start < MB_START_ID ||
      !tsm.is_free_sequence( handle, count, data, size )) {
    MBEntityHandle pstart = handleUtils.create_handle( type, MB_START_ID, proc );
    MBEntityHandle pend   = handleUtils.create_handle( type, MB_END_ID,  proc );
    handle = tsm.find_free_sequence( count, pstart, pend, data, size );
  }
  return handle;
}

MBErrorCode 
SequenceManager::create_entity_sequence( MBEntityType type,
                                         MBEntityID count,
                                         int size,
                                         MBEntityID start,
                                         int proc,
                                         MBEntityHandle& handle,
                                         EntitySequence*& sequence )
{
  if (proc == -1)
    proc = handleUtils.proc_rank();

  SequenceData* data = 0;
  handle = sequence_start_handle( type, count, size, start, proc, data );
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
    else 
      sequence = new VertexSequence( handle, count, count );
    
    break;
  
  default:
    if (size == 0)
      return MB_INDEX_OUT_OF_RANGE;

    if (data)
      sequence = new UnstructuredElemSeq( handle, count, size, data );
    else 
      sequence = new UnstructuredElemSeq( handle, count, size, count );

    break;
  }
  
  MBErrorCode result = typeData[type].insert_sequence( sequence );
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


MBErrorCode 
SequenceManager::create_meshset_sequence( MBEntityID count,
                                          MBEntityID start,
                                          int proc,
                                          const unsigned* flags,
                                          MBEntityHandle& handle,
                                          EntitySequence*& sequence )
{
  if (proc == -1)
    proc = handleUtils.proc_rank();

  SequenceData* data = 0;
  handle = sequence_start_handle( MBENTITYSET, count, 0, start, proc, data );
  if (!handle)
    return MB_MEMORY_ALLOCATION_FAILED;
  
  if (data)
    sequence = new MeshSetSequence( handle, count, flags, data );
  else
    sequence = new MeshSetSequence( handle, count, flags, count );
  
  
  
  MBErrorCode result = typeData[MBENTITYSET].insert_sequence( sequence );
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


MBErrorCode 
SequenceManager::create_meshset_sequence( MBEntityID count,
                                          MBEntityID start,
                                          int proc,
                                          unsigned flags,
                                          MBEntityHandle& handle,
                                          EntitySequence*& sequence )
{
  if (proc == -1)
    proc = handleUtils.proc_rank();

  SequenceData* data = 0;
  handle = sequence_start_handle( MBENTITYSET, count, 0, start, proc, data );
  if (!handle)
    return MB_MEMORY_ALLOCATION_FAILED;
  
  if (data)
    sequence = new MeshSetSequence( handle, count, flags, data );
  else
    sequence = new MeshSetSequence( handle, count, flags, count );
  
  
  
  MBErrorCode result = typeData[MBENTITYSET].insert_sequence( sequence );
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

MBErrorCode
SequenceManager::create_scd_sequence( int imin, int jmin, int kmin,
                                      int imax, int jmax, int kmax,
                                      MBEntityType type,
                                      MBEntityID start_id_hint,
                                      int processor_id,
                                      MBEntityHandle& handle,
                                      EntitySequence*& sequence )
{
  if (processor_id == -1)
    processor_id = handleUtils.proc_rank();

  int this_dim = MBCN::Dimension(type);

    // use > instead of != in the following assert to also catch cases where imin > imax, etc.
  assert((this_dim < 3 || kmax > kmin) &&
         (this_dim < 2 || jmax > jmin) &&
         (this_dim < 1 || imax > imin));

    // compute # entities; not as easy as it would appear...
  MBEntityID num_ent;
  if (MBVERTEX == type)
    num_ent = (MBEntityID)(imax-imin+1)*(MBEntityID)(jmax-jmin+1)*(MBEntityID)(kmax-kmin+1);
  else {
    num_ent = (imax-imin) *
      (this_dim >= 2 ? (jmax-jmin) : 1) *
      (this_dim >= 3 ? (kmax-kmin) : 1);
  }
  
    // get a start handle
  SequenceData* data = 0;
  handle = sequence_start_handle( type, num_ent, -1, start_id_hint, processor_id, data );
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
  
  MBErrorCode result = typeData[type].insert_sequence( sequence );
  if (MB_SUCCESS != result) {
    data = sequence->data();
    delete sequence;
    delete data;
    return result;
  }
  
  return MB_SUCCESS;
}

MBErrorCode
SequenceManager::create_scd_sequence( const HomCoord& coord_min,
                                      const HomCoord& coord_max,
                                      MBEntityType type,
                                      MBEntityID start_id_hint,
                                      int processor_id,
                                      MBEntityHandle& first_handle_out,
                                      EntitySequence*& sequence_out )
{
  return create_scd_sequence( coord_min.i(), coord_min.j(), coord_min.k(),
                              coord_max.i(), coord_max.j(), coord_max.k(),
                              type, start_id_hint,  processor_id,
                              first_handle_out, sequence_out );
}
 
MBErrorCode
SequenceManager::replace_subsequence( EntitySequence* new_seq, TagServer* ts )
{
  const MBEntityType type = TYPE_FROM_HANDLE(new_seq->start_handle());
  return typeData[type].replace_subsequence( new_seq, ts );
}

void SequenceManager::get_memory_use( unsigned long& total_entity_storage,
                                      unsigned long& total_storage ) const

{
  total_entity_storage = 0;
  total_storage = 0;
  unsigned long temp_entity, temp_total;
  for (MBEntityType i = MBVERTEX; i < MBMAXTYPE; ++i) {
    temp_entity = temp_total = 0;
    get_memory_use( i, temp_entity, temp_total );
    total_entity_storage += temp_entity;
    total_storage        += temp_total;
  }
}

void SequenceManager::get_memory_use( MBEntityType type,
                                      unsigned long& total_entity_storage,
                                      unsigned long& total_storage ) const
{
  typeData[type].get_memory_use( total_entity_storage, total_storage );
}

void SequenceManager::get_memory_use( const MBRange& entities,
                                      unsigned long& total_entity_storage,
                                      unsigned long& total_amortized_storage ) const
{
  total_entity_storage = 0;
  total_amortized_storage = 0;
  unsigned long temp_entity, temp_total;
  MBRange::const_pair_iterator i;
  for (i = entities.const_pair_begin(); i != entities.const_pair_end(); ++i) {
    const MBEntityType t1 = TYPE_FROM_HANDLE(i->first);
    const MBEntityType t2 = TYPE_FROM_HANDLE(i->second);
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

// These are meant to be called from the debugger (not declared in any header)
// so leave them out of release builds (-DNDEBUG).
#ifndef NDEBUG
#include <iostream>

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
  for (MBEntityType t = MBVERTEX; t < MBMAXTYPE; ++t) 
    if (!seq_man.entity_map(t).empty()) 
      s << std::endl 
        << "****************** " << MBCN::EntityTypeName( t ) << " ******************"
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
