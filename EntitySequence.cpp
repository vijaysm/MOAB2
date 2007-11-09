#include "EntitySequence.hpp"
#include "SequenceData.hpp"

bool EntitySequence::using_entire_data() const 
{
  return start_handle() == data()->start_handle()
        && end_handle() == data()->  end_handle();
}

int EntitySequence::values_per_entity() const
  { return 0; }

MBErrorCode EntitySequence::pop_back( MBEntityID count )
{
  MBEntityHandle new_end = endHandle - count;
  if (new_end < startHandle)
    return MB_FAILURE;
  
  endHandle = new_end;
  return MB_SUCCESS;
}

MBErrorCode EntitySequence::pop_front( MBEntityID count )
{
  MBEntityHandle new_start = startHandle + count;
  if (new_start > endHandle)
    return MB_FAILURE;
  
  startHandle = new_start;
  return MB_SUCCESS;
}


MBErrorCode EntitySequence::prepend_entities( MBEntityID count )
{
  MBEntityHandle new_start = startHandle - count;
  if (new_start < data()->start_handle())
    return MB_FAILURE;
  
  startHandle = new_start;
  return MB_SUCCESS;
}
 
MBErrorCode EntitySequence::append_entities( MBEntityID count )
{
  MBEntityHandle new_end = endHandle + count;
  if (new_end > data()->end_handle())
    return MB_FAILURE;
  
  endHandle = new_end;
  return MB_SUCCESS;
}

MBErrorCode EntitySequence::merge( EntitySequence& other )
{
  if (sequenceData != other.sequenceData)
    return MB_FAILURE;
  if (end_handle() + 1 == other.start_handle()) {
    endHandle = other.end_handle();
    other.startHandle = other.end_handle()+1;
  }
  else if (start_handle() == other.end_handle() + 1) {
    startHandle = other.start_handle();
    other.endHandle = other.start_handle()-1;
  }
  else
    return MB_FAILURE;
  return MB_SUCCESS;
}

unsigned long EntitySequence::get_per_entity_memory_use( MBEntityHandle,
                                                         MBEntityHandle ) const
  { return 0; }
