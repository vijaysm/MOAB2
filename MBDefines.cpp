

#include "MBInternals.hpp"

//! define non-inline versions of these functions for debugging
MBEntityHandle ifh(MBEntityHandle handle) 
{
  return ID_FROM_HANDLE(handle);
}

MBEntityType tfh(MBEntityHandle handle)
{
  return TYPE_FROM_HANDLE(handle);
}






