/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

#include "MBProcConfig.hpp"
#include "MBInternals.hpp"


/** Calculate ceiling of log 2 of a positive integer */
static unsigned ceil_log_2( unsigned n )
{
  unsigned result;
  for (result = 0; n > (((MBEntityHandle)1)<<result); ++result);
  return result;
}

//! Constructor
MBProcConfig::MBProcConfig( unsigned rank, unsigned size ) 
  : procRank( rank ),
    procSize( size ),
    procWidth( ceil_log_2( size ) ),
    idWidth( MB_ID_WIDTH - procWidth ),
    idMask( MB_ID_MASK >> procWidth ),
    procMask( ~(MB_TYPE_MASK|idMask) )
{}

MBEntityHandle MBProcConfig::handle( MBEntityType type, 
                                     unsigned sub_id, 
                                     unsigned proc ) const
{
  int junk;
  return CREATE_HANDLE( type, id( sub_id, proc ), junk );
}
