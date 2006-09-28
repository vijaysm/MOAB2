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
  int err;
  return CREATE_HANDLE( type, id( sub_id, proc ), err );
}

MBRange::const_iterator 
MBProcConfig::lower_bound( MBEntityType type, unsigned proc, const MBRange& range ) const
{
  int err;
  MBEntityHandle h = CREATE_HANDLE( type, id(0,proc), err );
  return err ? range.end() : MBRange::lower_bound(range.begin(), range.end(), h);
}

MBRange::const_iterator
MBProcConfig::upper_bound( MBEntityType type, unsigned proc, const MBRange& range ) const
{
  int err;
  MBEntityHandle h = CREATE_HANDLE( type, last_id(proc), err );
  return err ? range.end() : MBRange::upper_bound(range.begin(), range.end(), h);
}

std::pair<MBRange::const_iterator, MBRange::const_iterator>
MBProcConfig::equal_range( MBEntityType type, unsigned proc, const MBRange& range ) const
{
  std::pair<MBRange::const_iterator, MBRange::const_iterator> iters;
  int err;
  MBEntityHandle h;

  h = CREATE_HANDLE( type, id(0,proc), err );
  iters.first = err ? range.end() : MBRange::lower_bound(range.begin(), range.end(), h);  
  
  h = CREATE_HANDLE( type, last_id(proc), err );
  iters.second = err ? range.end() : MBRange::upper_bound( iters.first, range.end(), h );
  
  return iters;
}

MBRange MBProcConfig::subset( unsigned proc, const MBRange& range ) const
{
  int junk;
  MBRange result;
  MBRange::iterator insert_pos = result.begin();
  MBRange::const_pair_iterator iter;
  MBEntityHandle s, e;
  
  for (iter = range.const_pair_begin(); iter != range.const_pair_end(); ++iter)
  {
    const MBEntityType beg_type = TYPE_FROM_HANDLE(iter->first),
                       end_type = TYPE_FROM_HANDLE(iter->second);
    const unsigned beg_rank = rank(iter->first), end_rank = rank(iter->second);
    
    if (beg_type != end_type) {
      if (beg_rank <= proc) {
        s = beg_rank == proc ? iter->first : 
            CREATE_HANDLE( beg_type,    id(0,proc), junk );
        e = CREATE_HANDLE( beg_type, last_id(proc), junk );
        insert_pos = result.insert( insert_pos, s, e );
      }
      MBEntityType t = beg_type;
      for (++t; t != end_type; ++t) {
        s = CREATE_HANDLE( t,    id(0,proc), junk );
        e = CREATE_HANDLE( t, last_id(proc), junk );
        insert_pos = result.insert( insert_pos, s, e );
      }
      if (end_rank >= proc) {
        e = end_rank == proc ? iter->second :
            CREATE_HANDLE( end_type, last_id(proc), junk );
        s = CREATE_HANDLE( end_type,    id(0,proc), junk );
        insert_pos = result.insert( insert_pos, s, e );
      }
    }
    else if (beg_rank <= proc && end_rank >= proc) {
      s = (beg_rank == proc) ? iter->first  : CREATE_HANDLE( beg_type,    id(0,proc), junk );
      e = (end_rank == proc) ? iter->second : CREATE_HANDLE( beg_type, last_id(proc), junk );
      insert_pos = result.insert( insert_pos, s, e );
    }
  }
  
  return result;
}


