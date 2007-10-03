#include "MBHandleUtils.hpp"
#include "MBInternals.hpp"
#include "assert.h"

MBHandleUtils::MBHandleUtils(int proc_rank, int proc_size) : 
    procRank((int) proc_rank),
    procSize((int) proc_size),
    procWidth( ceil_log_2(procSize)),
    idWidth( MB_ID_WIDTH - procWidth ),
    idMask( MB_ID_MASK >> procWidth ),
    procMask( ~(MB_TYPE_MASK|idMask) )
{
  assert(0 <= proc_rank && 0 <= proc_size);
}

MBEntityHandle MBHandleUtils::create_handle( MBEntityType type, 
                                             MBEntityID sub_id, 
                                             unsigned proc ) const
{
  int err;
  return CREATE_HANDLE( type, create_id( sub_id, proc ), err );
}

MBRange MBHandleUtils::subset_by_proc( unsigned proc, 
                                       const MBRange &range) const
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
    const unsigned beg_rank = rank_from_handle(iter->first), 
      end_rank = rank_from_handle(iter->second);
    
    if (beg_type != end_type) {
      if (beg_rank <= proc) {
        s = beg_rank == proc ? iter->first : 
            CREATE_HANDLE( beg_type,    create_id(0,proc), junk );
        e = CREATE_HANDLE( beg_type, last_id(proc), junk );
        insert_pos = result.insert( insert_pos, s, e );
      }
      MBEntityType t = beg_type;
      for (++t; t != end_type; ++t) {
        s = CREATE_HANDLE( t,    create_id(0,proc), junk );
        e = CREATE_HANDLE( t, last_id(proc), junk );
        insert_pos = result.insert( insert_pos, s, e );
      }
      if (end_rank >= proc) {
        e = end_rank == proc ? iter->second :
            CREATE_HANDLE( end_type, last_id(proc), junk );
        s = CREATE_HANDLE( end_type,    create_id(0,proc), junk );
        insert_pos = result.insert( insert_pos, s, e );
      }
    }
    else if (beg_rank <= proc && end_rank >= proc) {
      s = (beg_rank == proc) ? iter->first  : 
        CREATE_HANDLE( beg_type, create_id(0,proc), junk );
      e = (end_rank == proc) ? iter->second : 
        CREATE_HANDLE( beg_type, last_id(proc), junk );
      insert_pos = result.insert( insert_pos, s, e );
    }
  }
  
  return result;
}

MBRange::const_iterator MBHandleUtils::lower_bound( MBEntityType type, 
                                                    unsigned proc, 
                                                    const MBRange &range) const
{
  int err;
  MBEntityHandle h = CREATE_HANDLE( type, create_id(0,proc), err );
  return err ? range.end() : 
    MBRange::lower_bound(range.begin(), range.end(), h);
}

MBRange::const_iterator MBHandleUtils::upper_bound( MBEntityType type, 
                                                    unsigned proc, 
                                                    const MBRange &range) const
{
  int err;
  MBEntityHandle h = CREATE_HANDLE( type, last_id(proc), err );
  return err ? range.end() : 
    MBRange::upper_bound(range.begin(), range.end(), h);
}

std::pair<MBRange::const_iterator,MBRange::const_iterator>
MBHandleUtils::equal_range( MBEntityType type, unsigned proc, 
                            const MBRange &range) const
{
  std::pair<MBRange::const_iterator, MBRange::const_iterator> iters;
  int err;
  MBEntityHandle h;

  h = CREATE_HANDLE( type, create_id(0,proc), err );
  iters.first = err ? range.end() : 
    MBRange::lower_bound(range.begin(), range.end(), h);  
  
  h = CREATE_HANDLE( type, last_id(proc), err );
  iters.second = err ? range.end() : 
    MBRange::upper_bound( iters.first, range.end(), h );
  
  return iters;
}

