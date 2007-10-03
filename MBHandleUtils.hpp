#ifndef MBHANDLEUTILS_HPP
#define MBHANDLEUTILS_HPP

#include "MBInterface.hpp"
#include "MBRange.hpp"

class MBHandleUtils 
{
public:
  MBHandleUtils(int proc_rank, int proc_size);

    //! Get processor rank
  unsigned proc_rank() const {return procRank;}
  
    //! Get processor size
  unsigned proc_size() const {return procSize;}
      
    //! Get CPU number from handle
  unsigned rank_from_handle( MBEntityHandle handle ) const
    { return (handle & procMask) >> idWidth; }
      
    //! Get CPU number from ID
  unsigned rank_from_id( MBEntityID id ) const
    { return id >> idWidth; }
      
    //! Get maximum entity ID that can be stored in a
    //! a handle, allowing for the processor number
  MBEntityID max_id() const
    { return idMask; }
      
    //! Create the ID portion of a handle by combining
    //! an actual ID and a processor number
  MBEntityID create_id( MBEntityID sub_id, unsigned proc ) const
    { return ((MBEntityHandle)proc << idWidth) | (MBEntityHandle)sub_id; }
      
    //! Extract non-rank portion of entity ID from handle
  MBEntityID id_from_handle( MBEntityHandle h ) const
    { return h & idMask; }
      
  MBEntityID first_id( unsigned proc ) const
    { return create_id( 1, proc ); }
    
  MBEntityID last_id( unsigned proc ) const
    { return create_id( max_id(), proc ); }
      
    //! Create an entity handle given type, rank, and id
  MBEntityHandle create_handle( MBEntityType type, 
                                MBEntityID sub_id, 
                                unsigned proc ) const;

    //! return a subset with corresponding proc values in handles
  MBRange subset_by_proc( unsigned proc, 
                          const MBRange &range) const;
  
    //! return a lower bound for handles with corresponding proc values
  MBRange::const_iterator lower_bound( MBEntityType type, 
                                       unsigned proc, 
                                       const MBRange &range) const;
  
    //! return an upper bound for handles with corresponding proc values
  MBRange::const_iterator upper_bound( MBEntityType type, 
                                       unsigned proc, 
                                       const MBRange &range) const;
  
    //! return an equal range for handles with corresponding proc values
  std::pair<MBRange::const_iterator,MBRange::const_iterator>
  equal_range( MBEntityType type, unsigned proc, 
               const MBRange &range) const;
  
private:
/** Calculate ceiling of log 2 of a positive integer */
  static unsigned ceil_log_2( unsigned n );

  unsigned procRank;    //!< ID of this processor
  unsigned procSize;    //!< Total number of processors
  unsigned procWidth;   //!< Number of bits in handle for processor ID
  unsigned idWidth;     //!< Number of bits in handle for entity ID
  MBEntityHandle idMask;
  MBEntityHandle procMask;
  
};

inline unsigned MBHandleUtils::ceil_log_2( unsigned n )
{
  unsigned result;
  for (result = 0; n > (((MBEntityHandle)1)<<result); ++result);
  return result;
}

#endif
