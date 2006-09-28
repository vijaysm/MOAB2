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

#ifndef MB_PROC_CONFIG_HPP
#define MB_PROC_CONFIG_HPP

#include "MBTypes.h"

/**\brief Multi-CPU information for parallel MOAB */
class MBProcConfig {
  public:

    MBProcConfig( unsigned rank, unsigned size );
    
      //! Get the current processor number
    unsigned rank() const 
      { return procRank; }
      
      //! Get the number of processors
    unsigned size() const 
      { return procSize; }
      
      //! Get CPU number from handle
    unsigned rank( MBEntityHandle handle ) const
      { return (handle & procMask) >> idWidth; }
      
      //! Get CPU number from ID
    unsigned rank_from_id( unsigned id ) const
      { return id >> idWidth; }
      
      //! Get maximum entity ID that can be stored in a
      //! a handle, allowing for the processor number
    unsigned max_id() const
      { return idMask; }
      
      //! Create the ID portion of a handle by combining
      //! an actual ID and a processor number
    MBEntityHandle id( unsigned sub_id, unsigned proc ) const
      { return (proc << idWidth) | sub_id; }
      
      //! Extract non-rank portion of entity ID from handle
    MBEntityHandle id( MBEntityHandle h ) const
      { return h & idMask; }
      
    MBEntityHandle first_id( unsigned proc ) const
      { return id( 1, proc ); }
    
    MBEntityHandle last_id( unsigned proc ) const
      { return id( max_id(), proc ); }
      
      //! Create an entity handle given type, rank, and id
    MBEntityHandle handle( MBEntityType type, 
                           unsigned sub_id, 
                           unsigned proc ) const;
  private:
  
    unsigned procRank;    //!< ID of this processor
    unsigned procSize;    //!< Total number of processors
    unsigned procWidth;   //!< Number of bits in handle for processor ID
    unsigned idWidth;     //!< Number of bits in handle for entity ID
    MBEntityHandle idMask;
    MBEntityHandle procMask;
};

#endif
