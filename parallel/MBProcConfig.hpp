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
#include "MBRange.hpp"

class MBInterface;


#ifdef USE_MPI
/* MPICH2 will fail if SEEK_* macros are defined
 * because they are also C++ enums. Undefine them
 * when including mpi.h and then redefine them
 * for sanity.
 */
#  ifdef SEEK_SET
#    define MB_SEEK_SET SEEK_SET
#    define MB_SEEK_CUR SEEK_CUR
#    define MB_SEEK_END SEEK_END
#    undef SEEK_SET
#    undef SEEK_CUR
#    undef SEEK_END
#  endif
#include "mpi.h"
#  ifdef MB_SEEK_SET
#    define SEEK_SET MB_SEEK_SET
#    define SEEK_CUR MB_SEEK_CUR
#    define SEEK_END MB_SEEK_END
#    undef MB_SEEK_SET
#    undef MB_SEEK_CUR
#    undef MB_SEEK_END
#  endif
extern "C" 
{
#include "types.h"
#include "errmem.h"
#include "crystal.h"
}
#else
typedef int MPI_Comm;
#define MPI_COMM_WORLD 0
typedef void* crystal_data;
#endif

/**\brief Multi-CPU information for parallel MOAB */
class MBProcConfig {
public:

  MBProcConfig(MPI_Comm proc_comm = MPI_COMM_WORLD);
  
  ~MBProcConfig();
  
    //! Get the current processor number
  unsigned proc_rank() const 
    { return procRank; }
      
    //! Get the number of processors
  unsigned proc_size() const 
    { return procSize; }
      
    //! get a crystal router for this parallel job
  crystal_data *crystal_router(bool construct_if_missing = true);

    //! get/set the communicator for this proc config
  const MPI_Comm proc_comm() const {return procComm;}
  void proc_comm(MPI_Comm this_comm) {procComm = this_comm;}
  
private:

    //! MPI communicator set for this instance
  MPI_Comm procComm;

    //! rank of this processor
  unsigned procRank;
  
    //! number of processors
  unsigned procSize;

    //! whether the crystal router's been initialized or not
  bool crystalInit;
  
    //! crystal router for this parallel job
  crystal_data crystalData;
  
};

#endif
