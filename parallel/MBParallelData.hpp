/**
 * \class MBParallelData
 * \brief Parallel data in MOAB
 * \author Tim Tautges
 *
 *  This class implements methods to retrieve information about 
 * the parallel mesh from MOAB.  Most of this data can be retrieved
 * directly from MOAB as sets and tags; this class provides convenience
 * methods implemented on top of other MOAB functions.
 *
 */

#ifndef MB_PARALLEL_DATA_HPP
#define MB_PARALLEL_DATA_HPP

#include "MBForward.hpp"
#include "MBRange.hpp"

class MBParallelComm;

class MBParallelData
{
public:

    //! constructor; if non-null parallelcomm, that is used to
    //! determine rank, otherwise rank is taken from impl
  MBParallelData(MBInterface *impl, MBParallelComm *pcomm = NULL);

    //! return partition sets; if tag_name is input, gets sets with
    //! that tag name, otherwise uses PARALLEL_PARTITION tag
  MBErrorCode get_partition_sets(MBRange &part_sets,
                                 const char *tag_name = NULL);

    //! get communication interface sets and the processors with which
    //! this processor communicates; sets are sorted by processor
  MBErrorCode get_interface_sets(std::vector<MBEntityHandle> &iface_sets,
                                 std::vector<int> &iface_procs);
  

private:

    //! interface instance to which this instance corresponds
  MBInterface *mbImpl;

    //! MBParallelComm object to which this is bound
  MBParallelComm *parallelComm;
  
};

inline MBParallelData::MBParallelData(MBInterface *impl, 
                                      MBParallelComm *pcomm) 
    : mbImpl(impl), parallelComm(pcomm) 
{}

#endif
