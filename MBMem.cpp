

/*  File      :   MBMem.cpp
 *  Purpose   :   Track Memory usage and memory debugging in MB
 *  Creator   :   Clinton Stimpson
 *  Date      :   28 Aug 2002
*/


/* Settings:
 *
 * MB_MEM_DEBUG   -- If defined, will dump MB.log file with memory errors.
 *
*/


#include "MBMem.hpp"
#include "MBInterface.hpp"


// Allocator for adjacency vectors in AEntityFactory
DEFINE_MB_ALLOCATOR_CLASS( MBAllocator, MBEntityHandle )


