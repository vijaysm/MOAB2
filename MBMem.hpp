
#ifndef MB_MEM_HPP
#define MB_MEM_HPP

/*  File      :   MBMem.hpp
 *  Purpose   :   Track Memory usage and memory debugging in MB
 *  Creator   :   Clinton Stimpson
 *  Date      :   28 Aug 2002
*/


/* Settings:
 *
 * MB_MEM_DEBUG   -- If defined, will dump MB.log file with memory errors.
 *
*/

#include "MBAlloc.hpp"
#include "MBInterface.hpp"

// define the allocators to use in various parts of MB

// Allocator for adjacency vectors in AEntityFactory
DECLARE_MB_ALLOCATOR_CLASS( MBAllocator, MBEntityHandle, MBAdjacencyVectorAllocator )
// the type of vector to use for storing adjacency information
typedef std::vector<MBEntityHandle, MBAdjacencyVectorAllocator > MBAdjacencyVector;


#endif  // MB_MEM_HPP
