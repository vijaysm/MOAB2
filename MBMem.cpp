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


