#/**
# * MOAB, a Mesh-Oriented datABase, is a software component for creating,
# * storing and accessing finite element mesh data.
# * 
# * Copyright 2004 Sandia Corporation.  Under the terms of Contract
# * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
# * retains certain rights in this software.
# * 
# * This library is free software; you can redistribute it and/or
# * modify it under the terms of the GNU Lesser General Public
# * License as published by the Free Software Foundation; either
# * version 2.1 of the License, or (at your option) any later version.
# * 
# */

MOAB_DIR = ${HOME}/MOAB
CHACO_DIR = ${HOME}/Chaco-2.2

include ${MOAB_DIR}/MB.lin
include ${MOAB_DIR}/MB.local

INCLUDES = -I${CHACO_DIR}/code/main -I${MOAB_DIR}

mbchaco:	main.o ComputePartition.o makefile
		${CXX} ${DEBUG_FLAG} main.o ComputePartition.o -L${MOAB_DIR} -lMOAB ${NETCDF_LIBS} -L${CHACO_DIR}/code -lchaco -lm -o mbchaco

lint:
		lint ${IFLAG} ${FILES.c} -lm

alint:
		alint ${IFLAG} ${FILES.c} -lm

clean:
		rm -f */*.o ${DEST_DIR}/core

.cpp.o:
		${CXX} -c ${DEBUG_FLAG} ${INCLUDES} -o $*.o $*.cpp
