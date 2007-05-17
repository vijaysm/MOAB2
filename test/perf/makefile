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
TSTT_DIR = ${HOME}/tstt/TSTTbind

include ${TSTT_DIR}/TSTTM/TSTTM-Defs.inc

# define REALTFI if you want code to compute node locations using TFI, otherwise xyz is just
# a linear interpolation of that coordinate over domain
CXXFLAGS = -DREALTFI -g
LDFLAGS = -g

CUBIT = /cubit

CUBIT_LIBS = \
	-L$(CGM_SRC)/showviz/lib -lshowviz_geom_mesh -lshowviz_mesh -lshowviz_bc -lshowviz_geom -lshowviz \
	-L${CGM_SRC} -lcubiti  \
	-L/usr/local/ParaViewComplete/bin -lvtkHybrid -lvtkIO -lvtkRendering -lvtkftgl -lvtkfreetype \
	    -lvtkImaging -lvtkGraphics -lvtkFiltering -lvtkCommon -lvtkexpat -lvtkjpeg -lvtkpng \
	    -lvtktiff -lvtkzlib \
	-L$(CGM_SRC)/df -lcbt_df  \
	-L$(CGM_SRC)/mdb -lcubit_mdb \
	-L$(CGM_SRC)/bc -lcubit_bc \
	-L${CGM_SRC}/camal/lib/Linux -lcamal_pave -lcamal_shr -lcamal_tri \
	-L${CGM_SRC}/cbt_mesquite -lcbt_mesquite \
	-L${HOME}/cubit/Mesquite_Platforms/lib -lmesquite08 \
	-L/cubit/camal/camal2.0.2/lib/Linux -lCMLTet202 \
	-L/usr/local/gcc-3.2.3/lib -lg2c \
	-L/cubit/lib -llpk \
	-L/cubit/exodus/exodus-4.01/lib -lexoIIv2c401 \
	-L/cubit/netcdf/lib -lnetcdf  \
	-L/cubit/verdict/verdict1.1.2/lib -lverdict112

CGM_LIBS = \
	-L$(CGM_SRC)/geom/virtual -lcubit_virtual \
	-L$(CGM_SRC)/geom/facet -lcubit_facet \
	-L$(CGM_SRC)/geom/facetbool -lcubit_facetboolstub \
	-L$(CGM_SRC)/geom -lcubit_geom \
	-L$(CGM_SRC)/geom/Cholla -lcholla \
	-L$(CGM_SRC)/util -lcubit_util 	\
	-L$(CGM_SRC)/geom/ACIS -lcubit_ACIS \
	-L/cubit/acis/acis13.2/bin/linux_so -lxstep -lxacis2k -lxcore2k -lxutil -lacisstep \
	    -lxiges -lacisiges -lSpaAWarp -lSpaAVis -lSpaASurf -lSpaALops -lSpaABlend -lSpaACIS -lSpaBase	

default: cubit_perf
all: cubit_perf perf tstt_perf

perf: perf.cpp
	g++ $(CXXFLAGS) $(DEBUG_FLAG) -o perf -DLINUX -DREALTFI -DIS_BUILDING_MB -I ../.. perf.cpp -L ../.. -lMOAB -L/cubit/netcdf/lib -lnetcdf_c++ -lnetcdf

cubit_perf: cubit_perf.o
	g++ ${DEBUG_FLAG} cubit_perf.o -o cubit_perf ${CUBIT_LIBS} ${CUBIT_LIBS} ${CGM_LIBS} \
        -L/usr/X11R6/lib -lXt -lXext -lX11 

cubit_perf.o: cubit_perf.cpp
	g++ -c -DLINUX -DREALTFI ${DEBUG_FLAG} -I${CGM_SRC}/mdb -I${CGM_SRC}/util cubit_perf.cpp

tstt_perf.o: tstt_perf.cpp
	${CXX} -c ${CXXFLAGS} -DREALTFI ${TSTTM_INCLUDES} tstt_perf.cpp

tstt_perf_binding.o: tstt_perf_binding.cpp
	${CXX} -c ${CXXFLAGS} -DREALTFI ${TSTTM_CBIND_SERVER_INCLUDES} tstt_perf_binding.cpp

tstt_perf: tstt_perf.o
	$(LD)  -o $@ ${DEBUG_FLAG} tstt_perf.o -Wl,-rpath,${SIDL_MESH_DIR} -Wl,-rpath,${MOAB_DIR} -L${SIDL_MESH_DIR} -lMOABserver -L$(SIDLDIR)/lib -lsidl -L${MOAB_DIR} -lMOAB ${NETCDF_LIBS}

tstt_perf_binding: tstt_perf_binding.o
	$(CXX)  -o $@ ${LDFLAGS} tstt_perf_binding.o ${TSTTM_CBIND_SERVER_LIBS}


