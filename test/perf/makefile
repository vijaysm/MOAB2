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
SIDL_MESH_DIR = ${MOAB_DIR}/TSTT
SIDLDIR = /usr/local/babel-0.8.8

include ${MOAB_DIR}/MB.lin
include ${MOAB_DIR}/MB.local

#DEBUG_FLAG = -g
DEBUG_FLAG = -O2
#DEBUG_FLAG = -pg -g 

# define REALTFI if you want code to compute node locations using TFI, otherwise xyz is just
# a linear interpolation of that coordinate over domain
#CXXFLAGS = -DREALTFI

LD = g++ ${DEBUG_FLAG}
CGM_SRC = ${HOME}/cubit/clean
CUBIT = /cubit

default: cubit_perf
all: cubit_perf perf tstt_perf

perf: perf.cpp
	g++ $(CXXFLAGS) $(DEBUG_FLAG) -o perf -DLINUX -DREALTFI -DIS_BUILDING_MB -I ../.. perf.cpp -L ../.. -lMOAB -L/cubit/netcdf/lib -lnetcdf_c++ -lnetcdf

cubit_perf: cubit_perf.o
	g++ ${DEBUG_FLAG} -I${CGM_SRC}/mdb -I${CGM_SRC}/util cubit_perf.o $(CGM_SRC)/libcubiti.so.1.0.1 -o cubit_perf -L$(CGM_SRC)/showviz/lib -lshowviz_geom_mesh -lshowviz_mesh -lshowviz_geom -lshowviz \
-L/cubit/VTK/VTK-4.2.2/lib/vtk -lvtkHybrid -lvtkIO -lvtkRendering -lvtkftgl -lvtkfreetype -lvtkImaging -lvtkGraphics -lvtkFiltering -lvtkCommon -lvtkexpat -lvtkjpeg -lvtkpng -lvtktiff -lvtkzlib -L$(CGM_SRC)/df -lcbt_df  -L$(CGM_SRC)/mdb -lcubit_mdb -L$(CGM_SRC)/geom/virtual -lcubit_virtual -L$(CGM_SRC)/geom/facet -lcubit_facet -L$(CGM_SRC)/geom/facetbool -lcubit_facetboolstub -L$(CGM_SRC)/geom -lcubit_geom -L$(CGM_SRC)/geom/Cholla -lcholla -L$(CGM_SRC)/util -lcubit_util -L$(CGM_SRC)/camal/lib/Linux -lCMLTri102 \
        -L$(CGM_SRC)/geom/ACIS -lcubit_ACIS -L/cubit/acis/acis13.2/bin/linux_so -lxstep -lxacis2k -lxcore2k -lxutil -lacisstep -lxiges -lacisiges -lSpaAWarp -lSpaAVis -lSpaASurf -lSpaALops -lSpaABlend -lSpaACIS -lSpaBase -L/cubit/lib -llpk -L/cubit/lib -lexoIIv2c326 -L/cubit/netcdf/lib -lnetcdf  -L/cubit/verdict/verdict1.0.1/lib -lverdict101-3.2  -L/usr/X11R6/lib -lXt -lXext -lX11 

cubit_perf.o: cubit_perf.cpp
	g++ -c -DLINUX -DREALTFI ${DEBUG_FLAG} -I${CGM_SRC}/mdb -I${CGM_SRC}/util cubit_perf.cpp

tstt_perf.o: tstt_perf.cpp
	g++ -c ${DEBUG_FLAG} -DLINUX -DREALTFI -I${SIDL_MESH_DIR} -I$(SIDLDIR)/include -I${SIDL_MESH_DIR}/server tstt_perf.cpp

tstt_perf: tstt_perf.o
	$(LD)  -o $@ ${DEBUG_FLAG} tstt_perf.o -Wl,-rpath,${SIDL_MESH_DIR} -Wl,-rpath,${MOAB_DIR} -L${SIDL_MESH_DIR} -lMOABserver -L$(SIDLDIR)/lib -lsidl -L${MOAB_DIR} -lMOAB ${NETCDF_LIBS}


