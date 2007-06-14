#
# Copyright 2006 Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
# retains certain rights in this software.
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
#
VPATH = mserver bserver Cclient F77client
vpath = mserver bserver Cclient F77client

include iMesh-Defs.inc

default: all

######## include make.Linux
CXX = g++
CXX_FLAGS_REQ = -fpic
CXX_FLAGS_DBG = -g -Wall
CXX_FLAGS_OPT = -O 

CC  = gcc
CC_FLAGS_REQ = -fpic
CC_FLAGS_DBG = -g -Wall
CC_FLAGS_OPT = -O

FC  = g77
#FC   = /usr/apps/pgi/3.3/linux86/bin/pgf77
FC_FLAGS_REQ = 
FC_FLAGS_DBG = -g
FC_FLAGS_OPT = -O

MKDEP = gcc -E -Wp,-MM

LD = $(CXX)
LD_FLAGS_PIC = -fpic -shared
LD_FLAGS_REQ = 
LD_FLAGS_DBG = -g
LD_FLAGS_OPT = -O

ARCHIVER = ar cr


INCPATH = -I. ${iMesh_SERVER_INCLUDES} -I$(MOAB_INCLUDEDIR)


#CXX_FLAGS = $(CXX_FLAGS_REQ) $(CXX_FLAGS_OPT) $(INCPATH)
CXX_FLAGS = $(CXX_FLAGS_REQ) $(CXX_FLAGS_DBG) $(INCPATH) ${iMesh_CLIENT_C_INCLUDES} 
#CC_FLAGS = $(CC_FLAGS_REQ) $(CC_FLAGS_OPT) $(INCPATH)
CC_FLAGS = $(CC_FLAGS_REQ) $(CC_FLAGS_DBG) $(INCPATH) ${iMesh_CLIENT_C_INCLUDES} 
FC_FLAGS = $(FC_FLAGS_REQ) $(FC_FLAGS_DBG) $(INCPATH)
LD_FLAGS = $(LD_FLAGS_REQ) $(LD_FLAGS_DBG) $(INCPATH)
#LD_FLAGS = $(LD_FLAGS_REQ) $(LD_FLAGS_OPT) $(INCPATH)
###########################################################

include mserver/babel.make
FACTORY_IMPL_HDRS := 
FACTORY_IMPL_SRCS := 
FACTORY_GEN_IOR_HDRS := 
FACTORY_GEN_IOR_SRCS := 
FACTORY_GEN_STUB_HDRS := 
FACTORY_GEN_STUB_SRCS := 
FACTORY_GEN_SKEL_SRCS := 
FACTORY_GEN_HDRS := $(FACTORY_GEN_IOR_HDRS) $(FACTORY_GEN_STUB_HDRS)
FACTORY_GEN_SRCS := $(FACTORY_GEN_IOR_SRCS) $(FACTORY_GEN_STUB_SRCS) $(FACTORY_GEN_SKEL_SRCS) 
FACTORY_OFILES := $(subst .cc,.o,$(FACTORY_GEN_SRCS) $(FACTORY_IMPL_SRCS))
FACTORY_OFILES := $(subst .c,.o, $(FACTORY_OFILES))

SERVER_IMPL_HDRS := $(filter-out iBase_Error_Impl.hh $(FACTORY_IMPL_HDRS),$(IMPLHDRS))
SERVER_IMPL_SRCS := $(filter-out iBase_Error_Impl.cc $(FACTORY_IMPL_SRCS),$(IMPLSRCS))
SERVER_GEN_IOR_HDRS := $(filter-out iBase_Error_IOR.h $(FACTORY_GEN_IOR_HDRS),$(IORHDRS))
SERVER_GEN_IOR_SRCS := $(filter-out iBase_Error_IOR.c $(FACTORY_GEN_IOR_SRCS),$(IORSRCS))
SERVER_GEN_STUB_HDRS := $(filter-out iBase_Error.hh $(FACTORY_GEN_STUB_HDRS),$(STUBHDRS))
SERVER_GEN_STUB_SRCS := $(filter-out iBase_Error.cc $(FACTORY_GEN_STUB_SRCS),$(STUBSRCS))
SERVER_GEN_SKEL_SRCS := $(filter-out iBase_Error_Skel.cc $(FACTORY_GEN_SKEL_SRCS),$(SKELSRCS))
SERVER_GEN_HDRS := $(SERVER_GEN_IOR_HDRS) $(SERVER_GEN_STUB_HDRS)
SERVER_GEN_SRCS := $(SERVER_GEN_IOR_SRCS) $(SERVER_GEN_STUB_SRCS) $(SERVER_GEN_SKEL_SRCS) 
SERVER_OFILES := $(subst .cc,.o,$(SERVER_GEN_SRCS) $(SERVER_IMPL_SRCS))
SERVER_OFILES := $(subst .c,.o, $(SERVER_OFILES))

SERVER_IMPL_OFILES = 


-include Cclient/babel.make
CCLIENT_GEN_HDRS := $(IORHDRS) $(STUBHDRS)
CCLIENT_GEN_SRCS := $(STUBSRCS)
CCLIENT_OFILES := $(subst .cc,.o,$(CCLIENT_GEN_SRCS))
CCLIENT_OFILES := $(subst .c,.o, $(CCLIENT_OFILES))

-include F77client/babel.make
F77CLIENT_GEN_HDRS := $(IORHDRS) $(STUBHDRS)
F77CLIENT_GEN_SRCS := $(STUBSRCS)
F77CLIENT_OFILES := $(subst .cc,.o,$(F77CLIENT_GEN_SRCS))
F77CLIENT_OFILES := $(subst .c,.o, $(F77CLIENT_OFILES))

all: testcxx
# testc_cbind

settings:
	@echo "FACTORY_OFILES = $(FACTORY_OFILES)"
	@echo "SERVER_OFILES  = $(SERVER_OFILES)"

#============= iBase =============
include bserver/babel.make
BSERVER_IMPL_HDRS := $(IMPLHDRS)
BSERVER_IMPL_SRCS := $(IMPLSRCS)
BSERVER_GEN_HDRS := $(IORHDRS) $(STUBHDRS)
BSERVER_GEN_SRCS := $(IORSRCS) $(SKELSRCS) $(STUBSRCS)
BSERVER_OFILES := $(subst .cc,.o,$(BSERVER_GEN_SRCS) $(BSERVER_IMPL_SRCS))
BSERVER_OFILES := $(subst .c,.o, $(BSERVER_OFILES))

brepo/.timestamp: iBase.sidl
	-rm -rf brepo
	$(BABEL) -tXML -o brepo iBase.sidl 
	touch $@ 

bserver/babel.make: brepo/.timestamp iBase.sidl
	-rm -f $@
	$(BABEL) -sC++ -o bserver iBase.sidl

${iBase_DIR}/libiBaseserver.so: bserver/babel.make ${BSERVER_OFILES}
	$(LD) $(LD_FLAGS) $(LD_FLAGS_PIC) -o $@ ${BSERVER_OFILES}

${iBase_DIR}/libiBaseserver.a: bserver/babel.make ${BSERVER_OFILES}
	$(ARCHIVER) $@ ${BSERVER_OFILES}

#============= iBase =============

mrepo/.timestamp: iMesh.sidl 
	-rm -rf mrepo
	$(BABEL) -tXML -o mrepo -Rbrepo iMesh.sidl
	touch $@ 

mserver/babel.make: mrepo/.timestamp brepo/.timestamp iMesh_SIDL.sidl iMesh.sidl 
	-rm -f $@
	$(BABEL) -R"mrepo;brepo" -sC++ -o mserver iMesh_SIDL.sidl

Cclient: mrepo/.timestamp brepo/.timestamp iMesh_SIDL.sidl
	-rm -rf Cclient
	$(BABEL) -R"mrepo;brepo" -cC -o Cclient iMesh_SIDL.sidl

F77client: mrepo/.timestamp brepo/.timestamp iMesh_SIDL.sidl
	-rm -rf F77client
	$(BABEL) -R"mrepo;brepo" -cF77 -o F77client iMesh_SIDL.sidl

${iMesh_DIR}/libiMesh.a: ${SERVER_IMPL_OFILES}
	$(ARCHIVER) $@ ${SERVER_IMPL_OFILES}

${iMesh_DIR}/libiMeshserver.so: mserver/babel.make ${SERVER_OFILES} ${iMesh_CBIND_SERVER_FILES} ${iBase_SERVER_FILES}
	$(LD) $(LD_FLAGS) $(LD_FLAGS_PIC) -o $@ ${SERVER_OFILES} ${iMesh_CBIND_SERVER_LIBS}

libiMeshserver.a: mserver/babel.make ${SERVER_OFILES} ${SERVER_IMPL_OFILES}
	$(ARCHIVER) $@ ${SERVER_OFILES} ${SERVER_IMPL_OFILES}

libiMeshCclient.so: Cclient ${CCLIENT_OFILES}
	$(LD) $(LD_FLAGS) $(LD_FLAGS_PIC) -o $@ ${CCLIENT_OFILES}

libiMeshF77client.so: F77client ${F77CLIENT_OFILES}
	$(LD) $(LD_FLAGS) $(LD_FLAGS_PIC) -o $@ ${F77CLIENT_OFILES}

Cclient/babel.make: Cclient

F77client/babel.make: F77client

testcxx: ${iMesh_SERVER_FILES} ${iBase_SERVER_FILES} testcxx.o 
	$(LD)  ${LD_FLAGS} -o $@ testcxx.o -Wl,-rpath,. -L. ${iMesh_SERVER_LIBS} ${iMesh_CBIND_SERVER_LIBS}
testc: libiMeshCclient.so ${iMesh_SERVER_FILES} ${iMesh_CLIENT_C_FILES} ${iBase_SERVER_FILES} testc.o 
	$(LD)  -o $@ testc.o -Wl,-rpath,. -L. ${iMesh_CLIENT_C_LIBS} ${iMesh_SERVER_LIBS}

testc_cbind.o : testc_cbind.c
	$(CC) -std=c99 -DIS_BUILDING_MB $(CC_FLAGS) -o $@ -c $<

testc_cbind:  ${iMesh_CBIND_SERVER_FILES} ${iBase_CBIND_SERVER_FILES} testc_cbind.o 
	$(LD)  -o $@ testc_cbind.o ${iMesh_CBIND_SERVER_LIBS}

clean:
	-rm -f *.o *.so testc testcxx *~ core.*

clean_all: clean
	-rm -rf F77client Cclient mrepo
	cd mserver; rm -f ${SERVER_GEN_HDRS} ${SERVER_GEN_SRCS}

depend:
	${MKDEP} -DIS_BUILDING_MB ${INCPATH} testc.c testcxx.cpp > make.dependencies


%.o : %.cxx
	$(CXX) -DIS_BUILDING_MB $(CXX_FLAGS) -o $@ -c $<

%.o : %.cpp
	$(CXX) -DIS_BUILDING_MB $(CXX_FLAGS) -o $@ -c $<

%.o : %.cc
	$(CXX) -DIS_BUILDING_MB $(CXX_FLAGS) -o $@ -c $<

%.o : %.c
	$(CC) -DIS_BUILDING_MB $(CC_FLAGS) -o $@ -c $<

%.o : %.f
	$(FC) -DIS_BUILDING_MB $(FC_FLAGS) -o $@ -c $<

#include make.dependencies
