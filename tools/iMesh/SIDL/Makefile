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
CXX = g++-4.2
CXX_FLAGS_REQ = -fpic
CXX_FLAGS_DBG = -g -Wall
CXX_FLAGS_OPT = -O 

CC  = gcc-4.2
CC_FLAGS_REQ = -fpic
CC_FLAGS_DBG = -g -Wall
CC_FLAGS_OPT = -O

FC  = gfortran-4.2
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


INCPATH = -I. ${iMesh_SIDL_SERVER_INCLUDES}


CXX_FLAGS = $(CXX_FLAGS_REQ) $(CXX_FLAGS_DBG) $(INCPATH) ${iMesh_SIDL_CCLIENT_INCLUDES} 
CC_FLAGS = $(CC_FLAGS_REQ) $(CC_FLAGS_DBG) $(INCPATH) ${iMesh_SIDL_CCLIENT_INCLUDES} 
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

repo/.btimestamp: iBase.sidl
	-rm -rf repo
	$(BABEL) -tXML -o repo iBase.sidl 
	touch $@ 

bserver/babel.make: repo/.btimestamp iBase.sidl
	-rm -f $@
	$(BABEL) -sC++ -o bserver iBase.sidl

${iBase_DIR}/libiBaseserver.so: bserver/babel.make ${BSERVER_OFILES}
	$(LD) $(LD_FLAGS) $(LD_FLAGS_PIC) -o $@ ${BSERVER_OFILES}

${iBase_DIR}/libiBaseserver.a: bserver/babel.make ${BSERVER_OFILES}
	$(ARCHIVER) $@ ${BSERVER_OFILES}

#============= iBase =============

repo/.mtimestamp: iMesh.sidl 
	$(BABEL) -tXML -o repo -Rrepo iMesh.sidl
	touch $@ 

mserver/babel.make: repo/.mtimestamp repo/.btimestamp iMesh_SIDL.sidl iMesh.sidl 
	-rm -f $@
	$(BABEL) -R"repo" -sC++ -o mserver iMesh_SIDL.sidl

Cclient: repo/.mtimestamp repo/.mtimestamp iMesh_SIDL.sidl
	-rm -rf Cclient
	$(BABEL) -R"repo" -cC -o Cclient iMesh_SIDL.sidl

F77client: repo/.mtimestamp repo/.btimestamp iMesh_SIDL.sidl
	-rm -rf F77client
	$(BABEL) -R"repo" -cF77 -o F77client iMesh_SIDL.sidl

${iMesh_DIR}/libiMesh.a: ${SERVER_IMPL_OFILES}
	$(ARCHIVER) $@ ${SERVER_IMPL_OFILES}

${iMesh_DIR}/libiMeshserver.so: mserver/babel.make ${SERVER_OFILES} ${iBase_SIDL_SERVER_FILES}
	$(LD) $(LD_FLAGS) $(LD_FLAGS_PIC) -o $@ ${SERVER_OFILES} ${iMesh_SIDL_SERVER_LIBS}

libiMeshCclient.so: Cclient ${CCLIENT_OFILES}
	$(LD) $(LD_FLAGS) $(LD_FLAGS_PIC) -o $@ ${CCLIENT_OFILES}

libiMeshF77client.so: F77client ${F77CLIENT_OFILES}
	$(LD) $(LD_FLAGS) $(LD_FLAGS_PIC) -o $@ ${F77CLIENT_OFILES}

Cclient/babel.make: Cclient

F77client/babel.make: F77client

testcxx: ${iMesh_SIDL_SERVER_FILES} ${iBase_SIDL_SERVER_FILES} testcxx.o 
	$(LD)  ${LD_FLAGS} -o $@ testcxx.o -Wl,-rpath,. -L. ${iMesh_SIDL_SERVER_LIBS}
testc: libiMeshCclient.so ${iMesh_SIDL_SERVER_FILES} ${iMesh_SIDL_CCLIENT_FILES} testc.o 
	$(LD)  -o $@ testc.o -Wl,-rpath,. -L. ${iMesh_SIDL_CCLIENT_LIBS} \
              ${iMesh_SIDL_SERVER_LIBS} 

clean:
	-rm -f *.o *.so testc testcxx *~ core.*

clean_all: clean
	-rm -rf F77client Cclient repo
	cd mserver; rm -f ${SERVER_GEN_HDRS} ${SERVER_GEN_SRCS}

depend:
	${MKDEP} -DIS_BUILDING_MB ${INCPATH} testc.c testcxx.cpp > make.dependencies


%.o : %.cxx
	$(CXX) $(CXX_FLAGS) -o $@ -c $<

%.o : %.cpp
	$(CXX) $(CXX_FLAGS) -o $@ -c $<

%.o : %.cc
	$(CXX) $(CXX_FLAGS) -o $@ -c $<

%.o : %.c
	$(CC) $(CC_FLAGS) -o $@ -c $<

%.o : %.f
	$(FC) $(FC_FLAGS) -o $@ -c $<

#include make.dependencies
