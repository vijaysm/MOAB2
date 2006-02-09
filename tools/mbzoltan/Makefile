##############################################################################
# Zoltan Library for Parallel Applications                                   #
# Copyright (c) 2000,2001,2002, Sandia National Laboratories.                #
# For more info, see the README file in the top-level Zoltan directory.      # 
##############################################################################
#	Author:  Vitus Leung

##############################################################################
# Makefile for simple example which uses the Zoltan and MOAB libraries
#
# To compile the example, type
#    gmake ZOLTAN_ARCH=<targetmachine> 
# where
#    targetmachine = the suffix of one of the configuration files Config.*
#                    in Zoltan/Utilities/Config.
##############################################################################

#
# relative to this directory
#
ZOLTAN_TOP = ../../../Zoltan
ZOLTAN_UTIL = ../../../Zoltan/Utilities
MOAB_TOP = ../..

#
# relative to the subdirectories
#
ZTOP = ../$(ZOLTAN_TOP)
ZUTIL = ../$(ZOLTAN_UTIL)
MTOP = ../$(MOAB_TOP)

##############################################################################
# Process architecture flag
##############################################################################

ifndef ZOLTAN_ARCH
ZOLTAN_ARCH	= generic
endif

ifndef ZOLTAN_CONFIG
ZOLTAN_CONFIG = $(ZOLTAN_TOP)/Utilities/Config/Config.$(ZOLTAN_ARCH)
endif

include $(ZOLTAN_CONFIG)

ifndef ZOLTAN_OBJ_DIR
ZOLTAN_OBJ_DIR = $(ZTOP)/Obj_$(ZOLTAN_ARCH)
endif

ifndef ZOLTAN_LIB_DIR
ZOLTAN_LIB_DIR = $(ZOLTAN_OBJ_DIR)
endif

############################################################################
#                 SUB DIRECTORIES
############################################################################
#

SUBDIRS = lib mbzoltan

############################################################################
#                 LIBRARIES
############################################################################
#
LIBDIR = -L../lib -L$(ZOLTAN_LIB_DIR) $(MPI_LIBPATH)

THIRD_PARTY_LIBS=

ifdef PARMETIS_LIBPATH
  THIRD_PARTY_LIBS += -lparmetis -lmetis
  LIBDIR += $(PARMETIS_LIBPATH)
endif

ifdef JOSTLE_LIBPATH
  THIRD_PARTY_LIBS += -ljostle
  LIBDIR += $(JOSTLE_LIBPATH)
endif

ifdef PATOH_LIBPATH
  THIRD_PARTY_LIBS += -lpatoh
  LIBDIR += $(PATOH_LIBPATH)
endif


LIBS =  -lmbzoltan -lzoltan $(THIRD_PARTY_LIBS) $(MPI_LIB) -lm\
	-lMOAB -lhdf5 -lz -lnetcdf_c++ -lnetcdf

ZOLTAN_INCLUDE= -I$(ZUTIL)/Communication \
                -I$(ZUTIL)/DDirectory \
                -I$(ZUTIL)/Memory \
                -I$(ZUTIL)/Timer \
                -I$(ZUTIL)/shared \
                -I$(ZTOP)/zz \
                -I$(ZTOP)/par \
                -I$(ZTOP)/lb \
                -I$(ZTOP)/include

MOAB_INCLUDE= -I$(MTOP)

INCLUDE = -I../lib $(ZOLTAN_INCLUDE) $(MOAB_INCLUDE) $(INCLUDE_PATH)

##############################################################################
# 
##############################################################################

TARGETS= mbzoltanLibrary mbzoltan

all: $(TARGETS)

mbzoltanLibrary:
	@cd lib; $(MAKE) all "CC=$(CC)" "AR=$(AR)" "CFLAGS=$(CFLAGS)"\
                            "INCLUDE=$(INCLUDE)" 

mbzoltan: mbzoltanLibrary
	@cd mbzoltan; $(MAKE) all "CC=$(CPPC)" "CFLAGS=$(CFLAGS)" \
			"INCLUDE=$(INCLUDE)" \
			  "LIBDIR=$(LIBDIR)" "LIBS=$(LIBS)"


clean:
	@cd lib; $(MAKE) clean
	@cd mbzoltan; $(MAKE) clean

