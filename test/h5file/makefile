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

MOAB_DIR=/home/jason/moab/rw
include $(MOAB_DIR)/MB.common
NP ?= 1

ifdef HDF5_DIR
  HDF5_RLINK=-Wl,-rpath,$(HDF5_DIR)/lib
endif

  # Run serial test
default: h5test.h5m

  # Run parallel test, too
  # NOTE: You probably want to specify the number of processors
  #       on the command line. e.g.:  make NP=4
  #       The default is 1, which is safe but uninteresting.
ifeq ($(PARALLEL),yes)
  default: ptest
endif

force:

h5test.h5m: h5file_test ../../libMOAB.so
	- h5file_test -k

h5file_test.o: h5file_test.cpp
	g++ -c $(CXXFLAGS) $(DEBUG_FLAG) -o $@ -I $(MOAB_DIR) $< 

h5file_test: h5file_test.o
	g++ $(CXXFLAGS) -Wl,-rpath,$(MOAB_DIR) -o $@ $< $(PCOM_APP_LIBS) $(HDF5_RLINK)

parallel.o: parallel.cpp 
	$(MPICXX) -c $(CXXFLAGS) $(DEBUG_FLAG) -I $(MOAB_DIR) \
	-I$(MOAB_DIR)/mhdf/include $(HDF5_INC) -o $@ $< 

parallel: parallel.o
	$(MPICXX) $(CXXFLAGS) $(DEBUG_FLAG) \
	-Wl,-rpath,$(MOAB_DIR) \
	$(PCOM_APP_LIBS) \
	$(HDF5_RLINK) \
	-o $@ $< 


ptest: parallel ptest.cub force
	- $(MPIRUN) $(DBG) -np $(NP) $< $@ ptest.cub
	- x=-1;                               \
	  while [ $$((++x)) -lt $(NP) ] ; do  \
	    num=`printf "%02d" $$x`;          \
	    file="$@.$${num}.h5m";            \
	    test -f $$file      && h5dump $$file        > dump.$$num; \
	    done
	- test -f $@            && h5dump $@            > dump.parallel
	- test -f $@.serial.h5m && h5dump $@.serial.h5m > dump.serial
	- test -f $@            && h5dump -A $@            > summary.parallel
	- test -f $@.serial.h5m && h5dump -A $@.serial.h5m > summary.serial

dump: h5test.h5m 
	- h5dump $< > $@

clean:
	rm -f dump gmon.out h5file_test h5file_test.o h5test.h5m parallel.o
