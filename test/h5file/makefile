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

MOAB_DIR=../..
include $(MOAB_DIR)/MB.common

dump: h5test.h5m
	- h5dump $< > $@

h5test.h5m: h5file_test ../../libMOAB.so
	- h5file_test -k

h5file_test: h5file_test.cpp
	g++ $(CXXFLAGS) $(DEBUG_FLAG) -Wl,-rpath,$(MOAB_DIR) \
	-o $@ -I $(MOAB_DIR) $< $(PCOM_APP_LIBS)

