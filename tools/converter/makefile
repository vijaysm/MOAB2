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

MOAB_DIR = ../..
include $(MOAB_DIR)/MB.common

convert: convert.cpp 
	$(CXX) $(DEBUG_FLAG) $(MACH_CXXFLAGS) \
	-I$(MOAB_DIR) -L$(MOAB_DIR) -lMOAB \
	convert.cpp -o convert $(MOAB_APP_LIBS)

../../libMOAB.so:
	cd ../.. && make libMOAB.so

clean:
	rm -f convert
