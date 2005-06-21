# 
# Copyright (c) 2005 Lawrence Livermore National Laboratory under
# contract number B545069 with the University of Wisconsin - Madison.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
# NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

MOAB_DIR = ../..
include $(MOAB_DIR)/MB.common

MOAB_IFACE_HEADERS = MBCore.hpp \
                     MBRange.hpp \
                     MBInterface.hpp

DEPENDS = $(MOAB_IFACE_HEADERS:%=../../%)

propagate_tags: propagate_tags.o parse.o
	$(CXX) $(DEBUG_FLAG) $(MACH_CXXFLAGS) \
	$^ -o $@ \
	-L$(MOAB_DIR) -lMOAB $(MOAB_APP_LIBS) \
	-Wl,-rpath,`pwd`/../..
	

propagate_tags.o: parse.hpp $(DEPENDS)

parse.o: parse.hpp $(DEPENDS)

.cpp.o:
	$(CXX) $(DEBUG_FLAG) $(MACH_CXXFLAGS) -I$(MOAB_DIR) \
	-c $< -o $@  

../../libMOAB.so:
	cd ../.. && make libMOAB.so

clean:
	rm -f convert
