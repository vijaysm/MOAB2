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

.phony: default all clean force_rebuild

default: lib/libmhdf-pic.a

all: lib/libmhdf.a lib/libmhdf.so lib/libmhdf-pic.a

lib/libmhdf.a: force_rebuild lib
	make -C obj ../$@

lib/libmhdf-pic.a: force_rebuild lib
	make -C obj-pic ../$@

lib/libmhdf.so: force_rebuild lib
	make -C obj-pic ../$@

lib:
	mkdir lib

clean:
	make -C obj clean 
	make -C obj-pic clean

force_rebuild:
