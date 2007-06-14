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
include Itaps-Defs.inc

iBase_CBIND_SERVER_INCLUDES = \
      -I${iBase_DIR}

iBase_CBIND_SERVER_LIBS = 

iBase_CBIND_SERVER_FILES = 

iBase_SERVER_INCLUDES = \
      -I${iBase_DIR} -I${iBase_DIR}/bserver -I${SIDLDIR}/include

iBase_SERVER_LIBS = \
      -L${iBase_DIR} -liBaseserver -L$(SIDLDIR)/lib -lsidl \
      -Wl,-rpath,${iBase_DIR}

iBase_SERVER_FILES = \
      ${iBase_DIR}/libiBaseserver.so

