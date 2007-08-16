include babel.make

TARGET_LIB = libiMeshserver.la
C_SRCS = $(IORSRCS)
CXX_SRCS = $(IMPLSRCS) $(SKELSRCS) $(STUBSRCS)
LIBLINK = $(CXXLINK)
bdir = mserver

include ../common.make
