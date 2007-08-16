include babel.make

TARGET_LIB = libiMeshFclient.la
C_SRCS = $(STUBSRCS)
CXX_SRCS =
LIBLINK = $(LINK)
bdir = Fclient

include ../common.make
