include babel.make

TARGET_LIB = libiMeshCclient.la
C_SRCS = $(STUBSRCS)
CXX_SRCS =
LIBLINK = $(LINK)
bdir = Cclient

include ../common.make
