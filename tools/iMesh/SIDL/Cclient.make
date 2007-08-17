include babel.make

TARGET_LIB = libiMeshCclient.la
LTLIBS = ../mserver/libiMeshserver.la
C_SRCS = $(STUBSRCS)
CXX_SRCS =
LIBLINK = $(LINK)
bdir = Cclient

include ../common.make
