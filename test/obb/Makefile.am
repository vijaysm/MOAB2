check_PROGRAMS = obb_test obb_time obb_tree_tool
TESTS = obb_test

AM_CPPFLAGS += -DSRCDIR=$(srcdir) \
               -I$(top_srcdir)/src -I$(top_builddir)/src -I.. \
               -I$(top_srcdir)/src/moab -I$(top_builddir)/src/moab 

LDADD = $(top_builddir)/src/libMOAB.la

obb_test_SOURCES = obb_test.cpp
obb_time_SOURCES = obb_time.cpp
obb_tree_tool_SOURCES = obb_tree_tool.cpp
