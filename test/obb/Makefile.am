check_PROGRAMS = obb_test obb_time obb_tree_tool
TESTS = obb_test

@CGM_CONFIG_OPTIONS@

DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_builddir)
LDADD = $(top_builddir)/libMOAB.la

if HAVE_CGM
  LDADD += $(CGM_LTFLAGS) 
endif

obb_test_SOURCES = obb_test.cpp
obb_time_SOURCES = obb_time.cpp
obb_tree_tool_SOURCES = obb_tree_tool.cpp
