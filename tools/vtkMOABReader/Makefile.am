AM_CPPFLAGS += -I$(top_srcdir)/src -I$(top_builddir)/src $(VTK_INCLUDES)
lib_LTLIBRARIES = libvtkMOABReader.la
libvtkMOABReader_la_SOURCES = vtkMOABReader.cxx
libvtkMOABReader_la_include_HEADERS = vtkMOABReader.h
libvtkMOABReader_la_includedir = $(includedir)

LDADD = libvtkMOABReader.la $(VTK_LIBS) $(top_builddir)/src/libMOAB.la 
AM_LDFLAGS = -R $(VTK_LIB_DIR)

