# NOTE:  The 'makefile' in the *installed* examples directory is generated
#        using the install-data-hook target below.  This logic assumes that:
#        1) All examples are listed in check_PROGRAMS
#        2) All examples have a single source file with the same name (case
#           sensitive) as the example executable with a .cpp suffix.

check_PROGRAMS = FileRead \
                 GeomSetHierarchy \
                 GetEntities \
                 KDTree \
		 ObbTree \
                 SetsNTags \
                 SkinMesh \
                 SurfArea

LDADD = ../src/libMOAB.la
AM_CPPFLAGS += -I../src \
               -I$(srcdir)/../src \
               -I$(srcdir)/../src/parallel \
               -I$(srcdir)/../src/oldinc

FileRead_SOURCES = FileRead.cpp
GeomSetHierarchy_SOURCES = GeomSetHierarchy.cpp
GetEntities_SOURCES = simple/GetEntities.cpp
SetsNTags_SOURCES = SetsNTags.cpp
SkinMesh_SOURCES = SkinMesh.cpp
SurfArea_SOURCES = SurfArea.cpp
KDTree_SOURCES = KDTree.cpp
ObbTree_SOURCES = ObbTree.cpp

exampledir = ${docdir}/examples
nobase_example_DATA = \
               examples.make \
               FileRead.cpp \
               GeomSetHierarchy.cpp \
               simple/GetEntities.cpp \
               simple/makefile \
               KDTree.cpp \
	       ObbTree.cpp \
               SetsNTags.cpp \
               SkinMesh.cpp \
               SurfArea.cpp 
               
if ENABLE_imesh
  imesh_DIR = itaps
else
  imesh_DIR = 
endif 

SUBDIRS = $(imesh_DIR)

ex_make = $(DESTDIR)$(exampledir)/makefile
install-data-hook:
	$(AM_V_at)rm -f $(ex_make)
	$(AM_V_at)mv $(DESTDIR)$(exampledir)/examples.make $(ex_make)
	$(AM_V_at)echo "all: $(check_PROGRAMS)" >>$(ex_make)
	$(AM_V_at)rule='	$$(CXX) -o $$@ $$< $$(CXXFLAGS) $$(MOAB_INCLUDES) $$(MOAB_LIBS_LINK)'; \
	for example in $(check_PROGRAMS); do \
	  echo >>$(ex_make); \
	  echo "$${example}: $${example}.cpp" >>$(ex_make); \
	  echo "$$rule" >>$(ex_make); \
	done

uninstall-hook:
	$(AM_V_at)rm -f $(DESTDIR)$(exampledir)/makefile
