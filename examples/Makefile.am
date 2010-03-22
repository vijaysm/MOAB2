check_PROGRAMS = FileRead \
                 GeomSetHierarchy \
                 GetEntities \
                 KDTree \
                 SetsNTags \
                 SkinMesh \
                 SurfArea

LDADD = ../src/libMOAB.la
AM_CPPFLAGS += -I../src \
               -I$(srcdir)/../src \
               -I$(srcdir)/../src/oldinc

FileRead_SOURCES = FileRead.cpp
GeomSetHierarchy_SOURCES = GeomSetHierarchy.cpp
GetEntities_SOURCES = simple/GetEntities.cpp
SetsNTags_SOURCES = SetsNTags.cpp
SkinMesh_SOURCES = SkinMesh.cpp
SurfArea_SOURCES = SurfArea.cpp
KDTree_SOURCES = KDTree.cpp

exampledir = ${docdir}/examples
nobase_example_DATA = \
               examples.make \
               FileRead.cpp \
               GeomSetHierarchy.cpp \
               simple/GetEntities.cpp \
               simple/makefile \
               KDTree.cpp \
               SetsNTags.cpp \
               SkinMesh.cpp \
               SurfArea.cpp 
               

install-data-hook:
	mv -f $(DESTDIR)$(exampledir)/examples.make $(DESTDIR)$(exampledir)/makefile

uninstall-hook:
	rm -f $(DESTDIR)$(exampledir)/makefile
