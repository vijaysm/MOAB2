include ${MOAB_DIR}/lib/moab.make

# MESH_DIR used to locate mesh files that come with MOAB source, 
# i.e. ${MOAB_DIR}/MeshFiles/unittest for in-source builds
MESH_DIR=../MeshFiles/unittest

structuredmesh : structuredmesh.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK} 

.cpp.o:
	${MOAB_CXX} -c ${MOAB_CXXFLAGS} ${MOAB_INCLUDES} $<
