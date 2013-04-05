# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ${MOAB_DIR}/lib/moab.make

# MESH_DIR is the top-level MOAB source directory, used to locate mesh files that come with MOAB source
MESH_DIR="../MeshFiles/unittest"

StructuredMeshSimple : StructuredMeshSimple.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK} 

.cpp.o:
	${MOAB_CXX} -c ${MOAB_CXXFLAGS} ${MOAB_INCLUDES} -DMESH_DIR=\"${MESH_DIR}\" -E $<
