# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ${MOAB_DIR}/lib/moab.make
include ${MOAB_DIR}/lib/iMesh-Defs.inc

.SUFFIXES: .o .cpp .F90

# MESH_DIR is the top-level MOAB source directory, used to locate mesh files that come with MOAB source
MESH_DIR="../MeshFiles/unittest"

StructuredMeshSimple : StructuredMeshSimple.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK} 

HelloMOAB : HelloMOAB.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

ReduceExchangeTags : ReduceExchangeTags.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

HelloMoabPar: HelloMoabPar.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK} 

TestExodusII: TestExodusII.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

DirectAccessWithHoles: DirectAccessWithHoles.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

DirectAccessNoHoles: DirectAccessNoHoles.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

DirectAccessNoHolesF90: DirectAccessNoHolesF90.o
	${MOAB_CXX} -o $@ $< ${IMESH_LIBS}

.cpp.o :
	${MOAB_CXX} ${MOAB_CXXFLAGS} ${MOAB_INCLUDES} -DMESH_DIR=\"${MESH_DIR}\" -c $<

.F90.o :
	${IMESH_FC} ${IMESH_FCFLAGS} ${IMESH_INCLUDES} ${IMESH_FCDEFS} -DMESH_DIR=\"${MESH_DIR}\" -c $<

