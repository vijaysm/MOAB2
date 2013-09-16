# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ${MOAB_DIR}/lib/moab.make
include ${MOAB_DIR}/lib/iMesh-Defs.inc

.SUFFIXES: .o .cpp .F90

# MESH_DIR is the directory containing mesh files that come with MOAB source
MESH_DIR="../MeshFiles/unittest"

EXAMPLES = HelloMOAB GetEntities SetsNTags StructuredMeshSimple DirectAccessWithHoles DirectAccessNoHoles 
PAREXAMPLES = HelloParMOAB ReduceExchangeTags LloydRelaxation
F90EXAMPLES = DirectAccessNoHolesF90
EXOIIEXAMPLES = TestExodusII

default: ${EXAMPLES}

HelloMOAB : HelloMOAB.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

GetEntities: GetEntities.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

SetsNTags: SetsNTags.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

LloydRelaxation: LloydRelaxation.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

StructuredMeshSimple : StructuredMeshSimple.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK} 

DirectAccessWithHoles: DirectAccessWithHoles.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

DirectAccessNoHoles: DirectAccessNoHoles.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

DirectAccessNoHolesF90: DirectAccessNoHolesF90.o
	${MOAB_CXX} -o $@ $< ${IMESH_LIBS}

ReduceExchangeTags : ReduceExchangeTags.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

HelloParMOAB: HelloParMOAB.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK} 

TestExodusII: TestExodusII.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

DeformMeshRemap: DeformMeshRemap.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

clean:
	rm -rf *.o ${EXAMPLES} ${PAREXAMPLES} ${EXOIIEXAMPLES}

.cpp.o :
	${MOAB_CXX} ${MOAB_CXXFLAGS} ${MOAB_INCLUDES} -DMESH_DIR=\"${MESH_DIR}\" -c $<

.F90.o :
	${IMESH_FC} ${IMESH_FCFLAGS} ${IMESH_INCLUDES} ${IMESH_FCDEFS} -DMESH_DIR=\"${MESH_DIR}\" -c $<

