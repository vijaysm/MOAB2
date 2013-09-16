# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ${MOAB_DIR}/lib/moab.make
include ${MOAB_DIR}/lib/iMesh-Defs.inc

.SUFFIXES: .o .cpp .F90

# MESH_DIR is the directory containing mesh files that come with MOAB source
MESH_DIR="../MeshFiles/unittest"

EXAMPLES = HelloMOAB GetEntities SetsNTags StructuredMeshSimple DirectAccessWithHoles DirectAccessNoHoles 
PAREXAMPLES = HelloParMOAB ReduceExchangeTags LloydRelaxation
F90EXAMPLES = DirectAccessNoHolesF90 PushParMeshIntoMoabF90
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

point_in_elem_search: point_in_elem_search.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

PushParMeshIntoMoabF90: PushParMeshIntoMoabF90.o
	${MOAB_CXX} -o $@ $< ${IMESH_LIBS} -lgfortran -L/usr/lib/openmpi/lib -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -ldl 

DeformMeshRemap: DeformMeshRemap.o
	${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

clean:
	rm -rf *.o ${EXAMPLES} ${PAREXAMPLES} ${EXOIIEXAMPLES}

.cpp.o :
	${MOAB_CXX} ${CXXFLAGS} ${MOAB_CXXFLAGS} ${MOAB_CPPFLAGS} ${MOAB_INCLUDES} -DMESH_DIR=\"${MESH_DIR}\" -c $<

 .F90.o :
	${IMESH_FC} ${FCFLAGS} ${IMESH_FCFLAGS} ${MOAB_CPPFLAGS} ${IMESH_INCLUDES} ${IMESH_FCDEFS} -DMESH_DIR=\"${MESH_DIR}\" -c $<

