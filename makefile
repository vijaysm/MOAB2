
# this makefile just contains custom build rules, and controls
# what gets built


# there are 3 main options for building the MB library
#  1. static libMOAB.a
#  2. shared libMOABxpcom.so  (uses xpcom framework)
#  3. shared libMOAB.so   (uses pcom framework)
#
#  to build the pcom files, just type "make pcom"
#  to build the xpcom files, just type "make xpcom"
#  
#  additional instructions for xpcom build:
#    1. get xpcom source and build
#       a. setenv CVSROOT :pserver:anonymous@cvs-mirror.mozilla.org:/cvsroot
#       b. cvs co mozilla/client.mk
#       c. cd mozilla
#       d. gmake -f client.mk pull_all BUILD_MODULES=xpcom
#       e. configure --enable-modules=xpcom
#       f. gmake
#    2. setenv MOZILLA_DEV_DIR path_to_mozilla_dir
#    3. add ${MOZILLA_DEV_DIR}/dist/lib to your LD_LIBRARY_PATH
#         !!!!! warning: this could disable your mozilla browser since the libs could conflict !!!!!
#         !!!!! warning: if you have another version of libxpcom.so, you'll get conflicts      !!!!!
#

MOAB_IMPL_VERSION = 0.99

include MB.$(MACHINE_TYPE)

SHELL = /bin/sh

include MB.local

# the source files for building the MB library
MB_LIB_SRCS = \
	AEntityFactory.cpp \
	DenseTagCollections.cpp \
	DualTool.cpp \
	EntitySequence.cpp \
	EntitySequenceManager.cpp \
	ExoIIUtil.cpp \
	GeomTopoTool.cpp \
	HigherOrderFactory.cpp \
	HomXform.cpp \
	MBBits.cpp \
	MBCN.cpp \
	MBCore.cpp \
	MBFactory.cpp \
	MBMeshSet.cpp \
	MBRange.cpp \
	MBReadUtil.cpp \
	MBSkinner.cpp \
	MBUtil.cpp \
	MBWriteUtil.cpp \
	MeshTopoUtil.cpp \
	PolyEntitySequence.cpp \
	ReadNCDF.cpp \
	ReadVtk.cpp \
	ScdElementSeq.cpp \
	ScdVertexSeq.cpp \
	SparseTagCollections.cpp \
	TagServer.cpp \
	Tqdcfr.cpp \
	WriteNCDF.cpp

UTEST_SRCS = merge_test.cpp scdseq_test.cpp test_adj.cpp test_exo.cpp test_rms.cpp

DUM_LIB_HDRS = ${MB_LIB_SRCS:.cpp=.hpp} 
MB_LIB_HDRS = ${DUM_LIB_HDRS:MBFactory.hpp=} 

# object files for building the static library
MB_LIB_OBJS = ${MB_LIB_SRCS:.cpp=.o} 
# object files for building the xpcom moab library
XPCOM_MB_LIB_OBJS = ${MB_LIB_SRCS:.cpp=.xpcom.o}

# moab library targets and link flags
MB_LIB_TARGET = libMOAB.a

XPCOM_MB_LIB_TARGET = libMOABxpcom.so
XPCOM_LIB_LFLAGS = -shared 

PCOM_MB_LIB_TARGET = libMOAB.so

# additional xpcom flags  -- xpcom requires MOZILLA_DEV_DIR set
XPCOM_INCLUDE = \
    -I${MOZILLA_DEV_DIR}/dist/include/xpcom \
    -I${MOZILLA_DEV_DIR}/dist/include/nspr \
    -I${MOZILLA_DEV_DIR}/dist/include/string

XPCOM_LIBS = -L${MOZILLA_DEV_DIR}/dist/lib -lxpcom -lnspr4
XPCOM_MACH_CXXFLAGS = -DXPCOM_MB ${MACH_CXXFLAGS} ${XPCOM_INCLUDE}

# by default, lets just build moab_test using a static moab library
default : pcom

static : libMOAB.a moab_test.static homxform_test.static scdseq_test.static

# build everything
all : moab_test test_exo test_tag_server test_ent_seq xpcom pcom

# build the static moab library
libMOAB.a : ${MB_LIB_OBJS}
	@ echo Linking ...  $@
	${PREFIX} ${ARCHIVER} $@ ${MB_LIB_OBJS}

# build the xpcom moab library
xpcom : ${XPCOM_MB_LIB_TARGET} moab_test.xpcom

${XPCOM_MB_LIB_TARGET} : components ${XPCOM_MB_LIB_OBJS}
	@ echo Linking...  libMOABxpcom.so
	${PREFIX} ${LINKER} ${XPCOM_LIB_LFLAGS} ${XPCOM_MB_LIB_OBJS} \
	-o ${XPCOM_MB_LIB_TARGET} ${NETCDF_LIBS}

# build the pcom moab library
pcom : moab_test ${PCOM_MB_LIB_TARGET}

${PCOM_MB_LIB_TARGET} : components ${MB_LIB_OBJS}
	@ echo Linking...  ${PCOM_MB_LIB_TARGET}
	${PREFIX} ${LINKER} ${PCOM_LIB_LFLAGS} ${MB_LIB_OBJS} \
	-o ${PCOM_MB_LIB_TARGET} ${NETCDF_LIBS} ${PCOM_LIB_LIBS} -lm


# make the components directory if it doesn't exist
# this is where the xpcom and pcom libraries go
components : 
	@ echo Making "components" directory...
	@ mkdir components

runtest: moab_test homxform_test.static scdseq_test.static
	./moab_test
	./homxform_test.static
	./scdseq_test.static

dist: clean_all ${MB_LIB_SRCS} ${MB_LIB_HDRS} test doc/MOAB-UG.doc doxygen/moab MBUnknownInterface.h MB.* makefile README ${UTEST_SRCS}
	tar czf MOAB-${MOAB_IMPL_VERSION}.tar.gz ${MB_LIB_SRCS} ${MB_LIB_HDRS} test doc/MOAB-UG.doc doxygen/moab MBUnknownInterface.h MB.* makefile ${UTEST_SRCS} TSTT

# build moab_test which uses static MB library
moab_test.static : MBTest.static.o ${MB_LIB_TARGET}  
	@ echo Linking...  moab_test.static
	${PREFIX} ${LINKER} ${DEBUG_FLAG} ${MACH_LFLAGS} MBTest.static.o -L./ -lMOAB \
	${NETCDF_LIBS} ${PCOM_APP_LIBS} -lm \
	-o moab_test.static

# build moab_test which uses xpcom MB library
moab_test.xpcom : MBTest.xpcom.o MBDefines.o
	@ echo Linking...  moab_test.xpcom
	${PREFIX} ${LINKER} ${DEBUG_FLAG} ${MACH_LFLAGS} ${XPCOM_LIBS} $? -o $@

# build moab_test which uses pcom MB library
moab_test : MBTest.o MBRange.o MBCN.o HomXform.o ${PCOM_MB_LIB_TARGET}
	@ echo Linking...  moab_test
	${PREFIX} ${LINKER} ${DEBUG_FLAG} ${MACH_LFLAGS} MBRange.o MBTest.o MBCN.o HomXform.o -o $@ ${PCOM_APP_LIBS} 

homxform_test.static : HomXform.o
	@ echo Linking...  homxform_test.static
	${PREFIX} ${LINKER} ${DEBUG_FLAG} ${MACH_LFLAGS} -DTEST HomXform.cpp -o homxform_test.static

scdseq_test.static : scdseq_test.o ${MB_LIB_TARGET}
	@ echo Linking...  scdseq_test.static
	${PREFIX} ${LINKER} ${DEBUG_FLAG} ${MACH_LFLAGS} scdseq_test.o -L./ -lMOAB \
	${NETCDF_LIBS} ${PCOM_APP_LIBS} -lm -o scdseq_test.static

scdseq_timing : scdseq_timing.o ${MB_LIB_TARGET}
	@ echo Linking...  scdseq_timing
	${PREFIX} ${LINKER} ${DEBUG_FLAG} ${MACH_LFLAGS} scdseq_timing.o -L./ -lMOAB \
	${NETCDF_LIBS} ${PCOM_APP_LIBS} -lm -o scdseq_timing

tqdcfr: Tqdcfr.o MBCN.o ${MB_LIB_TARGET}
	${PREFIX} ${LINKER} ${DEBUG_FLAG} ${MACH_LFLAGS} -DIS_BUILDING_MB -DTEST_TQDCFR Tqdcfr.cpp MBCN.o -L./ -lMOAB \
	${NETCDF_LIBS} -o tqdcfr
#
# build component tests
test_rms : libMOAB.a test_rms.o
	${LINKER} ${DEBUG_FLAG} test_rms.o -L. -lMOAB -o test_rms

merge_test : libMOAB.a merge_test.o
	${LINKER} ${DEBUG_FLAG} merge_test.o -L. -lMOAB ${NETCDF_LIBS} -o merge_test

delete_test : libMOAB.a delete_test.o 
	${LINKER} ${DEBUG_FLAG} delete_test.o -L. -lMOAB ${NETCDF_LIBS} -o delete_test

test_exo : libMOAB.a test_exo.o
	${LINKER} ${DEBUG_FLAG} test_exo.o -L. -lMOAB ${NETCDF_LIBS} -o test_exo

test_tag_server : libMOAB.a TagServer.o_test
	${CXX} ${DEBUG_FLAG} -o test_tag_server TagServer.o_test -L. -lMOAB

test_ent_seq : libMOAB.a MBEntitySequence.o_test
	${CXX} ${DEBUG_FLAG} -o test_ent_seq MBEntitySequence.o_test -L. -lMOAB ${NETCDF_INCLUDE} ${NETCDF_LIBS}

test_adj : test_adj.o libMOAB.a
	${LINKER} ${DEBUG_FLAG} test_adj.o -o test_adj -L. -lMOAB ${NETCDF_LIBS}

# build the dependencies
depend : 
	${MAKEDEPEND} -DIS_BUILDING_MB ${NETCDF_INCLUDE} ${MB_LIB_SRCS} MBTest.cpp > make.dependencies
	${MAKEDEPEND} -DIS_BUILDING_MB ${NETCDF_INCLUDE} -DXPCOM_MB ${MB_LIB_SRCS} MBTest.cpp >> make.dependencies
	${MAKEDEPEND} -DIS_BUILDING_MB ${NETCDF_INCLUDE} -DMB_STATIC MBTest.cpp >> make.dependencies
	${MAKEDEPEND} -DIS_BUILDING_MB -DTEST HomXform.cpp >> make.dependencies
	${MAKEDEPEND} -DIS_BUILDING_MB scdseq_test.cpp >> make.dependencies

# clean up intermediate files
clean_all : clean
	@ cd TSTT; ${MAKE} clean_all

clean : clean_pcom clean_xpcom
	@ rm -f *.o *.o_test
	@ rm -rf ${TEMPLATE_DIR}
	@ rm -f moab_test test_exo test_tag_server test_ent_seq 
	@ rm -f *.a *.so

clean_xpcom : 
	@ rm -f *.xpcom.o
	@ rm -f moab_test.xpcom libMOABxpcom.so 

clean_pcom : 
	@ rm -f *.pcom.o
	@ rm -f moab_test.pcom libMOAB.so


.SUFFIXES : .o .cpp .c .o_test .xpcom.o .static.o

# .cpp.o rule for building regular object files
.cpp.o:
	${ECHO_COMMAND} 
	${PREFIX} ${CXX} -DIS_BUILDING_MB ${DEBUG_FLAG} ${MACH_CXXFLAGS} ${NETCDF_INCLUDE} -c -o $@ $<

# .cpp.static.o rule for building regular object files
.cpp.static.o:
	${ECHO_COMMAND} 
	${PREFIX} ${CXX} -DIS_BUILDING_MB -DMB_STATIC ${DEBUG_FLAG} ${MACH_CXXFLAGS} ${NETCDF_INCLUDE} -c -o $@ $<

# .cpp.o_test for building object files which have an individual component test in them
# meaning they have their own main() defined
.cpp.o_test:
	${ECHO_COMMAND}
	${PREFIX} ${CXX} -DIS_BUILDING_MB ${DEBUG_FLAG} -DTEST ${MACH_CXXFLAGS} ${NETCDF_INCLUDE} -c -o $@ $<

# .cpp.xpcom.o rule for building object files to work with xpcom
.cpp.xpcom.o:
	${ECHO_COMMAND} 
	${PREFIX} ${CXX} ${DEBUG_FLAG} ${MACH_CXXFLAGS} ${XPCOM_MACH_CXXFLAGS} ${NETCDF_INCLUDE} -c -o $@ $<

include make.dependencies
