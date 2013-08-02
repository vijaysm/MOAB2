LDADD = $(top_builddir)/src/libMOAB.la

MESHDIR = $(top_srcdir)/MeshFiles/unittest

AM_CPPFLAGS += -DIS_BUILDING_MB \
               -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR) \
               -I. \
               -I.. -I$(srcdir)/.. \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/parallel \
               -I$(top_srcdir)/src/io/mhdf/include
               
if ENABLE_mbcoupler
  AM_CPPFLAGS += -I$(top_srcdir)/tools/mbcoupler \
                 -I$(top_srcdir)/itaps \
                 -I$(top_srcdir)/itaps/imesh \
                 -I$(top_builddir)/itaps/imesh \
                 -I$(top_srcdir)/src/moab/point_locater/lotte 
  
endif

if ENABLE_mbcslam
  AM_CPPFLAGS += -I$(top_srcdir)/tools/mbcslam
endif

# Run parallel tests in parallel
if USE_MPIEXEC
  OLD_ENVIRONMENT := ${TESTS_ENVIRONMENT}
  TESTS_ENVIRONMENT = ${MPIEXEC} ${MPIEXEC_NP} ${NP} ${OLD_ENVIRONMENT} 
endif

TESTS = pcomm_unit \
        parallel_unit_tests \
        uber_parallel_test \
        scdtest \
        pcomm_serial \
        $(NETCDF_TESTS) \
        $(HDF5_TESTS) \
        $(COUPLER_TESTS) \
        $(MBCSLAM_TESTS)

if PARALLEL_HDF5
  HDF5_TESTS = parallel_hdf5_test mhdf_parallel parallel_write_test
else
  HDF5_TESTS = 
endif

if NETCDF_FILE
  NETCDF_TESTS = scdpart read_nc_par ucdtrvpart mpastrvpart
else
  NETCDF_TESTS = 
endif

if PNETCDF_FILE
  if !NETCDF_FILE
    NETCDF_TESTS = scdpart read_nc_par ucdtrvpart mpastrvpart
  endif
endif

if ENABLE_mbcoupler
  COUPLER_TESTS = par_coupler_test
else
  COUPLER_TESTS = 
endif

if ENABLE_mbcslam
  MBCSLAM_TESTS = par_intx_sph 
else
  MBCSLAM_TESTS = 
endif

check_PROGRAMS = $(TESTS) mbparallelcomm_test partcheck structured3 parmerge

pcomm_unit_SOURCES = pcomm_unit.cpp
parallel_hdf5_test_SOURCES = parallel_hdf5_test.cc
mhdf_parallel_SOURCES = mhdf_parallel.c
parallel_unit_tests_SOURCES = parallel_unit_tests.cpp
parallel_write_test_SOURCES = parallel_write_test.cc
uber_parallel_test_SOURCES = uber_parallel_test.cpp
pcomm_serial_SOURCES = pcomm_serial.cpp
mbparallelcomm_test_SOURCES = mbparallelcomm_test.cpp
scdtest_SOURCES = scdtest.cpp
partcheck_SOURCES = partcheck.cpp
structured3_SOURCES = structured3.cpp
parmerge_SOURCES = parmerge.cpp
scdpart_SOURCES = scdpart.cpp
read_nc_par_SOURCES = ../io/read_nc.cpp
ucdtrvpart_SOURCES = ucdtrvpart.cpp
mpastrvpart_SOURCES = mpastrvpart.cpp

if ENABLE_mbcoupler
  par_coupler_test_SOURCES = par_coupler_test.cpp
  par_coupler_test_LDADD = $(LDADD) ../../tools/mbcoupler/libmbcoupler.la
endif

if ENABLE_mbcslam
  par_intx_sph_SOURCES = par_intx_sph.cpp
  par_intx_sph_LDADD = $(LDADD) ../../tools/mbcoupler/libmbcoupler.la ../../tools/mbcslam/libmbcslam.la  
endif
# Other files to clean up (e.g. output from tests)
MOSTLYCLEANFILES = mhdf_ll.h5m tmp0.h5m tmp1.h5m tmp2.h5m tmp3.h5m \
                   partial0.vtk \
                   partial1.vtk  \
                   partialConvex1.vtk \
                   partialConvex0.vtk \
                   initial0.vtk \
                   intersect1.h5m \
                   intersect0.h5m \
                   test.h5m \
                   initial1.vtk \
                   test_mpas.h5m \
                   dum.h5m
                   
