# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign subdir-objects

LDADD = $(top_builddir)/src/libMOAB.la

MESHDIR = $(abs_top_srcdir)/MeshFiles/unittest

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
  AM_CPPFLAGS += -I$(top_srcdir)/tools/mbcoupler
endif

if ENABLE_mbcslam
  AM_CPPFLAGS += -I$(top_srcdir)/tools/mbcslam
endif

# Run parallel tests in parallel
if USE_MPIEXEC
  LOG_COMPILER = ${MPIEXEC}
  AM_LOG_FLAGS = ${MPIEXEC_NP} ${NP}
endif

TESTS = pcomm_unit \
        scdtest \
        pcomm_serial \
	par_spatial_locator_test \
        $(NETCDF_TESTS) \
        $(HDF5_TESTS) \
        $(MBCSLAM_TESTS) $(IMESH_TESTS)

if HAVE_HDF5_PARALLEL
  HDF5_TESTS = parallel_hdf5_test mhdf_parallel parallel_write_test \
        parallel_unit_tests \
        uber_parallel_test parallel_adj parmerge_test augment_with_ghosts
else
  HDF5_TESTS = 
endif

if HAVE_NETCDF
  NETCDF_TESTS = scdpart read_nc_par ucdtrvpart mpastrvpart gcrm_par
else
  NETCDF_TESTS = 
endif

if HAVE_PNETCDF
if !HAVE_NETCDF
  NETCDF_TESTS += scdpart read_nc_par ucdtrvpart mpastrvpart write_nc_par
else
  NETCDF_TESTS += write_nc_par
endif
endif


if ENABLE_mbcslam
if HAVE_HDF5_PARALLEL
  MBCSLAM_TESTS = par_intx_sph 
else
  MBCSLAM_TESTS =
endif
else
  MBCSLAM_TESTS = 
endif

IMESH_TESTS = 
if ENABLE_imesh
if HAVE_HDF5_PARALLEL
if ENABLE_FORTRAN
   IMESH_TESTS += imeshp_test
   AM_FCFLAGS += -I$(top_srcdir)/itaps \
                 -I$(top_srcdir)/itaps/imesh \
                 -I$(top_builddir)/itaps \
                 -I$(top_builddir)/itaps/imesh
endif
endif 
endif

check_PROGRAMS = $(TESTS) mbparallelcomm_test partcheck structured3 
if HAVE_HDF5_PARALLEL
  check_PROGRAMS += parmerge
endif


pcomm_unit_SOURCES = pcomm_unit.cpp
parallel_hdf5_test_SOURCES = parallel_hdf5_test.cc
mhdf_parallel_SOURCES = mhdf_parallel.c
mhdf_parallel_LDADD = $(LDADD) $(HDF5_LIBS)
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
gcrm_par_SOURCES = gcrm_par.cpp
write_nc_par_SOURCES = ../io/write_nc.cpp
par_spatial_locator_test_SOURCES = par_spatial_locator_test.cpp
parallel_adj_SOURCES = ../adj_moab_test.cpp
parmerge_test_SOURCES = parmerge_test.cpp
augment_with_ghosts_SOURCES = augment_with_ghosts.cpp

if ENABLE_imesh
if HAVE_HDF5_PARALLEL
if ENABLE_FORTRAN
  imeshp_test_SOURCES = imeshp_test.F90
  imeshp_test_LDADD =   $(top_builddir)/itaps/imesh/libiMesh.la $(LDADD) ${MOAB_CXX_LINKFLAGS} ${MOAB_CXX_LIBS}
endif
endif
endif


if ENABLE_mbcslam
  par_intx_sph_SOURCES = par_intx_sph.cpp
  par_intx_sph_LDADD = $(LDADD) ../../tools/mbcslam/libmbcslam.la  
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
                   dum.h5m \
                   test_gcrm.h5m \
                   test_gcrm_rcbzoltan.h5m \
                   test_mpas_no_mixed_elements.h5m \
                   test_mpas_rcbzoltan.h5m \
                   test_mpas_rcbzoltan_no_mixed_elements.h5m \
                   test_par_eul_across_files.nc \
                   test_par_eul_append.nc \
                   test_par_eul_ghosting.nc \
                   test_par_eul_T2.nc \
                   test_par_eul_T.nc \
                   test_par_fv_T.nc \
                   test_par_gcrm_vars.nc \
                   test_par_homme_T.nc \
                   test_par_mpas_vars.nc \
                   test2.h5m  \
                   testpm.h5m \
                   newfile.h5m
