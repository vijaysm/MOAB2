if ENABLE_imesh
  IMESH_TEST = imesh
else
  IMESH_TEST =
endif

if HAVE_CGM
  IGEOM_TEST = igeom
else
  IGEOM_TEST =
endif

if HAVE_NETCDF
  EXO_TEST = convert 
else
  EXO_TEST = 
endif

if HAVE_CPPUNIT
  UNIT_TEST = unit
else
  UNIT_TEST =
endif

if HAVE_MPI
  PARALLEL_UNTANGLE_SHAPE = parallel_untangle_shape
  PARALLEL_SMOOTH_LAPLACE = parallel_smooth_laplace
else
  PARALLEL_UNTANGLE_SHAPE = 
  PARALLEL_SMOOTH_LAPLACE = 
endif

MESHDIR = $(abs_top_srcdir)/MeshFiles/mesquite

AM_CPPFLAGS  = -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR) \
               $(MSQ_INCLUDES) \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/test

LDADD = $(top_builddir)/src/mesquite/libmesquite.la
check_PROGRAMS = algorithm_test activeset analytical_grad3D paraboloid_domain high_aspect_ratio idft_time jacobi
check_PROGRAMS += laplacian laplacian_polygon nongradient pyramid wedge simple_hybrid synchronous transform untangle
check_PROGRAMS += headers_basic_test benchmark_tests feasible_newton 2d_target 2d_metrics 2d_formulation 

headers_basic_test_SOURCES = headers_basic_test.cpp
2d_target_SOURCES = 2d_target_test.cpp
2d_metrics_SOURCES = 2d_metrics_test.cpp
2d_formulation_SOURCES = 2d_formulation_test.cpp
algorithm_test_SOURCES = algorithm_test.cpp
activeset_SOURCES = activeset_test.cpp
analytical_grad3D_SOURCES = analytical_grad3D_test.cpp
benchmark_tests_SOURCES = benchmark_tests.cpp
feasible_newton_SOURCES = feasible_newton_test.cpp
high_aspect_ratio_SOURCES = high_aspect_ratio_test.cpp
idft_time_SOURCES = idft_time_test.cpp
jacobi_SOURCES = jacobi_test.cpp
laplacian_SOURCES = laplacian_test.cpp
laplacian_polygon_SOURCES = laplacian_polygon_test.cpp
nongradient_SOURCES = nongradient_test.cpp
paraboloid_domain_SOURCES = paraboloid_domain_test.cpp
pyramid_SOURCES = pyramid_test.cpp
simple_hybrid_SOURCES = simple_hybrid_test.cpp
synchronous_SOURCES = synchronous_test.cpp
transform_SOURCES = transform_test.cpp
untangle_SOURCES = untangle_test.cpp
wedge_SOURCES = wedge_test.cpp

TESTS = $(check_PROGRAMS)
SUBDIRS = . \
    $(EXO_TEST) \
    higher_order \
    $(IGEOM_TEST) \
    $(IMESH_TEST) \
    $(PARALLEL_SMOOTH_LAPLACE) \
    $(PARALLEL_UNTANGLE_SHAPE) \
    slaved \
    wrapper_tests\
    $(UNIT_TEST)

EXTRA_DIST = 2d_formulation_initial.vtk 2d_metrics_reference.vtk msqcppflags.make
MOSTLYCLEANFILES = original_mesh.vtk smoothed_mesh.vtk smooth_ho.vtk *results.vtk transformed_mesh.vtk *.gpt Exp*.vtk out.vtk \
									 *.gnuplot synchronous.vtk high_aspect.vtk msqcppflags.make

include $(top_srcdir)/src/mesquite/make.inc

