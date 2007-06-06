if ENABLE_mbconvert
  converter_DIR = converter
else
  converter_DIR=
endif

if ENABLE_hexmodops
  hexmodops_DIR = hexmodops
else
  hexmodops_DIR=
endif
if ENABLE_mbperf
  mbperf_DIR = mbperf
else
  mbperf_DIR =
endif
if ENABLE_mbchaco
  mbchaco_DIR = mbchaco
else
  mbchaco_DIR=
endif
if ENABLE_mbzoltan
  mbzoltan_DIR = mbzoltan
else
  mbzoltan_DIR=
endif

if ENABLE_qvdual
  qvdual_DIR = qvdual
else
  qvdual_DIR=
endif

if ENABLE_mbsize
  size_DIR = size
else
  size_DIR=
endif

if ENABLE_mbskin
  skin_DIR = skin
else
  skin_DIR=
endif

if ENABLE_mbtagprop
  tagprop_DIR = tagprop
else
  tagprop_DIR=
endif

if ENABLE_spheredecomp
  spheredecomp_DIR = spheredecomp
else
  spheredecomp_DIR=
endif

if ENABLE_mbsurfplot
  mbsurfplot_DIR = surfplot
else
  mbsurfplot_DIR=
endif

if ENABLE_dagmc
  dagmc_DIR = dagmc
else
  dagmc_DIR = 
endif 

if ENABLE_imesh
  imesh_DIR = iMesh
else
  imesh_DIR = 
endif 

SUBDIRS = $(converter_DIR) \
          $(hexmodops_DIR) \
          $(mbchaco_DIR) \
          $(mbzoltan_DIR) \
          $(mbperf_DIR) \
          $(qvdual_DIR) \
          $(size_DIR) \
          $(skin_DIR) \
          $(tagprop_DIR) \
          $(spheredecomp_DIR) \
          $(mbsurfplot_DIR) \
          $(dagmc_DIR) \
          $(imesh_DIR)

EXTRA_DIST = mbzoltan/README 
