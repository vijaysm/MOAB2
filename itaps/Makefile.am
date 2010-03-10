
if ENABLE_imesh
  imesh_DIR = imesh
else
  imesh_DIR = 
endif 

if ENABLE_igeom
  igeom_DIR = igeom
else
  igeom_DIR = 
endif 

include_HEADERS = iBase.h iBase_f.h

SUBDIRS = $(imesh_DIR) $(igeom_DIR)

EXTRA_DIST = mkprotos.sh
