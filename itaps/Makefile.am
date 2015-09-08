
if ENABLE_imesh
  imesh_DIR = imesh
else
  imesh_DIR = 
endif 

if ENABLE_fbigeom
  fbigeom_DIR = fbigeom
else
  fbigeom_DIR = 
endif 

include_HEADERS = iBase.h iBase_f.h

SUBDIRS = $(imesh_DIR) $(fbigeom_DIR)

EXTRA_DIST = mkprotos.sh
