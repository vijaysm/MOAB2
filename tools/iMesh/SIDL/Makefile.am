SUBDIRS = . mserver Cclient Fclient

# $(RECURSIVE_TARGETS): mserver/babel.make Cclient/babel.make Fclient/babel.make
# 
# BABEL = @BABEL_DIR@/bin/babel
# if USE_BABEL
#   SIDL_FILES = $(srcdir)/iBase.sidl $(srcdir)/iMesh.sidl $(srcdir)/iMesh_SIDL.sidl
#   server_babel_cmd = $(BABEL) -sC++ -omserver 
#   cclient_babel_cmd = $(BABEL) -cC -oCclient
#   fclient_babel_cmd = $(BABEL) -cF77 -oFclient
# else
#   SIDL_FILES=
#   server_babel_cmd = touch mserver/babel.make
#   cclient_babel_cmd = touch Cclient/babel.make
#   fclient_babel_cmd = touch Fclient/babel.make
# endif
# 
# mserver/babel.make: $(SIDL_FILES)
# 	$(server_babel_cmd) $(SIDL_FILES)
# 
# Cclient/babel.make: $(SIDL_FILES)
# 	$(cclient_babel_cmd) $(SIDL_FILES)
# 
# Fclient/babel.make: $(SIDL_FILES)
# 	$(fclient_babel_cmd) $(SIDL_FILES)

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = iMesh-SIDL-Defs.inc

EXTRA_DIST = iBase.sidl iMesh.sidl iMesh_SIDL.sidl iBase_SIDL_defs.h

# By default, iMesh-SIDL-Defs.inc will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does 
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
install-data-hook:
	echo "iMesh_SIDL_includes=-I$(includedir) -I$(includedir)/SIDL" >> $(DESTDIR)$(cfgdir)/iMesh-SIDL-Defs.inc
	echo "iMesh_SIDL_LFDFLAGS=-L${libdir}" >> $(DESTDIR)$(cfgdir)/iMesh-SIDL-Defs.inc
