##############################################################################
# Zoltan Library for Parallel Applications                                   #
# Copyright (c) 2000,2001,2002, Sandia National Laboratories.                #
# For more info, see the README file in the top-level Zoltan directory.      # 
##############################################################################
##############################################################################
# CVS File Information
#    $RCSfile$
#    $Author$
#    $Date$
#    $Revision$
##############################################################################

##############################################################################
# Makefile for simple example which uses the Zoltan and MOAB libraries
#
# To compile the example, type
#    gmake 
##############################################################################

TARGETS = mbzoltanLibrary mbzoltan
SUBDIRS = lib mbzoltan
DISTFILES = Makefile Makefile.zoltan README
DIST_SUBDIRS = lib mbzoltan
PATH_SEPARATOR = :
am__cd = CDPATH="$${ZSH_VERSION+.}$(PATH_SEPARATOR)" && cd

all: $(TARGETS)

mbzoltanLibrary:
	@$(MAKE) -f Makefile.zoltan $@

mbzoltan: 
	@$(MAKE) -f Makefile.zoltan $@

install: mbzoltan
	@$(MAKE) -f Makefile.zoltan $@

distdir: $(DISTFILES)
	@list='$(DISTFILES)'; for file in $$list; do \
	  cp -p ./$$file $(distdir)/$$file \
	  || exit 1; \
	done
	list='$(DIST_SUBDIRS)'; for subdir in $$list; do \
	  mkdir -p -- "$(distdir)/$$subdir" \
	  || exit 1; \
	  distdir=`$(am__cd) $(distdir) && pwd`; \
	  top_distdir=`$(am__cd) $(top_distdir) && pwd`; \
	  (cd $$subdir && \
	    $(MAKE) $(AM_MAKEFLAGS) \
	      top_distdir="$$top_distdir" \
	      distdir="$$distdir/$$subdir" \
	      distdir) \
	    || exit 1; \
	done

clean:
	@cd lib; $(MAKE) clean
	@cd mbzoltan; $(MAKE) clean

distclean: clean
	rm -f Config.moab

