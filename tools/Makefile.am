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

SUBDIRS = $(converter_DIR) \
          $(hexmodops_DIR) \
          $(mbchaco_DIR) \
          $(mbperf_DIR) \
          $(qvdual_DIR) \
          $(size_DIR) \
          $(skin_DIR) \
          $(tagprop_DIR) \
          $(spheredecomp_DIR)

# For targets that require a Makefile even for tools that 
# are not enabled, and if that tool has it's own configure
# script, create the Makefile if necessary
.PHONY: tool-makefiles

distdir mostlyclean-recursive clean-recursive distclean-recursive \
maintainer-clean-recursive: tool-makefiles

tool-makefiles: qvdual/Makefile mbchaco/Makefile

qvdual/Makefile: qvdual/Makefile.in qvdual/config.status
	(cd qvdual && ./config.status Makefile)
qvdual/config.status: qvdual/configure
	eval `grep -e ^ac_configure_extra_args= $(top_builddir)/config.status`; \
	cd qvdual; \
	./configure $$ac_configure_extra_args
qvdual/Makefile.in: qvdual/configure.in qvdual/Makefile.am
	(cd qvdual && $(AUTOMAKE))
qvdual/configure: qvdual/configure.in qvdual/aclocal.m4 qvdual/acsite.m4
	(cd qvdual && $(AUTOCONF))

mbchaco/Makefile: mbchaco/Makefile.in mbchaco/config.status
	(cd mbchaco && ./config.status Makefile)
mbchaco/config.status: mbchaco/configure
	eval `grep -e ^ac_configure_extra_args= $(top_builddir)/config.status`; \
	cd mbchaco; \
	./configure $$ac_configure_extra_args
mbchaco/Makefile.in: mbchaco/configure.in mbchaco/Makefile.am
	(cd mbchaco && $(AUTOMAKE))
mbchaco/configure: mbchaco/configure.in mbchaco/aclocal.m4 mbchaco/acsite.m4
	(cd mbchaco && $(AUTOCONF))
