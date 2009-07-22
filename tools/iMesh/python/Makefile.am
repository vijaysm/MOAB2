BUILD_EXT_FLAGS = -I"$(PYTHON_DIR):$(NUMPY_DIR):$(abs_top_builddir)" \
                  --imesh-dir="$(abs_builddir)/.."

INSTALL_FLAGS = --skip-build \
                --prefix="$(prefix)" \
                --exec-prefix="$(exec_prefix)" \
                --install-lib="$(libdir)" \
                --install-headers="$(includedir)" \
                --install-data="$(datadir)" \
                --install-scripts="$(bindir)"

all:
	$(PYTHON) $(abs_srcdir)/setup.py build_ext $(BUILD_EXT_FLAGS)
	$(PYTHON) $(abs_srcdir)/setup.py build

install: all
	$(PYTHON) $(abs_srcdir)/setup.py install $(INSTALL_FLAGS)

clean:
	$(PYTHON) $(abs_srcdir)/setup.py clean --all

EXTRA_DIST = doc perf pkg test tools \
             common.h \
             errors.h \
             iBase.c \
             iBase_Python.h \
             iMesh.c \
             iMesh_entSet.inl \
             iMesh_iter.inl \
             iMesh_Python.h \
             iMesh_tag.inl \
             Makefile.am \
             numpy_extensions.h \
             numpy_extensions.inl \
             setup.cfg.example \
             setup.py
