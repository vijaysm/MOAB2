BUILD_EXT_FLAGS = --include-dirs="$(PYTHON_DIR):$(NUMPY_DIR):$(abs_top_srcdir)"\
                  --imesh-dir="$(abs_srcdir)/.."

INSTALL_FLAGS = --skip-build \
                --prefix="$(prefix)" \
                --exec-prefix="$(exec_prefix)" \
                --root="$(DESTDIR)" \
                --install-lib="$(libdir)" \
                --install-headers="$(incluedir)" \
                --install-data="$(datadir)" \
                --install-scripts="$(bindir)"

all:
	python setup.py build_ext $(BUILD_EXT_FLAGS)
	python setup.py build

install:
	python setup.py install $(INSTALL_FLAGS)

clean:
	python setup.py clean --all
