=================
 Building PyTAPS
=================

Like most Python packages, PyTAPS uses Distutils for installation, so in general
setup consists simply of typing ``python setup.py install`` at the root
directory for PyTAPS. However, certain iMesh interfaces may require some
additional setup.

The PyTAPS setup script supports importing definitions from the
`iMesh-Defs.inc` file. In order to make use of this, specify the root
location of your iMesh installation in the environment variable ``IMESHPATH``.
For example, if your `iMesh-Defs.inc` is located in
`/usr/local/iMesh/lib/iMesh-Defs.inc`, then ``IMESHPATH`` should be
`/usr/local/iMesh`. If the `iMesh-Defs.inc` file was not installed, or you
don't wish to use it, you can manually specify the build options as described
below in `Non-standard Library Locations`_.

Requirements
============

In order to build PyTAPS, several external libraries are required:

* Python 2.5+
* `Numpy <http://numpy.scipy.org/>`_ 1.3.0+
* An iMesh implementation (currently only `MOAB
  <http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>`_ is supported)

To run the performance tests or tools, `Matplotlib
<http://matplotlib.sourceforge.net/>`_ is also required. Finally, to build the
documentation (you're reading it!), `Sphinx <http://sphinx.pocoo.org/>`_ is
required.

Testing
=======

To run unit/regression tests on the package, you may specify ``test`` as an
argument to `setup.py`. The test command also accepts a verbosity level (an
integer), specified with ``-v`` or ``--verbosity``.

Building Documentation
======================

The documentation that you're currently reading can be built from `setup.py` by
specifying ``doc`` as an argument. This command supports the following options:

+-----------------------+--------+---------------------------------------+
| ``--builder=BUILDER`` | ``-b`` | documentation builder (default: html) |
+-----------------------+--------+---------------------------------------+
| ``--target=TARGET``   | ``-t`` | target directory for output           |
+-----------------------+--------+---------------------------------------+

With MOAB
=========

Since Python modules are loaded as shared libraries, MOAB's iMesh library needs
to be built with position-independent code. When configuring MOAB, specify
``--with-pic`` on the command line.

Non-standard Library Locations
==============================

In some cases, required objects for building PyTAPS aren't in the expected
directories. One solution to this is to include the appropriate directories in
the environment variables ``CPATH`` and ``PYTHONPATH``. Another, more flexible
method is to use the ``setup.cfg`` file. Information on how to use this file can
be found in the official Python `documentation <http://docs.python.org/install/index.html#distutils-configuration-files>`_.
