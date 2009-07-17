=================
 Building PyTAPS
=================

Requirements
============

In order to build PyTAPS, several external libraries are required:

* Python 2.5+
* `Numpy <http://numpy.scipy.org/>`_ 1.3.0+
* An iMesh implementation (currently only `MOAB
  <http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>`_ is explicitly supported)

To run the performance tests or tools, `Matplotlib
<http://matplotlib.sourceforge.net/>`_ is also required. Finally, to build the
documentation (you're reading it!), `Sphinx <http://sphinx.pocoo.org/>`_ is
required.

Building Within MOAB
====================

If you downloaded PyTAPS as a part of MOAB, you can use the standard UNIX
``./configure; make; make install`` process to build it. To do so, add
``--enable-pytaps`` to your configure command. If your Python or Numpy libraries
are in a non-standard location, also specify ``--with-python=DIR`` and
``--with-numpy=DIR``, respectively.

Building Independently
======================

Like most Python packages, PyTAPS uses Distutils for installation, so in general
setup consists simply of typing ``python setup.py install`` at the root
directory for PyTAPS. However, certain iMesh interfaces may require some
additional setup.

The PyTAPS setup script supports importing definitions from the
`iMesh-Defs.inc` file. In order to make use of this, specify the command-line
options ``--imesh-dir=PATH`` to the `build_ext` and `perf_build` commands. For
example, if your `iMesh-Defs.inc` is located in
`/usr/local/iMesh/lib/iMesh-Defs.inc`, then ``--imesh-dir`` should be
`/usr/local/iMesh/lib`.  This options may also be specified in the `setup.cfg`
file:

.. literalinclude:: ../setup.cfg.example
   :language: ini

If the `iMesh-Defs.inc` file was not installed, or you don't wish to use it,
you can manually specify the build options as described below in `Non-standard
Library Locations`_.

Testing
-------

To run unit/regression tests on the package, you may specify ``test`` as an
argument to `setup.py`. The test command also accepts a verbosity level (an
integer), specified with ``-v`` or ``--verbosity``.

Performance Testing
-------------------

To run performance tests on the package comparing the speed of a pure-C
usage of ITAPS with PyTAPS, you may pass ``perf`` as a command to `setup.py`.
The performance test requires an input file to test. You may also specify the
number of times to repeat the tests::

  python setup.py perf --file=/path/to/file --count=N

The performance tests consist of two sub-commands, described below.

Building the Performance Test
_____________________________

This command builds the C portion of the performance tests. The following
options are allowed:

+--------------------------+--------+---------------------------------------+
| ``--include-dirs=PATHS`` | ``-I`` | list of directories to search for     |
|                          |        | include files                         |
+--------------------------+--------+---------------------------------------+
| ``--library-dirs=PATHS`` | ``-L`` | list of directories to search for     |
|                          |        | library files                         |
+--------------------------+--------+---------------------------------------+
| ``--libraries=LIBS``     | ``-l`` | list of libraries to link to          |
+--------------------------+--------+---------------------------------------+
| ``--imesh-dir=PATH``     | N/A    | root location of iMesh libraries      |
+--------------------------+--------+---------------------------------------+

Running the Performance Test
____________________________

This command executes the performance tests. The following options are allowed
(as with ``perf``):

+-----------------+--------+-------------------------------------------+
| ``--file=PATH`` | ``-F`` | file or directory containing test file(s) |
+-----------------+--------+-------------------------------------------+
| ``--count=N``   | ``-c`` | number of times to repeat each test       |
+-----------------+--------+-------------------------------------------+

Building Documentation
----------------------

The documentation that you're currently reading can be built from `setup.py` by
specifying ``doc`` as an argument. This command supports the following options:

+-----------------------+--------+---------------------------------------+
| ``--builder=BUILDER`` | ``-b`` | documentation builder (default: html) |
+-----------------------+--------+---------------------------------------+
| ``--target=TARGET``   | ``-t`` | target directory for output           |
+-----------------------+--------+---------------------------------------+

With MOAB
---------

Since Python modules are loaded as shared libraries, MOAB's iMesh library needs
to be built with position-independent code. When configuring MOAB, specify
``--with-pic`` on the command line.

Non-standard Library Locations
------------------------------

In some cases, required objects for building PyTAPS aren't in the expected
directories. One solution to this is to include the appropriate directories in
the environment variables ``CPATH`` and ``PYTHONPATH``. Another, more flexible
method is to use the `setup.cfg` file. Information on how to use this file can
be found in the official Python `documentation
<http://docs.python.org/install/index.html#distutils-configuration-files>`_.
