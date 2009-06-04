=================
 Building PyTAPS
=================

Like most Python packages, PyTAPS uses Distutils for installation, so in general
setup consists simply of typing ``python setup.py install`` at the root
directory for PyTAPS. However, certain iMesh interfaces may require some
additional setup.

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
