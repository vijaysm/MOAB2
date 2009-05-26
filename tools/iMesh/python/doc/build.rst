=================
 Building PyTAPS
=================

Like most Python packages, PyTAPS uses Distutils for installation, so in general
setup consists simply of typing ``python setup.py install`` at the root
directory for PyTAPS. However, certain iMesh interfaces may require some
additional setup.

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