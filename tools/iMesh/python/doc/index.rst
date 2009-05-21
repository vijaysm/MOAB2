========
 PyTAPS
========

About PyTAPS
============

PyTAPS is a Python package designed to interface with ITAPS, focusing on the
iMesh interface.

Building PyTAPS
===============

Like most Python packages, PyTAPS uses Distutils for installation, so in general
setup consists simply of typing ``python setup.py install`` at the root
directory for PyTAPS. However, certain iMesh interfaces may require some
additional setup.

With MOAB
---------

Since Python modules are loaded as shared libraries, MOAB's iMesh library needs
to be built with position-independent code. When configuring MOAB, specify
``--with-pic`` on the command line.

