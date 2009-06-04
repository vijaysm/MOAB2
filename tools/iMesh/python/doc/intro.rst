==============
 Introduction
==============

From C
======

Compared to the C interface for ITAPS, PyTAPS has some important differences
designed to take advantage of Python language features. Users familiar with the
C interface should review the changes below.

Return Values
-------------

Instead of relying on out-values to return data from interface functions, PyTAPS
uses the simpler method of merely returning the data. In cases where multiple
return values are needed, PyTAPS returns a tuple of the values. As a notable
exception to this, return values from adjacency queries are wrapped into helper
classes: :mod:`itaps.helpers`.

Error Handling
--------------

All instances of the out-value ``int *err`` in functions have been replaced with
Python ``RuntimeError``\ s that are raised if an error of some kind occurs.

Arrays
------

In the C interface, arrays are passed as pairs of parameters when used as input
(``type *data, int size``) and as triples of parameters when used as output
(``type **data, int *alloc, int *size``). The PyTAPS interface uses `Numpy
<http://numpy.scipy.org/>`_ for all array data, and so the size is passed
implicitly along with the array. Since Python programs (generally) do not
manually manage their own memory, the ``alloc`` parameter is similarly
eliminated.

Overloads
---------

To reduce the number of methods that must be memorized when using these
interfaces, PyTAPS coalesces all entity, entity array, and entity set functions
into a single function where appropriate. For instance, the iMesh functions
``iMesh_addEntToSet``, ``iMesh_addEntArrToSet``, and ``iMesh_addEntSet`` are all
called by :meth:`itaps.iMesh.EntitySet.add`. Which C API function is ultimately
called is determined by the type of the argument passed to ``add``.

One notable exception to this rule is :meth:`itaps.iMesh.createEnt` and
:meth:`itaps.iMesh.createEntArr`, which have the same signatures and so cannot
be overloaded based on the types of the arguments.

An Example
==========

Below is a short (but complete) example of how one might use PyTAPS. The example
accepts two mesh files (one pre- and one post-deformation) and calculates the
volumes of the regions in each mesh. From there, it then determines the ratios
of the volumes of each region and graphs them (or optionally prints the raw data
to the console).

.. literalinclude :: ../tools/volume.py
