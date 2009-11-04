================
 PyTAPS Helpers
================

.. module:: itaps.helpers
   :synopsis: Helper classes to simplify common operations.

AdjacencyList
=============

.. class:: AdjacencyList(adj, offsets)

   .. attribute:: adj

      A one-dimensional array of entities adjacent to the queried entities

   .. attribute:: offsets

      An array of offsets into :attr:`adj` for each of the queried entities

   .. method:: __getitem__(i[, j])

      Return the entities adjacent to the `i`\ th entity. If `j` is specified,
      returns only the `j`\ th entity of the preceding array. As you would
      expect, this is generally called as ``list[i, j]``.

      .. note::
         This method is equivalent to ``adj[ offsets[i]:offsets[i+1] ]`` when
         only `i` is specified, and is equivalent to ``adj[ offsets[i]+j ]``
         when both `i` and `j` are specified.

      :param i: Index of the entity to query for adjacencies
      :param j: Index into the `i`\ th entity's adjacency array
      :return: If `j` is specified, a single entity. Otherwise, an array of
               entities.

   .. method:: length([i])

      Return the number of entities whose adjacencies are stored in this object.
      If `i` is specified, return the number of adjacencies for the `i`\ th
      entity.

      :param i: Index of the entity to query
      :return: If `i` is ``None``, the number of entities whose adjacencies
               are stored. Otherwise, the number of adjacencies for the
               `i`\ th entity.


IndexedAdjacencyList
====================

.. class:: IndexedAdjacencyList(entities, adj, indices, offsets)

   .. attribute:: entities

      A one-dimensional array of entities

   .. attribute:: adj

      A one-dimensional array of all entities adjacent to the elements of
      :attr:`entities`

   .. attribute:: indices

      An index buffer into :attr:`adj`

   .. attribute:: offsets

      An array of offsets into :attr:`indices` for each of the queried entities

   .. method:: __getitem__(i[, j])

      Return the entities adjacent to the `i` th entity. If `j` is specified,
      returns only the `j`\ th entity of the preceding array. As you would
      expect, this is generally called as ``list[i, j]``.

      .. note::
         This method is equivalent to ``adj[ index(i, j) ]``, and relies on the
         special indexing features of Numpy arrays.

      :param i: Index of the entity to query for adjacencies
      :param j: Index into the `i` th entity's adjacency array
      :return: If `j` is specified, a single entity. Otherwise, an array of
               entities.

   .. method:: index(i[, j])

      Return the indices of the entities adjacent to the `i`\ th entity. If `j`
      is specified, returns only the `j`\ th index of the preceding array.

      :param i: Index of the entity to query for adjacencies
      :param j: Index into the `i`\ th entity's adjacency array
      :return: If `j` is specified, a single index. Otherwise, an array of
               indices.

      .. note::
         This method is equivalent to ``indices[ offsets[i]:offsets[i+1] ]``
         when only `i` is specified, and is equivalent to
         ``indices[ offsets[i]+j ]`` when both `i` and `j` are specified.

   .. method:: length([i])

      Return the number of entities whose adjacencies are stored in this object.
      If `i` is specified, return the number of adjacencies for the `i`\ th
      entity.

      :param i: Index of the entity to query
      :return: If `i` is ``None``, the number of entities whose adjacencies
               are stored. Otherwise, the number of adjacencies for the
               `i`\ th entity.
