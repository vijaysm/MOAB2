=================
 iBase Interface
=================

.. module:: itaps.iBase

.. class:: itaps.iBase.Type

   An enumeration of entity types corresponding to ``iBase_EntityType``.

   .. data:: vertex

      A zero-dimensional entity

   .. data:: edge

      A one-dimensional entity

   .. data:: face

      A two-dimensional entity

   .. data:: region

      A three-dimensional entity

   .. data:: all

      Allows the user to request information about all the types


.. class:: itaps.iBase.AdjCost

   An enumeration of entity types corresponding to ``iBase_AdjacencyCost``.

   .. data:: unavailable

      Adjacency information not supported

   .. data:: all_order_1

      No more than local mesh traversal required

   .. data:: all_order_logn

      Global tree search

   .. data:: all_order_n

      Global exhaustive search

   .. data:: some_order_1

      Only some adjacency info, local

   .. data:: some_order_logn

      Only some adjacency info, tree

   .. data:: some_order_n

      Only some adjacency info, exhaustive


.. class:: itaps.iBase.StorageOrder

   An enumeration of entity types corresponding to ``iBase_StorageOrder``.

   .. data:: interleaved

      Coordinates are interleaved, e.g. ``[ x0, y0, z0, x1, y1, z1, ... ]``.

   .. data:: blocked

      Coordinates are blocked, e.g. ``[ x0, x1, ..., y0, y1, ..., z0, z1,
      ...]``.


.. class:: itaps.iBase.CreationStatus

   An enumeration of entity types corresponding to ``iBase_CreationStatus``.

   .. data:: new

      New entity was created

   .. data:: exists

      Entity already exists

   .. data:: duplicated

      Duplicate entity created

   .. data:: failed

      Creation failed