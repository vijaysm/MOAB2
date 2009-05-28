=================
 iMesh Interface
=================

.. module:: itaps

.. class:: itaps.iMesh([options])

   Return a new ``iMesh`` object with any implementation-specific options
   defined in ``options``.

   :param options: Implementation-specific options string

   .. attribute:: rootSet

      Get the handle of the root set for this instance. The entire mesh in this
      instance can be accessed from this set.

   .. attribute:: geometricDimension

      Get/set the geometric dimension of mesh represented in this instance.
      When setting the dimension, an application should not expect this
      function to succeed unless the mesh database is empty (no vertices
      created, no files read, etc.)

   .. attribute:: dfltStorage

      Get the default storage order used by this implementation.

   .. attribute:: adjTable

      Get the adjacency table for this implementation.  This table is a 4x4
      array, with indices 0-based, where A(i,j) (i=row, j=column) represents
      the relative cost of retrieving adjacencies between entities of dimension
      i to entities of dimension j.

   .. method:: load(filename[, options])

      Load a mesh from a file. Equivalent to ``rootSet.load(filename,
      options)``.

      :param filename: File name from which the mesh is to be loaded
      :param options: Implementation-specific options string

   .. method:: save(filename[, options])

      Save the mesh to a file. Equivalent to ``rootSet.load(filename,
      options)``.

      :param filename: File name to which the mesh is to be saved
      :param options: Implementation-specific options string

   .. method:: areEHValid(doReset)

      Return whether entity handles have changed since last reset or since
      instance construction. If true, it is not guaranteed that a handle from
      before the last call to this function represents the same entity as the
      same handle value does now. If ``doReset`` is true, resets the starting
      point for this function.

      :param doReset: If true, perform a reset on the starting point after
                      which handles are invariant.
      :return: True iff entity handles have changed

   .. method:: getVtxCoords(entities[, storageOrder])

      Get coordinates of specified vertices. If ``entitites`` is an array of
      entity handles, ``storageOrder`` is required. Otherwise, it is unused.

      :param entities: Entity or array of entities being queried
      :param storageOrder: Storage order of vertices to be returned
      :return: Array of vertices in the specified storage order. One-dimensional
               array if ``entities`` is a single element, two-dimesional
               otherwise

   .. method:: getEntType(entities)

      Get the entity type for the specified entities.

      :param entities: Entity or array of entities being queried
      :return: If ``entities`` is a single element, the type of the entity.
               Otherwise, an array of the entity types.

   .. method:: getEntTopo(entities)

      Get the entity topology for the specified entities.

      :param entities: Entity or array of entities being queried
      :return: If ``entities`` is a single element, the topology of the entity.
               Otherwise, an array of the entity topologies.

   .. method:: getEntAdj(entities, typeReq)

      Get entities of the specified type adjacent to elements of ``entities``.
      If ``entities`` is a single entity handle, returns an array of adjacent
      entities.

      If ``entities`` is an array of entities, returns a tuple type, with the
      first element being an array of offsets into the second element such that
      ``ret[0][ ret[1][i]:ret[1][i+1] ]`` is a list of entities adjacent to
      ``entities[i]``.

      :param entities: Entity or array of entities being queried
      :param typeReq: Type of adjacent entities being requested
      :return: If ``entities`` is a single element, an array of adjacent
               entities. Otherwise, a tuple containing an array of offsets and
               an array of adjacent entities.

   .. method:: getEnt2ndAdj(entities, bridgeType, typeReq)

      Get "2nd order" adjacencies to an array of entities, that is, from each 
      entity, through other entities of a specified "bridge" dimension, to
      other entities of another specified "to" dimension. If ``entities`` is a
      single entity handle, returns an array of adjacent entities.

      If ``entities`` is an array of entities, returns a tuple type, with the
      first element being an array of offsets into the second element such that
      ``ret[0][ ret[1][i]:ret[1][i+1] ]`` is a list of entities adjacent to
      ``entities[i]``.

      :param entities: Entity or array of entities being queried
      :param brideType: Type of bridge entity for 2nd order adjacencies
      :param typeReq: Type of adjacent entities being requested
      :return: If ``entities`` is a single element, an array of adjacent
               entities. Otherwise, a tuple containing an array of offsets and
               an array of adjacent entities.

   .. method:: createEntSet(isList)

      Create an entity set, either ordered (``isList == True``) or unordered 
      (``isList == False``). Unordered entity sets can contain a given entity
      or set only once.

      :param isList: True if the list should be ordered, false otherwise
      :return: The newly-created entity set

   .. method:: destroyEntSet(entSet)

      Destroy an entity set.

      :param entSet: Entity set to be destroyed

   .. method:: setVtxCoords(entities, coords[, storageOrder])

      Set the coordinates for the specified vertex or array of vertices. If
      ``entities`` is an array of vertices, ``storageOrder`` must be specified;
      otherwise it is ignored.

      :param entities: Vertex handle or array of vertex handles being set
      :param coords: New coordinates to assign to vertices
      :param storageOrder: Storage order of coordinates to be assigned

   .. method:: createVtx(coords[, storageOrder])

      Create a vertex or array of vertices with the specified coordinates. If
      creating multiple vertices, ``storageOrder`` must be specified; otherwise
      it is ignored.

      :param coords: Coordinates of new vertices to create
      :param storageOrder: Storage order of coordinates

   .. method:: createEnt(topo, entities)

      Create a new entity with the specified lower-order topology.

      :param topo: Topology of the entity to be created
      :param entities: Array of lower order entity handles used to construct
                       new entity
      :return: Tuple containing the created entity and its creation status

   .. method:: createEntArr(topo, entitites)

      Create an array of new entities with the specified lower-oder topology.

      :param topo: Topology of the entities to be created
      :param entities: Array of lower order entity handles used to construct
                       new entities
      :return: Tuple containing the created entities and their creation statuses

   .. method:: deleteEnt(entities)

      Delete the specified entity or array of entities.

      :param entities: An entity or array of entities to delete

   .. method:: createTag(name, size, type)

      Create a tag with specified ``name``, ``size``, and ``type``. The tag
      size is the number of values of type ``type`` that can be held. ``type``
      is one of the following:

      +-------+---------------+
      | ``i`` | Integer       |
      +-------+---------------+
      | ``d`` | Double        |
      +-------+---------------+
      | ``E`` | Entity handle |
      +-------+---------------+
      | ``b`` | Binary data   |
      +-------+---------------+

      :param name: Tag name
      :param size: Size of tag in number of values
      :param type: Character representing the tag's type
      :return: The created tag

   .. method:: destroyTag(tag, forced)

      Destroy a tag. If ``forced`` is true and entities still have values set
      for this tag, the tag is deleted anyway and those values disappear.
      Otherwise the tag is not deleted if entities still have values set for it.

      :param tag: Tag to delete
      :param forced: True if the tag should be deleted even if there are values
                     set for it

   .. method:: getTagHandle(name)

      Get the handle of an existing tag with the specified ``name``.

      :param name: The name of the tag to find
      :return: The tag with the specified name

   .. method:: getAllTags(entities)

      Get all the tags associated with a specified entity or entity set.

      :param entities: Entity or entity set being queried
      :return: Array of tags associated with ``entities``


.. class:: itaps.iMesh.topology

   An enumeration of mesh element topologies corresponding to
   ``iMesh_EntityTopology``.

   .. data:: point

      A general zero-dimensional entity

   .. data:: line_segment

      A general one-dimensional entity

   .. data:: polygon

      A general two-dimensional element

   .. data:: triangle

      A three-sided, two-dimensional element

   .. data:: quadrilateral

      A four-sided, two-dimensional element

   .. data:: polyhedron

      A general three-dimensional element

   .. data:: tetrahedron

      A four-sided, three-dimensional element whose faces are triangles

   .. data:: hexahedron

      A six-sided, three-dimensional element whose faces are quadrilaterals

   .. data:: prism

      A five-sided, three-dimensional element which has three quadrilateral
      faces and two triangular faces

   .. data:: pyramid

      A five-sided, three-dimensional element which has one quadrilateral face
      and four triangular faces

   .. data:: septahedron

      A hexahedral entity with one collapsed edge

   .. data:: all

      Allows the user to request information about all the topology types


.. class:: itaps.iMesh.iterator(set, type, topology[, size=1])

   Return a new iterator on the entity set ``set`` to iterate over entities of
   the specified ``type`` and ``topology``. If ``size`` is greater than 1, each
   step of the iteration will return an array of ``size`` entities. All
   entities of a given type or topology are requested by specifying 
   ``iBase.type.all`` or  `iMesh.topology.all``, respectively.

   :param set: Entity set to iterate over
   :param type: Type of entities being requested
   :param topo: Topology of entities being requested
   :param count: Number of entities to return on each step of iteration

   .. method:: reset()

      Resets the iterator to the beginning.


.. class:: itaps.iMesh.entitySet

   .. attribute:: isList

      Return whether this entity set is ordered.

   .. method:: load(entSet, filename[, options])

      Load a mesh from a file, adding it to this entity set.

      :param filename: File name from which the mesh is to be loaded
      :param options: Implementation-specific options string

   .. method:: save(filename[, options])

      Save the subset of the mesh contained in this entity set to a file.

      :param filename: File name to which the mesh is to be saved
      :param options: Implementation-specific options string

   .. method:: getNumOfType(type)

      Get the number of entities with the specified type in this entity set.

      :param type: Type of entity requested
      :return: The number of entities in entity set of the requested type

   .. method:: getNumOfTopo(topo)

      Get the number of entities with the specified topology in this entity set.

      :param type: Topology of entity requested
      :return: The number of entities in the entity set of the requested
               topology

   .. method:: getEntities(type, topo)

      Get entities of a specific type and/or topology in this entity set. All 
      entities of a given type or topology are requested by specifying
      ``iBase.type.all`` or ``iMesh.topology.all``, respectively.

      :param entSet: Entity set being queried
      :param type: Type of entities being requested
      :param topo: Topology of entities being requested
      :return: Array of entity handles from ``entSet`` meeting the requirements
               of ``type`` and ``topo``.

   .. method:: getAdjEntIndices(type, topo, adjType)

      Given an entity set and optionally a type or topology, return a tuple
      containing the following:

      * The entities in the set of the specified ``type`` and/or ``topology``
      * The entities adjacent to those entities with the specified type
        ``adjType``, as a list of unique handles
      * An index buffer containing, for each entity in the first list,
        the indices of the entities adjacent to it
      * An array of offsets into the index buffer for each entity in the first
        list

      That is, given an entity located in ``ret[0][i]``, the list of entities to
      which it is adjacent is::

        ret[1][  ret[2][ ret[3][i]:ret[3][i+1] ]  ]

      :param type: Type of entities being requested
      :param topo: Topology of entities being requested
      :param adjType: Type of adjacent entities being requested
      :return: 4-tuple containing the adjacency information

   .. method:: getNumEntSets(numHops)

      Get the number of sets contained in this entity set. If this entity set is
      not the root set, ``numHops`` indicates the maximum number of contained
      sets from ``self`` to one of the contained sets, inclusive of ``self``.

      :param numHops: Maximum number of contained sets from ``self`` to a
                      contained set, including ``self``.
      :return: Number of entity sets found

   .. method:: getEntSets(numHops)

      Get the sets contained in this entity set. If this entity set is not the
      root set, ``numHops`` indicates the maximum number of contained sets from
      ``self`` to one of the contained sets, inclusive of ``self``.

      :param numHops: Maximum number of contained sets from ``self`` to a
                      contained set, including ``self``.
      :return: Array of entity sets found      

   .. method:: add(entities)

      Add an entity, entity set, or array of entities to this entity set.

      :param entities: The entity, entity set, or array of entities to add

   .. method:: remove(entities)

      Remove an entity, entity set, or array of entities from this entity set.

      :param entities: The entity, entity set, or array of entities to remove

   .. method:: contains(entities)

      Return whether an entity, entity set, or array of entities is contained
      in this entity set.

      :param entities: The entity, entity set, or array of entities to query
      :return: If ``entities`` is an array of entities, an array of booleans
               corresponding to each element of ``entities``. Otherwise, a
               single boolean.

   .. method:: addChild(entSet)

      Add ``entSet`` as a child to this entity set.

      :param entSet: The entity set to add

   .. method:: removeChild(entSet)

      Remove ``entSet`` as a child from this entity set.

      :param entSet: The entity set to remove

   .. method:: isChild(entSet)

      Return whether an entity set is a child of this entity set.

      :param entSet: The entity set to query:
      :return: True if ``entSet`` is a child of this entity set, false otherwise

   .. method:: getNumChildren(numHops)

      Get the number of child sets linked from this entity set. If ``numHops``
      is non-zero, this represents the maximum hops from this entity set to any
      child in the count.

      :param numHops: Maximum hops from this entity set to a child set,
                      inclusive of the child set
      :return: Number of children

   .. method:: getNumParents(numHops)

      Get the number of parent sets linked from this entity set. If ``numHops``
      is non-zero, this represents the maximum hops from this entity set to any
      parents in the count.

      :param numHops: Maximum hops from this entity set to a parent set,
                      inclusive of the parent set
      :return: Number of parents

   .. method:: getChildren(numHops)

      Get the child sets linked from this entity set. If ``numHops`` is
      non-zero, this represents the maximum hops from this entity set to any
      child in the result.

      :param numHops: Maximum hops from this entity set to a child set,
                      inclusive of the child set
      :return: Array of children

   .. method:: getParents(numHops)

      Get the parents sets linked from this entity set. If ``numHops`` is
      non-zero, this represents the maximum hops from this entity set to any
      parent in the result.

      :param numHops: Maximum hops from this entity set to a parent set,
                      inclusive of the parent set
      :return: Array of parents

   .. method:: iterate(type, topo[, count=1])

      Initialize an iterator over the specified entity type and topology for
      this entity set. If ``count`` is greater than 1, each step of the
      iteration returns an array of ``count`` entities. Equivalent to::

        itaps.iMesh.iterator(self, type, topo, count)

      :param type: Type of entities being requested
      :param topo: Topology of entities being requested
      :param count: Number of entities to return on each step of iteration
      :return: An ``itaps.iMesh.iterator`` instance

   .. method:: difference(entSet)

      Subtract contents of an entity set from this set. Equivalent to
      ``self - entSet``.

      :param entSet: Entity set to subtract
      :return: Resulting entity set

   .. method:: intersection(entSet)

      Intersect contents of an entity set with this set. Equivalent to
      ``self & entSet``.

      :param entSet: Entity set to intersect
      :return: Resulting entity set

   .. method:: union(entSet)

      Unite contents of an entity set with this set. Equivalent to
      ``self | entSet``.

      :param entSet: Entity set to unite
      :return: Resulting entity set


.. class:: itaps.iMesh.tag

   .. attribute:: name

      Get the name for this tag.

   .. attribute:: sizeValues

      Get the size in number of values for this tag.

   .. attribute:: sizeBytes

      Get the size in bytes for this tag.

   .. attribute:: type

      Get the data type for this tag as a character code (see above).

   .. method:: setData(entities, data[, type])

      Set value(s) for the tag on an entity, entity set, or array of entities.
      If ``type`` is not specified, this function will retrieve the tag type
      automatically.

      :param entities: Entity, entity set, or array of entities on which tag is
                       being set
      :param data: Data to set
      :param type: Character representing the tag's type (as above)

   .. method:: getData(entities, [, type])

      Get value(s) for the tag on an entity, entity set, or array of entities.
      If ``type`` is not specified, this function will retrieve the tag type
      automatically.

      :param entities: Entity, entity set, or array of entities on which tag is
                       being retrieved
      :param type: Character representing the tag's type (as above)
      :return: The retrieved data

   .. method:: rmvTag(entities)

      Remove the tag value from an entity, entity set, or array of entities.

      :param entities: Entity, entity set, or array of entities from which tag
                       is being removed
