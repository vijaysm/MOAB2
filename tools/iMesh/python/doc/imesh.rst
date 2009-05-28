=================
 iMesh Interface
=================

.. module:: itaps

.. class:: itaps.iMesh

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

   .. method:: getEntities(entSet, type, topo)

      Get entities of a specific type and/or topology in ``entSet``. All 
      entities of a given type or topology are requested by specifying
      ``iBase.type.all`` or ``iMesh.topology.all``, respectively.

      :param entSet: Entity set being queried
      :param type: Type of entities being requested
      :param topo: Topology of entities being requested
      :return: Array of entity handles from ``entSet`` meeting the requirements
               of ``type`` and ``topo``.      

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

   .. method:: getAdjEntIndices(entSet, typeRequestor, topoRequestor, typeRequested)

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

      +---+---------------+
      | i | Integer       |
      +---+---------------+
      | d | Double        |
      +---+---------------+
      | E | Entity handle |
      +---+---------------+
      | b | Binary data   |
      +---+---------------+

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

   .. method:: setData(entities, tag, data[, type])

      Set value(s) for a tag on an entity, entity set, or array of entities.
      If ``type`` is not specified, this function will retrieve the tag type
      automatically.

      :param entities: Entity, entity set, or array of entities on which tag is
                       being set
      :param tag: Tag being set
      :param data: Data to set
      :param type: Character representing the tag's type (as above)

   .. method:: getData(entities, tag[, type])

      Get value(s) for a tag on an entity, entity set, or array of entities.
      If ``type`` is not specified, this function will retrieve the tag type
      automatically.

      :param entities: Entity, entity set, or array of entities on which tag is
                       being retrieved
      :param tag: Tag being retrieved
      :param type: Character representing the tag's type (as above)
      :return: The retrieved data

   .. method:: getAllTags(entities)

      Get all the tags associated with a specified entity or entity set.

      :param entities: Entity or entity set being queried
      :return: Array of tags associated with ``entities``

   .. method:: rmvTag(entities, tag)

      Remove a tag value from an entity, entity set, or array of entities.

      :param entities: Entity, entity set, or array of entities from which tag
                       is being removed
      :param tag: Tag to be removed


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


.. class:: itaps.iMesh.iterator

   .. method:: reset()


.. class:: itaps.iMesh.entitySet

   .. attribute:: isList

      Returns whether the entity set is ordered.

   .. method:: load(entSet, filename[, options])

      Load a mesh from a file, adding it to the entity set.

      :param filename: File name from which the mesh is to be loaded
      :param options: Implementation-specific options string

   .. method:: save(filename[, options])

      Save the subset of the mesh contained in the entity set to a file.

      :param filename: File name to which the mesh is to be saved
      :param options: Implementation-specific options string

   .. method:: getNumOfType(type)

      Get the number of entities with the specified type in the entity set.

      :param type: Type of entity requested
      :return: The number of entities in entity set of the requested type

   .. method:: getNumOfTopo(topo)

      Get the number of entities with the specified topology in the entity set.

      :param type: Topology of entity requested
      :return: The number of entities in the entity set of the requested
               topology

   .. method:: getNumEntSets(numHops)

   .. method:: getEntSets(numHops)

   .. method:: add(entities)

   .. method:: remove(entities)

   .. method:: contains(entities)

   .. method:: addChild(entSet)

   .. method:: removeChild(entSet)

   .. method:: isChild(entSet)

   .. method:: getNumChildren(numHops)

   .. method:: getNumParents(numHops)

   .. method:: getChildren(numHops)

   .. method:: getParents(numHops)

   .. method:: iterate(type, topo[, count=1])

   .. method:: difference(entSet)

   .. method:: intersection(entSet)

   .. method:: union(entSet)


.. class:: itaps.iMesh.tag

   .. attribute:: name

   .. attribute:: sizeValues

   .. attribute:: sizeBytes

   .. attribute:: type