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

   .. method:: load(entSet, filename[, options])

      Load a mesh from a file.

      :param entSet: Set to which loaded mesh will be added, root set if not
                      desired
      :param filename: File name from which the mesh is to be loaded
      :param options: Implementation-specific options string

   .. method:: save(entSet, filename[, options])

      Save the mesh to a file.

      :param entSet: Save a mesh to a file. If entity set is specified, save
                     only the mesh contained in that set.
      :param filename: File name to which the mesh is to be saved
      :param options: Implementation-specific options string

   .. method:: getNumOfType(entSet, type)

      Get the number of entities with the specified type in ``entSet``.

      :param entSet: Entity set being queried
      :param type: Type of entity requested
      :return: The number of entities in ``entSet`` of the requested type

   .. method:: getNumOfTopo(entSet, topo)

      Get the number of entities with the specified topology in ``entSet``.

      :param entSet: Entity set being queried
      :param type: Topology of entity requested
      :return: The number of entities in ``entSet`` of the requested topology

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
      ``ret[1][ ret[0][i]:ret[0][i+1] ]`` is a list of entities adjacent to
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
      ``ret[1][ ret[0][i]:ret[0][i+1] ]`` is a list of entities adjacent to
      ``entities[i]``.

      :param entities: Entity or array of entities being queried
      :param brideType: Type of bridge entity for 2nd order adjacencies
      :param typeReq: Type of adjacent entities being requested
      :return: If ``entities`` is a single element, an array of adjacent
               entities. Otherwise, a tuple containing an array of offsets and
               an array of adjacent entities.

   .. method:: getAdjEntIndices(entSet, typeRequestor, topoRequestor, typeRequested)

   .. method:: createEntSet(isList)

   .. method:: destroyEntSet(entSet)

   .. method:: setVtxCoords(entities, coords[, storageOrder])

   .. method:: createVtx(coords[, storageOrder])

   .. method:: createEnt(topo, entities)

   .. method:: createEntArr(topo, entitites)

   .. method:: deleteEnt(entities)

   .. method:: createTag(name, size, type)

   .. method:: destroyTag(tag, forced)

   .. method:: getTagHandle(name)

   .. method:: setData(entities, tag, data[, type])

   .. method:: getData(entities, tag[, type])

   .. method:: getAllTags(entities)

   .. method:: rmvTag(entities, tag)


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