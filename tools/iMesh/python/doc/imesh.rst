=================
 iMesh Interface
=================

.. module:: itaps

.. class:: itaps.iMesh

   .. attribute:: rootSet

   .. attribute:: geometricDimension

   .. attribute:: dfltStorage

   .. attribute:: adjTable

   .. method:: load(set, filename[, options])

      Load a mesh

   .. method:: save(set, filename[, options])

      Save the mesh

   .. method:: getNumOfType(type)

   .. method:: getNumOfTopo(topo)

   .. method:: areEHValid(doReset)

   .. method:: getEntities(set, type, topo)

   .. method:: getVtxCoords(entities[, storageOrder])

   .. method:: getEntType(entities)

   .. method:: getEntTopo(entities)

   .. method:: getEntAdj(entities, typeReq)

   .. method:: getEnt2ndAdj(entities, bridgeType, typeReq)

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