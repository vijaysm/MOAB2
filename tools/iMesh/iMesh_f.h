#define iMesh_Instance integer
#define iMesh_EntityIterator integer
#define iMesh_EntityArrIterator integer
#define iBase_EntityHandle integer
#define iBase_EntitySetHandle integer
#define iBase_TagHandle integer



      parameter (iMesh_POINT = 0)
      parameter (iMesh_LINE_SEGMENT = 1)
      parameter (iMesh_POLYGON = 2)
      parameter (iMesh_TRIANGLE = 3)
      parameter (iMesh_QUADRILATERAL = 4)
      parameter (iMesh_POLYHEDRON = 5)
      parameter (iMesh_TETRAHEDRON = 6)
      parameter (iMesh_HEXAHEDRON = 7)
      parameter (iMesh_PRISM = 8)
      parameter (iMesh_PYRAMID = 9)
      parameter (iMesh_SEPTAHEDRON = 10)
      parameter (iMesh_ALL_TOPOLOGIES = 11)

      parameter (iMesh_UNAVAILABLE = 0)
      parameter (iMesh_IMMEDIATE = 1)
      parameter (iMesh_LOCAL_TRAVERSAL = 2)
      parameter (iMesh_GLOBAL_TRAVERSAL = 3)
