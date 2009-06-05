from itaps import *
import unittest
import tempfile

class TestBasic(unittest.TestCase):
    def testMinimal(self):
        mesh = iMesh.Mesh("hello there")

        mesh.geometricDimension = 2
        self.assertEqual(mesh.geometricDimension, 2)
        root = mesh.rootSet
        mesh.dfltStorage

        self.assert_(mesh.areEHValid(True))
        self.assertEqual(root.getNumOfType(iBase.Type.all),     0)
        self.assertEqual(root.getNumOfTopo(iMesh.Topology.all), 0)
        self.assertEqual(mesh.adjTable.shape, (4,4))

    def testVertex(self):
        mesh = iMesh.Mesh()
        ent = mesh.createVtx([1,2,3])
        root = mesh.rootSet

        self.assertEqual(root.getNumOfType(iBase.Type.vertex),    1)
        self.assertEqual(root.getNumOfTopo(iMesh.Topology.point), 1)

        self.assert_( (mesh.getVtxCoords(ent) == [1,2,3]).all() )
        self.assert_( (mesh.getVtxCoords([ent], iBase.StorageOrder.interleaved)
                       == [[1,2,3]]).all() )

        self.assertEqual(mesh.getEntType(ent), iBase.Type.vertex)
        self.assertEqual(mesh.getEntTopo(ent), iMesh.Topology.point)

        mesh.setVtxCoords(ent,[4,5,6])
        self.assert_( (mesh.getVtxCoords(ent) == [4,5,6]).all() )

    def testVertices(self):
        mesh = iMesh.Mesh()
        verts = [[1,2,3], [4,5,6], [7,8,9], [10,11,12]]
        ents = mesh.createVtx(verts, iBase.StorageOrder.interleaved)
        root = mesh.rootSet

        self.assertEqual(root.getNumOfType(iBase.Type.vertex),    4)
        self.assertEqual(root.getNumOfTopo(iMesh.Topology.point), 4)

        coords = mesh.getVtxCoords(ents, iBase.StorageOrder.interleaved)
        self.assert_( (coords == verts).all())

        self.assertEqual(mesh.getEntType(ents[0]), iBase.Type.vertex)
        self.assert_( (mesh.getEntType(ents) == 4*[iBase.Type.vertex]).all() )

        self.assertEqual(mesh.getEntTopo(ents[0]), iMesh.Topology.point)
        self.assert_( (mesh.getEntTopo(ents) == 
                       4*[iMesh.Topology.point]).all() )

        verts = [[12,11,10], [9,8,7], [6,5,4], [3,2,1]]
        mesh.setVtxCoords(ents, verts, iBase.StorageOrder.interleaved)
        
        coords = mesh.getVtxCoords(ents, iBase.StorageOrder.interleaved)
        self.assert_( (coords == verts).all())

    def testCreateEntArr(self):
        mesh = iMesh.Mesh()
        verts = [[0,0,0], [0,0,1], [0,1,0], [0,1,1]]
        ents = mesh.createVtx(verts, iBase.StorageOrder.interleaved)
        root = mesh.rootSet
        topo = iMesh.Topology

        lines = mesh.createEntArr(topo.line_segment,ents)[0]
        self.assertEqual(root.getNumOfType(iBase.Type.vertex),  4)
        self.assertEqual(root.getNumOfType(iBase.Type.edge),    2)

        self.assertEqual(root.getNumOfTopo(topo.point),         4)
        self.assertEqual(root.getNumOfTopo(topo.line_segment),  2)

    def testSave(self):
        file = tempfile.NamedTemporaryFile()

        mesh = iMesh.Mesh()
        verts = [1,2,3]
        mesh.createVtx(verts)
        
        mesh.save(file.name)
        
        mesh = iMesh.Mesh()
        root = mesh.rootSet
        mesh.load(file.name)
        ents = root.getEntities(iBase.Type.all, iMesh.Topology.all)

        self.assertEqual(root.getNumOfType(iBase.Type.vertex),    1)
        self.assertEqual(root.getNumOfTopo(iMesh.Topology.point), 1)

        coords = mesh.getVtxCoords(ents, iBase.StorageOrder.interleaved)
        self.assert_( (coords == [1,2,3]).all() )

    def testAltSave(self):
        file = tempfile.NamedTemporaryFile()

        mesh = iMesh.Mesh()
        verts = [1,2,3]
        mesh.createVtx(verts)
        
        mesh.rootSet.save(file.name)
        
        mesh = iMesh.Mesh()
        root = mesh.rootSet
        root.load(file.name)
        ents = root.getEntities(iBase.Type.all, iMesh.Topology.all)

        self.assertEqual(root.getNumOfType(iBase.Type.vertex),    1)
        self.assertEqual(root.getNumOfTopo(iMesh.Topology.point), 1)

        coords = mesh.getVtxCoords(ents, iBase.StorageOrder.interleaved)
        self.assert_( (coords == [1,2,3]).all() )

if __name__ == '__main__':
    unittest.main()
