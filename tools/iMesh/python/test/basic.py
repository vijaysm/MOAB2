from itaps import *
import unittest
import tempfile

class TestBasic(unittest.TestCase):
    def testMinimal(self):
        mesh = iMesh("hello there")

        mesh.geometricDimension = 2
        self.assertEqual(mesh.geometricDimension, 2)
        root = mesh.rootSet
        mesh.dfltStorage

        self.assert_(mesh.areEHValid(True))
        self.assertEqual(root.getNumOfType(iBase.type.all),     0)
        self.assertEqual(root.getNumOfTopo(iMesh.topology.all), 0)
        self.assertEqual(mesh.adjTable.shape, (4,4))

    def testVertex(self):
        mesh = iMesh()
        ent = mesh.createVtx([1,2,3])
        root = mesh.rootSet

        self.assertEqual(root.getNumOfType(iBase.type.vertex),    1)
        self.assertEqual(root.getNumOfTopo(iMesh.topology.point), 1)

        self.assert_( (mesh.getVtxCoords(ent) == [1,2,3]).all() )
        self.assert_( (mesh.getVtxCoords([ent], iBase.storageOrder.interleaved)
                       == [[1,2,3]]).all() )

        self.assertEqual(mesh.getEntType(ent), iBase.type.vertex)
        self.assertEqual(mesh.getEntTopo(ent), iMesh.topology.point)

        mesh.setVtxCoords(ent,[4,5,6])
        self.assert_( (mesh.getVtxCoords(ent) == [4,5,6]).all() )

    def testVertices(self):
        mesh = iMesh()
        verts = [[1,2,3], [4,5,6], [7,8,9], [10,11,12]]
        ents = mesh.createVtx(verts, iBase.storageOrder.interleaved)
        root = mesh.rootSet

        self.assertEqual(root.getNumOfType(iBase.type.vertex),    4)
        self.assertEqual(root.getNumOfTopo(iMesh.topology.point), 4)

        coords = mesh.getVtxCoords(ents, iBase.storageOrder.interleaved)
        self.assert_( (coords == verts).all())

        self.assertEqual(mesh.getEntType(ents[0]), iBase.type.vertex)
        self.assert_( (mesh.getEntType(ents) == 4*[iBase.type.vertex]).all() )

        self.assertEqual(mesh.getEntTopo(ents[0]), iMesh.topology.point)
        self.assert_( (mesh.getEntTopo(ents) == 
                       4*[iMesh.topology.point]).all() )

        verts = [[12,11,10], [9,8,7], [6,5,4], [3,2,1]]
        mesh.setVtxCoords(ents, verts, iBase.storageOrder.interleaved)
        
        coords = mesh.getVtxCoords(ents, iBase.storageOrder.interleaved)
        self.assert_( (coords == verts).all())

    def testCreateEntArr(self):
        mesh = iMesh()
        verts = [[0,0,0], [0,0,1], [0,1,0], [0,1,1]]
        ents = mesh.createVtx(verts, iBase.storageOrder.interleaved)
        root = mesh.rootSet
        topo = iMesh.topology

        lines = mesh.createEntArr(topo.line_segment,ents)[0]
        self.assertEqual(root.getNumOfType(iBase.type.vertex),  4)
        self.assertEqual(root.getNumOfType(iBase.type.edge),    2)

        self.assertEqual(root.getNumOfTopo(topo.point),         4)
        self.assertEqual(root.getNumOfTopo(topo.line_segment),  2)

    def testSave(self):
        file = tempfile.NamedTemporaryFile()

        mesh = iMesh()
        verts = [1,2,3]
        mesh.createVtx(verts)
        
        mesh.save(file.name)
        
        mesh = iMesh()
        root = mesh.rootSet
        mesh.load(file.name)
        ents = root.getEntities(iBase.type.all, iMesh.topology.all)

        self.assertEqual(root.getNumOfType(iBase.type.vertex),    1)
        self.assertEqual(root.getNumOfTopo(iMesh.topology.point), 1)

        coords = mesh.getVtxCoords(ents, iBase.storageOrder.interleaved)
        self.assert_( (coords == [1,2,3]).all() )

    def testAltSave(self):
        file = tempfile.NamedTemporaryFile()

        mesh = iMesh()
        verts = [1,2,3]
        mesh.createVtx(verts)
        
        mesh.rootSet.save(file.name)
        
        mesh = iMesh()
        root = mesh.rootSet
        root.load(file.name)
        ents = root.getEntities(iBase.type.all, iMesh.topology.all)

        self.assertEqual(root.getNumOfType(iBase.type.vertex),    1)
        self.assertEqual(root.getNumOfTopo(iMesh.topology.point), 1)

        coords = mesh.getVtxCoords(ents, iBase.storageOrder.interleaved)
        self.assert_( (coords == [1,2,3]).all() )


suite = unittest.TestLoader().loadTestsFromTestCase(TestBasic)
unittest.TextTestRunner(verbosity=2).run(suite)
