from itaps import *
import unittest

topo = iMesh.topology # shorthand

class TestAdj(unittest.TestCase):
    def setUp(self):
        self.mesh = iMesh()
        self.verts = [[0,0,0], [0,0,1], [0,1,0], [0,1,1]]
        self.ents = self.mesh.createVtx(self.verts,
                                        iBase.storageOrder.interleaved)

        self.lines = [ 
            self.mesh.createEnt(topo.line_segment, self.ents[0:2] )[0],
            self.mesh.createEnt(topo.line_segment, self.ents[1:3] )[0],
            self.mesh.createEnt(topo.line_segment, self.ents[2:4] )[0],
            self.mesh.createEnt(topo.line_segment, self.ents[::-3])[0] ]

    def tearDown(self):
        self.mesh  = None
        self.verts = None
        self.ents  = None
        self.root  = None
        self.lines = None

    @staticmethod
    def jaggify(arr,offsets):
        prev = offsets[0]
        jag = []
        for i in offsets[1:]:
            jag.append(arr[prev:i].tolist())
            prev = i
        return jag

    def testSquare(self):
        quad = self.mesh.createEnt(topo.quadrilateral, self.lines)[0]
        root = self.mesh.rootSet

        self.assertEqual(root.getNumOfType(iBase.type.vertex),  4)
        self.assertEqual(root.getNumOfType(iBase.type.edge),    4)
        self.assertEqual(root.getNumOfType(iBase.type.face),    1)

        self.assertEqual(root.getNumOfTopo(topo.point),         4)
        self.assertEqual(root.getNumOfTopo(topo.line_segment),  4)
        self.assertEqual(root.getNumOfTopo(topo.quadrilateral), 1)

        self.mesh.deleteEnt(quad)
        self.assertEqual(root.getNumOfType(iBase.type.face),    0)
        self.assertEqual(root.getNumOfTopo(topo.quadrilateral), 0)

        self.mesh.deleteEnt(self.lines)
        self.assertEqual(root.getNumOfType(iBase.type.edge),    0)
        self.assertEqual(root.getNumOfTopo(topo.line_segment),  0)

        self.mesh.deleteEnt(self.ents)
        self.assertEqual(root.getNumOfType(iBase.type.vertex),  0)
        self.assertEqual(root.getNumOfTopo(topo.point),         0)

    def testAdj(self):
        adj = self.mesh.getEntAdj(self.ents[1], iBase.type.all)
        self.assertEqual(adj.tolist(), self.lines[0:2])

        adj = self.mesh.getEntAdj(self.ents, iBase.type.all)
        self.assertEqual(TestAdj.jaggify(adj[0],adj[1]),
                         [ self.lines[::3],
                           self.lines[0:2],
                           self.lines[1:3],
                           self.lines[2:4] ])

    def test2ndAdj(self):
        quad = self.mesh.createEnt(topo.quadrilateral, self.lines)[0]

        adj = self.mesh.getEnt2ndAdj(self.ents[1], iBase.type.edge,
                                     iBase.type.vertex)
        self.assertEqual(adj.tolist(), self.ents[0:3:2].tolist())

        adj = self.mesh.getEnt2ndAdj(self.ents, iBase.type.edge,
                                     iBase.type.vertex)
        self.assertEqual(TestAdj.jaggify(adj[0],adj[1]),
                         [ self.ents[1::2].tolist(),
                           self.ents[0::2].tolist(),
                           self.ents[1::2].tolist(),
                           self.ents[0::2].tolist() ])

    def testAdjIndices(self):
        set = self.mesh.createEntSet(True)
        set.add(self.ents)
        adj = self.mesh.getAdjEntIndices(set, iBase.type.all, topo.all,
                                         iBase.type.all)

        self.assertEqual(adj[0].tolist(), self.ents.tolist())
        self.assertEqual(adj[1].tolist(), self.lines)
        self.assertEqual(TestAdj.jaggify(adj[2],adj[3]),
                         [ [0,3], [0,1], [1,2], [2,3] ])


suite = unittest.TestLoader().loadTestsFromTestCase(TestAdj)
unittest.TextTestRunner(verbosity=2).run(suite)
