from itaps import *
import unittest

topo = iMesh.Topology # shorthand

class TestAdj(unittest.TestCase):
    def setUp(self):
        self.mesh = iMesh.Mesh()
        self.verts = [[0,0,0], [0,0,1], [0,1,0], [0,1,1]]
        self.ents = self.mesh.createVtx(self.verts,
                                        iBase.StorageOrder.interleaved)

        self.lines = [ 
            self.mesh.createEnt(topo.line_segment, self.ents[0:2] )[0],
            self.mesh.createEnt(topo.line_segment, self.ents[1:3] )[0],
            self.mesh.createEnt(topo.line_segment, self.ents[2:4] )[0],
            self.mesh.createEnt(topo.line_segment, self.ents[::-3])[0] ]

    def testSquare(self):
        quad = self.mesh.createEnt(topo.quadrilateral, self.lines)[0]

        self.assertEqual(self.mesh.getNumOfType(iBase.Type.vertex),  4)
        self.assertEqual(self.mesh.getNumOfType(iBase.Type.edge),    4)
        self.assertEqual(self.mesh.getNumOfType(iBase.Type.face),    1)

        self.assertEqual(self.mesh.getNumOfTopo(topo.point),         4)
        self.assertEqual(self.mesh.getNumOfTopo(topo.line_segment),  4)
        self.assertEqual(self.mesh.getNumOfTopo(topo.quadrilateral), 1)

        self.mesh.deleteEnt(quad)
        self.assertEqual(self.mesh.getNumOfType(iBase.Type.face),    0)
        self.assertEqual(self.mesh.getNumOfTopo(topo.quadrilateral), 0)

        self.mesh.deleteEnt(self.lines)
        self.assertEqual(self.mesh.getNumOfType(iBase.Type.edge),    0)
        self.assertEqual(self.mesh.getNumOfTopo(topo.line_segment),  0)

        self.mesh.deleteEnt(self.ents)
        self.assertEqual(self.mesh.getNumOfType(iBase.Type.vertex),  0)
        self.assertEqual(self.mesh.getNumOfTopo(topo.point),         0)

    def testAdj(self):
        adj = self.mesh.getEntAdj(self.ents[1], iBase.Type.all)
        self.assertEqual(adj.tolist(), self.lines[0:2])

        adj = self.mesh.getEntAdj(self.ents, iBase.Type.all)

        for i in range( len(adj) ):
            self.assertEqual(adj.length(i), 2)

        self.assertEqual(adj[0].tolist(), self.lines[::3])
        self.assertEqual(adj[1].tolist(), self.lines[0:2])
        self.assertEqual(adj[2].tolist(), self.lines[1:3])
        self.assertEqual(adj[3].tolist(), self.lines[2:4])

        self.assertEqual(adj[0,0], self.lines[0])
        self.assertEqual(adj[2,1], self.lines[2])

        self.assertRaises(IndexError, adj.__getitem__, (0,2))

    def test2ndAdj(self):
        quad = self.mesh.createEnt(topo.quadrilateral, self.lines)[0]

        adj = self.mesh.getEnt2ndAdj(self.ents[1], iBase.Type.edge,
                                     iBase.Type.vertex)
        self.assertEqual(adj.tolist(), self.ents[0:3:2].tolist())

        adj = self.mesh.getEnt2ndAdj(self.ents, iBase.Type.edge,
                                     iBase.Type.vertex)

        for i in range( len(adj) ):
            self.assertEqual(adj.length(i), 2)

        self.assertEqual(adj[0].tolist(), self.ents[1::2].tolist())
        self.assertEqual(adj[1].tolist(), self.ents[0::2].tolist())
        self.assertEqual(adj[2].tolist(), self.ents[1::2].tolist())
        self.assertEqual(adj[3].tolist(), self.ents[0::2].tolist())

        self.assertEqual(adj[0,0], self.ents[1])
        self.assertEqual(adj[2,1], self.ents[3])

        self.assertRaises(IndexError, adj.__getitem__, (0,2))

    def testAdjIndices(self):
        set = self.mesh.createEntSet(True)
        set.add(self.ents)
        adj = set.getAdjEntIndices(iBase.Type.all, topo.all,
                                   iBase.Type.all)

        self.assertEqual(adj.entities.tolist(), self.ents.tolist())
        self.assertEqual(adj.adj.tolist(), self.lines)

        self.assert_( (adj.index(0) == [0,3]).all() )
        self.assert_( (adj.index(1) == [0,1]).all() )
        self.assert_( (adj.index(2) == [1,2]).all() )
        self.assert_( (adj.index(3) == [2,3]).all() )

if __name__ == '__main__':
    unittest.main()

