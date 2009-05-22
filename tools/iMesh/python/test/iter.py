from itaps import *
import unittest

class TestIter(unittest.TestCase):
    def setUp(self):
        self.mesh = iMesh()
        self.empty = self.mesh.createEntSet(True)

        self.ent = self.mesh.createVtx([1,2,3])
        self.set = self.mesh.createEntSet(True)
        self.set.add(self.ent)

    def testEmpty(self):
        for i in self.empty.iterate(iBase.type.all, iMesh.topology.all):
            self.fail('empty iterator has >0 elements')

    def testArrEmpty(self):
        for i in self.empty.iterate(iBase.type.all, iMesh.topology.all, 16):
            self.fail('empty iterator has >0 elements')

    def testSingle(self):
        count = 0
        for i in self.set.iterate(iBase.type.all, iMesh.topology.all):
            count += 1
            self.assertEqual(i, self.ent)
        self.assertEqual(count, 1)

    def testArr(self):
        count = 0
        for i in self.set.iterate(iBase.type.all, iMesh.topology.all, 16):
            count += 1
            self.assertEqual(i[0], self.ent)
        self.assertEqual(count,1)

    def testAlternate(self):
        count = 0
        iter = iMesh.iterator(self.mesh, self.set, iBase.type.all,
                              iMesh.topology.all)
        for i in iter:
            count += 1
            self.assertEqual(i, self.ent)
        self.assertEqual(count, 1)
        iter.reset()


suite = unittest.TestLoader().loadTestsFromTestCase(TestIter)
unittest.TextTestRunner(verbosity=2).run(suite)
