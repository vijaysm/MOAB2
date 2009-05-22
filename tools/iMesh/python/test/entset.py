from itaps import *
import unittest

class TestEntSet(unittest.TestCase):
    def setUp(self):
        self.mesh = iMesh()
        self.set = self.mesh.createEntSet(True)

    def tearDown(self):
        self.mesh = None

    def testCreation(self):
        self.assertEqual(self.set.isList, True)
        self.assertEqual(self.set.getNumEntSets(1), 0)
        self.assertEqual(self.mesh.rootSet.getNumEntSets(0), 1)

        foo = self.mesh.rootSet.getEntSets(0)[0]
        self.assertEqual(self.set, foo)

        self.mesh.destroyEntSet(self.set)
        self.assertEqual(self.mesh.rootSet.getNumEntSets(0), 0)

    def testEnt(self):
        ent = self.mesh.createVtx([1,2,3])

        self.assertFalse(self.set.contains(ent))

        self.set.add(ent)
        self.assert_(self.set.contains(ent))

        self.set.remove(ent)
        self.assertFalse(self.set.contains(ent))

    def testEntArr(self):
        ents = self.mesh.createVtx([[1,2,3], [4,5,6], [7,8,9]],
                                   iBase.storageOrder.interleaved)
        
        self.assertFalse(self.set.contains(ents).any())

        self.set.add(ents)
        self.assert_(self.set.contains(ents).all())

        self.set.remove(ents)
        self.assertFalse(self.set.contains(ents).any())

    def testEntSet(self):
        sub = self.mesh.createEntSet(True)

        self.assertFalse(self.set.contains(sub))

        self.set.add(sub)
        self.assert_(self.set.contains(sub))

        self.set.remove(sub)
        self.assertFalse(self.set.contains(sub))

    def testChildren(self):
        sub = self.mesh.createEntSet(True)
        self.set.addChild(sub)

        self.assert_(self.set.isChild(sub))
        self.assertEqual(self.set.getNumChildren(1), 1)
        self.assertEqual(sub.getNumParents(1),  1)

        self.assertEqual(self.set.getChildren(1)[0], sub)
        self.assertEqual(sub.getParents(1)[0], self.set)

        self.set.removeChild(sub)

        self.assert_(not self.set.isChild(sub))
        self.assertEqual(self.set.getNumChildren(1), 0)
        self.assertEqual(sub.getNumParents(1),  0)

    def testSubtract(self):
        set2 = self.mesh.createEntSet(True)
        ents = self.mesh.createVtx([[1,2,3], [4,5,6]],
                                   iBase.storageOrder.interleaved)
        self.set.add(ents)
        set2.add(ents[0])

        diff = self.set - set2
        self.assertFalse(diff.contains(ents[0]))
        self.assertTrue (diff.contains(ents[1]))

        diff = self.set.difference(set2)
        self.assertFalse(diff.contains(ents[0]))
        self.assertTrue (diff.contains(ents[1]))

    def testIntersect(self):
        set2 = self.mesh.createEntSet(True)
        ents = self.mesh.createVtx([[1,2,3], [4,5,6], [7,8,9]],
                                   iBase.storageOrder.interleaved)
        self.set.add(ents[0:2])
        set2.add(ents[1:3])

        sect = self.set & set2
        self.assertFalse(sect.contains(ents[0]))
        self.assertTrue (sect.contains(ents[1]))
        self.assertFalse(sect.contains(ents[2]))

        sect = self.set.intersection(set2)
        self.assertFalse(sect.contains(ents[0]))
        self.assertTrue (sect.contains(ents[1]))
        self.assertFalse(sect.contains(ents[2]))

    def testUnite(self):
        set2 = self.mesh.createEntSet(True)
        ents = self.mesh.createVtx([[1,2,3], [4,5,6]],
                                   iBase.storageOrder.interleaved)
        self.set.add(ents[0])
        set2.add(ents[1])

        union = self.set | set2
        self.assertTrue(union.contains(ents[0]))
        self.assertTrue(union.contains(ents[1]))

        union = self.set.union(set2)
        self.assertTrue(union.contains(ents[0]))
        self.assertTrue(union.contains(ents[1]))


suite = unittest.TestLoader().loadTestsFromTestCase(TestEntSet)
unittest.TextTestRunner(verbosity=2).run(suite)

