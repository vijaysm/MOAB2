from itaps import *
import unittest
from numpy import *

class TestTags(unittest.TestCase):
    def setUp(self):
        self.mesh = iMesh()

        self.itag = self.mesh.createTag('int',    1, 'i')
        self.dtag = self.mesh.createTag('double', 1, 'd')
        self.etag = self.mesh.createTag('handle', 1, 'E')
        self.btag = self.mesh.createTag('bytes',  3, 'b')

        self.ents = self.mesh.createVtx([[1,2,3], [4,5,6], [7,8,9]],
                                        iBase.storageOrder.interleaved)
        self.ent = self.ents[0]
        self.set = self.mesh.createEntSet(True)

    def testCreation(self):
        self.assertEqual(self.itag.name, 'int')
        self.assertEqual(self.itag.type, 'i')
        self.assertEqual(self.itag.sizeValues, 1)
        self.assertEqual(self.itag.sizeBytes, 4)

    def testFind(self):
        t = self.mesh.getTagHandle('int')
        self.assertEqual(t.name, self.itag.name)

        self.assertRaises(RuntimeError, self.mesh.getTagHandle, 'potato')

    def testIntData(self):
        self.mesh.setData(self.ent, self.itag, 42)
        self.assertEqual(self.mesh.getData(self.ent, self.itag), 42)
        self.assertEqual(self.mesh.getData(self.ent, self.itag, 'i'), 42)

        self.mesh.rmvTag(self.ent, self.itag)
        self.assertRaises(RuntimeError, self.mesh.getData,
                          self.ent,self.itag)

    def testDblData(self):
        self.mesh.setData(self.ent, self.dtag, 42.0)
        self.assertEqual(self.mesh.getData(self.ent, self.dtag), 42.0)
        self.assertEqual(self.mesh.getData(self.ent, self.dtag, 'd'), 42.0)

        self.mesh.rmvTag(self.ent, self.dtag)
        self.assertRaises(RuntimeError, self.mesh.getData,
                          self.ent,self.dtag)

    def testEHData(self):
        self.mesh.setData(self.ent, self.etag, self.ent)
        self.assertEqual(self.mesh.getData(self.ent, self.etag), self.ent)
        self.assertEqual(self.mesh.getData(self.ent, self.etag, 'E'), self.ent)

        self.mesh.rmvTag(self.ent, self.etag)
        self.assertRaises(RuntimeError, self.mesh.getData,
                          self.ent,self.etag)        

    def testRawData(self):
        data = array([1,2,3], int8)
        self.mesh.setData(self.ent, self.btag, data)
        self.assert_( (self.mesh.getData(self.ent, self.btag) == data).all() )
        self.assert_( (self.mesh.getData(self.ent, self.btag, 'b') == data)
                      .all() )

        self.mesh.rmvTag(self.ent, self.btag)
        self.assertRaises(RuntimeError, self.mesh.getData,
                          self.ent,self.btag)

    def testGetAll(self):
        self.mesh.setData(self.ent, self.itag, 42)
        self.mesh.setData(self.ent, self.dtag, 42)

        tags = self.mesh.getAllTags(self.ent)
        self.assertEqual(tags[0].name, self.itag.name) # TODO: ignore order?
        self.assertEqual(tags[1].name, self.dtag.name)


    def testIntArrData(self):
        self.mesh.setData(self.ents, self.itag, 3*[42])

        self.assert_((self.mesh.getData(self.ents, self.itag, 'i') ==
                      3*[42]).all())
        self.assertEqual(self.mesh.getData(self.ents[0], self.itag, 'i'), 42)

        self.mesh.rmvTag(self.ents, self.itag)
        self.assertRaises(RuntimeError, self.mesh.getData,
                          self.ents,self.itag)

    def testDblArrData(self):
        self.mesh.setData(self.ents, self.dtag, 3*[42.0])

        self.assert_((self.mesh.getData(self.ents, self.dtag, 'd') ==
                      3*[42.0]).all())
        self.assertEqual(self.mesh.getData(self.ents[0], self.dtag, 'd'), 42.0)

        self.mesh.rmvTag(self.ents, self.dtag)
        self.assertRaises(RuntimeError, self.mesh.getData,
                          self.ents,self.dtag)

    def testEHArrData(self):
        self.mesh.setData(self.ents, self.etag, self.ents)

        self.assertEqual(str(self.mesh.getData(self.ents, self.etag, 'E')),
                         str(self.ents))
        self.assertEqual(self.mesh.getData(self.ents[0], self.etag, 'E'),
                         self.ents[0])

        self.mesh.rmvTag(self.ents, self.etag)
        self.assertRaises(RuntimeError, self.mesh.getData,
                          self.ents,self.etag)

    def testRawArrData(self):
        data = array(3*[1,2,3], int8)
        self.mesh.setData(self.ents, self.btag, data)

        self.assert_((self.mesh.getData(self.ents, self.btag, 'b') ==
                      data).all())
        self.assert_((self.mesh.getData(self.ents[0], self.btag, 'b') == 
                      data[0:3]).all())

        self.mesh.rmvTag(self.ents, self.btag)
        self.assertRaises(RuntimeError, self.mesh.getData,
                          self.ents,self.btag)



    def testIntSetData(self):
        self.mesh.setData(self.set, self.itag, 42)
        self.assertEqual(self.mesh.getData(self.set, self.itag), 42)
        self.assertEqual(self.mesh.getData(self.set, self.itag, 'i'), 42)

        self.mesh.rmvTag(self.set, self.itag)
        self.assertRaises(RuntimeError, self.mesh.getData,
                          self.set,self.itag)

    def testDblSetData(self):
        self.mesh.setData(self.set, self.dtag, 42)
        self.assertEqual(self.mesh.getData(self.set, self.dtag), 42)
        self.assertEqual(self.mesh.getData(self.set, self.dtag, 'd'), 42)

        self.mesh.rmvTag(self.set, self.dtag)
        self.assertRaises(RuntimeError, self.mesh.getData,
                          self.set,self.dtag)
        
    def testEHSetData(self):
        self.mesh.setData(self.set, self.etag, self.ent)
        self.assertEqual(self.mesh.getData(self.set, self.etag), self.ent)
        self.assertEqual(self.mesh.getData(self.set, self.etag, 'd'), self.ent)

        self.mesh.rmvTag(self.set, self.etag)
        self.assertRaises(RuntimeError, self.mesh.getData,
                          self.set,self.etag)

    def testRawSetData(self):
        data = array([1,2,3], int8)
        self.mesh.setData(self.set, self.btag, data)

        self.assert_((self.mesh.getData(self.set, self.btag) == data).all())
        self.assert_((self.mesh.getData(self.set, self.btag, 'b') == data)
                     .all())

        self.mesh.rmvTag(self.set, self.btag)
        self.assertRaises(RuntimeError, self.mesh.getData,
                          self.set,self.btag)

    def testGetAllSet(self):
        self.mesh.setData(self.set, self.itag, 42)
        self.mesh.setData(self.set, self.dtag, 42)

        tags = self.mesh.getAllTags(self.set)
        self.assertEqual(tags[0].name, self.itag.name) # TODO: ignore order?
        self.assertEqual(tags[1].name, self.dtag.name)
        

suite = unittest.TestLoader().loadTestsFromTestCase(TestTags)
unittest.TextTestRunner(verbosity=2).run(suite)
