from petsc4py import PETSc
from slepc4py import SLEPc
import unittest

# --------------------------------------------------------------------

class BaseTestObject(object):

    CLASS, FACTORY = None, 'create'
    TARGS, KARGS = (), {}
    BUILD = None
    def setUp(self):
        self.obj = self.CLASS()
        getattr(self.obj,self.FACTORY)(*self.TARGS, **self.KARGS)
        if not self.obj: self.obj.create()

    def tearDown(self):
        self.obj = None

    def testTypeRegistry(self):
        type_reg = PETSc.__type_registry__
        cookie = self.obj.getCookie()
        self.assertTrue(type_reg[cookie] is self.CLASS )

    def testLogClass(self):
        name = self.CLASS.__name__
        logcls = PETSc.Log.Class(name)
        cookie = self.obj.getCookie()
        self.assertEqual(logcls.id, cookie)

    def testClass(self):
        self.assertTrue(isinstance(self.obj, self.CLASS))
        self.assertTrue(type(self.obj) is self.CLASS)

    def testNonZero(self):
        self.assertTrue(bool(self.obj))

    def testDestroy(self):
        self.assertTrue(bool(self.obj))
        self.obj.destroy()
        self.assertFalse(bool(self.obj))
        ## self.assertRaises(PETSc.Error, self.obj.destroy)
        ## self.assertTrue(self.obj.this is this)

    def testOptions(self):
        self.assertFalse(self.obj.getOptionsPrefix())
        prefix1 = 'my_'
        self.obj.setOptionsPrefix(prefix1)
        self.assertEqual(self.obj.getOptionsPrefix(), prefix1)
        prefix2 = 'opt_'
        self.obj.setOptionsPrefix(prefix2)
        self.assertEqual(self.obj.getOptionsPrefix(), prefix2)
        ## self.obj.appendOptionsPrefix(prefix1)
        ## self.assertEqual(self.obj.getOptionsPrefix(),
        ##                  prefix2 + prefix1)
        ## self.obj.prependOptionsPrefix(prefix1)
        ## self.assertEqual(self.obj.getOptionsPrefix(),
        ##                  prefix1 + prefix2 + prefix1)
        self.obj.setFromOptions()

    def testName(self):
        oldname = self.obj.getName()
        newname = '%s-%s' %(oldname, oldname)
        self.obj.setName(newname)
        self.assertEqual(self.obj.getName(), newname)
        self.obj.setName(oldname)
        self.assertEqual(self.obj.getName(), oldname)

    def testRefCount(self):
        self.assertEqual(self.obj.getRefCount(), 1)
        self.obj.incRef()
        self.assertEqual(self.obj.getRefCount(), 2)
        self.obj.incRef()
        self.assertEqual(self.obj.getRefCount(), 3)
        self.obj.decRef()
        self.assertEqual(self.obj.getRefCount(), 2)
        self.obj.decRef()
        self.assertEqual(self.obj.getRefCount(), 1)
        self.obj.decRef()
        self.assertFalse(bool(self.obj))

    def testComposeQuery(self):
        self.assertEqual(self.obj.getRefCount(), 1)
        self.obj.compose('myobj', self.obj)
        self.assertTrue(type(self.obj.query('myobj')) is self.CLASS)
        self.assertEqual(self.obj.query('myobj'), self.obj)
        self.assertEqual(self.obj.getRefCount(), 2)
        self.obj.compose('myobj', None)
        self.assertEqual(self.obj.getRefCount(), 1)
        self.assertEqual(self.obj.query('myobj'), None)

    def testComm(self):
        comm = self.obj.getComm()
        self.assertTrue(isinstance(comm, PETSc.Comm))
        self.assertTrue(comm in [PETSc.COMM_SELF, PETSc.COMM_WORLD])

    def testProperties(self):
        self.assertEqual(self.obj.getCookie(),    self.obj.cookie)
        self.assertEqual(self.obj.getClassName(), self.obj.klass)
        self.assertEqual(self.obj.getType(),      self.obj.type)
        self.assertEqual(self.obj.getName(),      self.obj.name)
        self.assertEqual(self.obj.getComm(),      self.obj.comm)
        self.assertEqual(self.obj.getRefCount(),  self.obj.refcount)

# --------------------------------------------------------------------

class TestObjectEPS(BaseTestObject, unittest.TestCase):
    CLASS = SLEPc.EPS

class TestObjectSVD(BaseTestObject, unittest.TestCase):
    CLASS = SLEPc.SVD

class TestObjectST(BaseTestObject, unittest.TestCase):
    CLASS = SLEPc.ST

class TestObjectIP(BaseTestObject, unittest.TestCase):
    CLASS = SLEPc.IP

# --------------------------------------------------------------------

if __name__ == '__main__':
    unittest.main()
