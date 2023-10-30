import unittest
import Hypercube
import SepVector
import numpy as np

class TestSepVector(unittest.TestCase):

    def testDiff(self):
        vec1 = SepVector.getSepVector(ns=[10, 10])
        vec1.rand()
        vec2 = vec1.clone()
        self.assertFalse(vec1.isDifferent(vec2))

        vec2.rand()
        self.assertTrue(vec1.isDifferent(vec2))

    def testWindow(self):
        vec1 = SepVector.getSepVector(ns=[30, 30], os=[0, 0], ds=[1., 1.])
        vec2 = vec1.window(f1=3, n2=10)
        axes = vec2.getHyper().axes
        self.assertEqual(axes[0].n, 27)
        self.assertEqual(axes[0].o, 3.)
        self.assertEqual(axes[1].n, 10)

    def testCreateFromNumpy(self):
        x=np.zeros((3,))
        x[0]=4
        x[1]=2
        x[2]=1
        y=SepVector.getSepVector(x)
        m=y.getNdArray()
        self.assertEqual(m[0], 4)
        self.assertEqual(m[1], 2)
        self.assertEqual(m[2], 1)


if __name__ == '__main__':
    unittest.main()
