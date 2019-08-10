import unittest
import pyisopach
import numpy as np

class TestElement(unittest.TestCase):

    def test_incorrect_element(self):
        # Test to see if nonsense element throws an exception
        self.assertRaises(KeyError, pyisopach.Element, "spooky scary skeletons", 69)

    def test_molecular_weight(self):
        # Check whether the calculation of molecular weight is correct.
        elem = pyisopach.Element("O", 1)
        self.assertAlmostEqual(elem.molecular_weight, 15.999, 3)

    def test_isotopic_ratios_sum_one(self):
        # Isotopic ratios should always equal 1.0, regardless of how many elements are
        # passed.
        elem = pyisopach.Element("O", 1)
        self.assertEquals(sum(elem.isotopic_ratios), 1.0)

        elem = pyisopach.Element("O", 2)
        self.assertEquals(sum(elem.isotopic_ratios), 1.0)

    def test_isotopic_ratios_values(self):
        # Isotopic Ratios Taken from IUPAC handbook.
        elem = pyisopach.Element("O", 1)
        self.assertTrue(np.array_equal(elem.isotopic_ratios, [0.99757, 0.00038, 0.00205]))

    def test_atomic_charge(self):
        # The atomic charge of Oxygen == -2
        elem = pyisopach.Element("O", 12)
        self.assertEquals(elem.atomic_charge, -2)

    def test_atomic_weight(self):
        elem = pyisopach.Element("O", 1)
        self.assertAlmostEqual(elem.atomic_weight, 15.99, 1)

    def test_isotopic_weight(self):
        elem = pyisopach.Element("O", 1)
        self.assertAlmostEqual(sum(elem.isotopic_weight), 50.99, 1)

if __name__ == '__main__':
    unittest.main()