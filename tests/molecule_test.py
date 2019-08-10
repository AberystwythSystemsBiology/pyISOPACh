import unittest
import numpy as np
import pyisopach

class MoleculeTest(unittest.TestCase):

    def check_molecule_dict_simple(self):
        mol = pyisopach.Molecule("H1")
        print(mol._generate_structure_dict())
        self.assertEquals(1,1)


if __name__ == '__main__':
    unittest.main()