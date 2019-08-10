import unittest
import pyisopach

class TestMolecule(unittest.TestCase):

    def test_structure_dict_easy(self):
        # Just create a molecule of two hydrogen, and ensure
        # that the number of hydrogen is equal to two.
        molecule = pyisopach.Molecule("H2")
        self.assertEqual(molecule._structure_dict["H"], 2)

    def test_structure_dict_intermediate(self):
        # Create a water molecule, and check if the number
        # of oxygen is equal to one.
        molecule = pyisopach.Molecule("H2O")
        self.assertEqual(molecule._structure_dict["O"], 1)

    def test_structure_dict_hard(self):
        # This is unlikely to ever happen, but given a weirdly
        # formulated chemical formula, the class should still
        # be able to determine quantities

        molecule = pyisopach.Molecule("H2OH67O419Ch")
        self.assertEqual(molecule._structure_dict["H"], 69)
        self.assertEqual(molecule._structure_dict["O"], 420)
        self.assertEqual(molecule._structure_dict["Ch"], 1)

    def test_num_atoms(self):
        molecule = pyisopach.Molecule("H2O")
        # Two hydrogen plus one oxygen
        self.assertEqual(molecule.num_atoms, 2+1)

    def test_molecular_weight(self):
        # This is Sucrose.
        molecule = pyisopach.Molecule("C12H22O11")

        # According to wikipedia, it has a molecular weight of
        # 342.29 g/mol
        self.assertAlmostEqual(molecule.molecular_weight, 342.29, 1)

    def test_isotopic_distribution(self):
        # Again, this is Sucrose.
        molecule = pyisopach.Molecule("C12H22O11")

        # Using an expensive piece of software that I borrowed,
        # I calculate the following isotopic pattern
        expected_masses = [342.11, 343.12, 344.12]
        expected_intensities = [100.00, 12.97, 2.26]

        for index, (mass, inten) in enumerate(zip(*molecule.isotopic_distribution())):
            self.assertAlmostEqual(mass, expected_masses[index], 1)
            self.assertAlmostEqual(inten, expected_intensities[index], 1)

if __name__ == '__main__':
    unittest.main()