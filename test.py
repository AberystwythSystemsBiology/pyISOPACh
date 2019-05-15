from pyisopach import Molecule

mol = Molecule("C12H22O11")


print(mol.molecular_weight)

print(mol.isotopic_distribution())
