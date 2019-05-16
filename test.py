import pyisopach
mol = pyisopach.Molecule("C12H22O11")
# Return molecular weight
mol.molecular_weight
342.2970125766493
# Calculate isotopic distribution/pattern
print(mol.isotopic_distribution())
