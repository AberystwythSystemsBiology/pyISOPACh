# pyISOPACh - ISOtope PAttern Calculator

This is a sister package of the DIMEdb project. This program calculates the isotopic distribution/pattern of a given chemical species.

## Example Usage

```python
# Import the package into python
>>> import pyisopach
# Create Molecule object for Sucrose
>>> mol = pyisopach.Molecule("C12H22O11")
# Return molecular weight
>>> mol.molecular_weight
342.2970125766493
# Calculate isotopic distribution/pattern
>>> mol.isotopic_distribution()
(array([342.11621155, 343.11956639, 344.12045733]), array([100.        ,  12.97887395,   2.260493  ]))
```


## License
Code released under [the MIT license](https://github.com/AberystwythSystemsBiology/pyisopach/blob/master/LICENSE).
