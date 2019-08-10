# pyISOPACh - a "fairly fast" ISOtope PAttern Calculator for Python

![PyPI - License](https://img.shields.io/pypi/l/pyISOPACh)
![PyPI](https://img.shields.io/pypi/v/pyisopach)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/pyisopach)
![PyPI - Status](https://img.shields.io/pypi/status/pyisopach)

This is a sister package of the DIMEdb project. This program calculates the isotopic distribution/pattern of a given chemical species.

## Installation

pyISOPACh requires Python 3+ and is unfortunately not compatible with Python 2. If you are still using Python 2, a clever workaround is to install Python 3 and use that instead.

The easiest way of installing pyISOPACh is using ```pip```:

```
pip install pyisopach
```

Alternatively, you can use ```git``` and ```pip``` in unison to get the development branch:

```
pip install https://github.com/AberystwythSystemsBiology/pyISOPACh
```

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
Code released under [the MIT license](https://github.com/AberystwythSystemsBiology/pyISOPACh/blob/master/LICENSE).
