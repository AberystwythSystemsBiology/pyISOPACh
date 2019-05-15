from .element import Element
from .molecule import Molecule
from .periodic_table import get_periodic_table

# Only load in once, and once only.
periodic_table = get_periodic_table()
