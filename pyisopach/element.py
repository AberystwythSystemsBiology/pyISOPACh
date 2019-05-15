import json
import numpy as np
from .periodic_table import get_periodic_table


class Element(object):

    periodic_table = get_periodic_table()

    def __init__(self, symbol, count):
        self.symbol = symbol
        self.count = count
        self.atomic_weight = self._get_atomic_weight(symbol)
        self.isotopic_weights = self._get_isotopic_weight(symbol)
        self.atomic_charge = self._get_atomic_charge(symbol)
        self.ratios = self._get_ratios(symbol)

    @property
    def molecular_weight(self):
        return self.atomic_weight * float(self.count)

    @property
    def isotopic_ratios(self):
        return self.periodic_table[self.symbol]["isotopic_ratio"]

    @property
    def atomic_charge(self, symbol):
        return self.periodic_table[symbol]["atomic_charge"]

    @property
    def atomic_weight(self, symbol):
        isotopic_weight = self._get_isotopic_weight(self.symbol)
        ratios = self.periodic_table[symbol]["isotopic_ratio"]
        return float(np.matrix(ratios) * np.transpose(np.matrix(isotopic_weight)))

    @property
    def isotopic_weight(self, symbol):
        return self.periodic_table[self.symbol]["isotopic_weight"]
