import numpy as np
from .periodic_table import get_periodic_table
from typing import List

class Element(object):
    def __init__(self, symbol: str, count: int):
        self.symbol = symbol
        self.count = count
        self._periodic_table = get_periodic_table()[symbol]

    @property
    def molecular_weight(self) -> float:
        return self.atomic_weight * float(self.count)

    @property
    def isotopic_ratios(self) -> List[float]:
        return self._periodic_table["isotopic_ratio"]

    @property
    def atomic_charge(self) -> int:
        return self._periodic_table["atomic_charge"]

    @property
    def atomic_weight(self) -> float:
        isotopic_weight = self.isotopic_weight
        ratios = self._periodic_table["isotopic_ratio"]
        return float(np.matrix(ratios) * np.transpose(np.matrix(isotopic_weight)))

    @property
    def isotopic_weight(self) -> List[float]:
        return self._periodic_table["isotopic_weight"]
