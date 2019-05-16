import itertools
from collections import Counter
import operator
from re import findall
from .element import Element
from typing import List, Tuple

ELECTRON_WEIGHT = 0.0005484

class Molecule:
    def __init__(self, molecular_formula: str):
        self.molecular_formula = molecular_formula
        self._structure_dict = self._generate_structure_dict()

    @property
    def num_atoms(self) -> int:
        return sum(self._structure_dict.values())

    @property
    def molecular_weight(self) -> float:
        return sum([elem.molecular_weight for elem in self._as_elements])

    @property
    def _as_elements(self) -> List[Element]:
        return ([Element(s, c) for (s, c) in self._structure_dict.items()])

    def _generate_structure_dict(self) -> dict:
         parsed = findall(r"([A-Z][a-z]*)(\d*)|(\()|(\))(\d*)", self.molecular_formula)
         structure_dict = {}
         for element_details in parsed:
             element = element_details[0]
             if element not in structure_dict:
                 structure_dict[element] = 0
             element_count = sum([int(x) for x in element_details[1:] if x != ""])
             if element_count > 0:
                 structure_dict[element] += element_count
             else:
                structure_dict[element] = 1
         return structure_dict


    @property
    def isotopic_distribution(self):

        def _get_weights_and_ratios() -> Tuple[int, int]:
            weights, ratios = [], []

            for elem in self._as_elements:
                for _ in range(elem.count):
                    weights.append(elem.isotopic_weight)
                    ratios.append(elem.isotopic_ratios)

            return weights, ratios

        def _cartesian_product(
                weights : list,
                ratios : list,
                calc_weights: list = [],
                calc_ratios: list = [],
                count: int = 1,
                threshold: float = 0.05):

            calc_ratios = []
            calc_weights = []

            if calc_ratios == []:
                calc_ratios = ratios[0]
                calc_weights = weights[0]

            _normalised_ratios = [n / max(calc_ratios) for n in calc_ratios]

            for ratio_indx in range(len(ratios[count])):
                _ratio = ratios[count][ratio_indx]
                _weight = weights[count][ratio_indx]

                for norm_ratio_index in range(len(_normalised_ratios)):
                    _norm_ratio = _normalised_ratios[norm_ratio_index] * 100
                    _norm_weight = calc_weights[norm_ratio_index]

                    _transformed_ratio = _norm_ratio * _ratio

                    if _transformed_ratio > threshold:
                        calc_ratios += [_transformed_ratio]
                        calc_weights += [_norm_weight + _weight]

            count += 1
            if count < len(ratios) and len(calc_ratios) < 10000:
                calc_weights, calc_ratios = _cartesian_product(weights, ratios, calc_weights, calc_ratios, count)
            return calc_weights, calc_ratios

        def _filter_cartesian_product(calc_weight, calc_ratios):
            pass

        def _generate_isotope_patterns():
            pass

        weights, ratios = _get_weights_and_ratios()
        calc_weights, calc_ratios = _cartesian_product(weights, ratios)
