import numpy as np

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
        parsed = findall(r"([A-Z][a-z]*)(\d*)|(\()|(\))(\d*)",
                         self.molecular_formula)
        structure_dict = {}
        for element_details in parsed:
            element = element_details[0]
            if element not in structure_dict:
                structure_dict[element] = 0
            element_count = sum(
                [int(x) for x in element_details[1:] if x != ""])
            if element_count > 0:
                structure_dict[element] += element_count
            else:
                structure_dict[element] = 1
        return structure_dict

    def isotopic_distribution(self, electrons: int = 1, charge: int = 0):
        def _get_weights_and_ratios() -> Tuple[int, int]:
            weights, ratios = [], []

            for elem in self._as_elements:
                for _ in range(elem.count):
                    weights.append(elem.isotopic_weight)
                    ratios.append(elem.isotopic_ratios)

            return weights, ratios

        def _cartesian_product(weights: list,
                               ratios: list,
                               f_weights: list = [],
                               f_ratios: list = [],
                               count: int = 1,
                               cartesian_threshold: float = 0.05
                               ) -> Tuple[np.array, np.array]:

            new_ratios = []
            new_weights = []

            if count == 1:
                f_ratios = ratios[0]
                f_weights = weights[0]

            normalised_ratio = [r / max(f_ratios) for r in f_ratios]

            for ratio_indx, _ in enumerate(ratios[count]):
                _ratio = ratios[count][ratio_indx]
                _weight = weights[count][ratio_indx]
                for norm_ratio_indx, _ in enumerate(normalised_ratio):
                    _norm_ratio = normalised_ratio[norm_ratio_indx] * 100  #
                    _norm_weight = f_weights[norm_ratio_indx]  #
                    _transformed_ratio = _norm_ratio * _ratio
                    if _transformed_ratio > cartesian_threshold:
                        new_ratios += [_transformed_ratio]
                        new_weights += [_norm_weight + _weight]
            count = count + 1
            if count < len(ratios) and len(new_ratios) < 10000:
                new_weights, new_ratios = _cartesian_product(
                    weights, ratios, new_weights, new_ratios, count)
            return np.array(new_weights), np.array(new_ratios)

        def _filter_low_ratios(
                calc_weights: np.array,
                calc_ratios: np.array,
                weight_limit: float = 1e-60) -> Tuple[np.array, np.array]:

            indx = calc_ratios > weight_limit
            return calc_weights[indx], calc_ratios[indx]

        def _generate_output(calc_weights: np.array, calc_ratios: np.array
                             ) -> Tuple[List[float], List[float]]:

            adj_weights = calc_weights + (ELECTRON_WEIGHT * electrons) / abs(charge)
            calc_dict = {x: 0 for x in np.unique(adj_weights)}


            for weight in calc_dict:
                calc_dict[weight] = np.sum(calc_ratios[adj_weights == weight]) * 100 / np.max(calc_ratios)

            return list(calc_dict.keys()), list(calc_dict.values())


        if charge == 0:
            charge = 1

        weights, ratios = _get_weights_and_ratios()
        calc_weights, calc_ratios = _cartesian_product(weights, ratios)
        calc_weights, calc_ratios = _filter_low_ratios(calc_weights,
                                                       calc_ratios)

        masses, intensities = _generate_output(calc_weights, calc_ratios)

        return np.array(masses), np.array(intensities)
