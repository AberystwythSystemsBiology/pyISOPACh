import itertools
from collections import Counter
import operator
from re import findall
from .element import Element

ELECTRON_WEIGHT = 0.0005484

class Molecule:
    def __init__(self, molecular_formula):
        self.molecular_formula = molecular_formula
        self._structure_dict = self._generate_structure_dict()

    @property
    def num_atoms(self):
        return sum(self._structure_dict.values())

    @property
    def molecular_weight(self):
        return sum([elem.molecular_weight for elem in self._as_elements])

    @property
    def _as_elements(self):
        return ([Element(s, c) for (s, c) in self._structure_dict.items()])

    def _generate_structure_dict(self):
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


    def _apply_rule_dict(self, rule_dict):
        _structure_dict = self.structure_dict

        if "multiply" in rule_dict.keys():
            for element in _structure_dict.keys():
                _structure_dict[element] = _structure_dict[
                    element] * rule_dict["multiply"]

        elif "divide" in rule_dict.keys():
            for element in _structure_dict.keys():
                _structure_dict[element] = _structure_dict[
                    element] / rule_dict["divide"]

        for element in rule_dict["add"]:
            try:
                _structure_dict[element] += rule_dict["add"][element]
            except KeyError:
                _structure_dict[element] = rule_dict["add"][element]

        for element in rule_dict["remove"]:
            try:
                _structure_dict[element] -= rule_dict["remove"][element]
            except KeyError:
                continue

        return self._element_list
