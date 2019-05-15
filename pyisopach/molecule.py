import itertools
from collections import Counter
import operator
from re import findall
from .element import Element

ELECTRON_WEIGHT = 0.0005484


class Molecule:
    def __init__(self, molecular_formula):
        self.molecular_formula = molecular_formula

    @property
    def num_atoms(self):
        return sum(self.structure_dict.values())

    @property
    def accurate_mass(self):
        element_list = self.generate_element_list(self.structure_dict)
        return sum([x.calculate_weight() for x in element_list])

    @property
    def _structure_dict(self):
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

    def generate_element_list(self, chemical_structure=None):
        element_list = [
            Element(symbol, self.structure_dict[symbol])
            for symbol in self.structure_dict.keys()
        ]
        return element_list

    def _rule_dict_elements_list(self, rule_dict):
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

        return self.generate_element_list(_structure_dict)

    def isotopic_distribution(self, rule_dict=None, electrons=1, charge=0):

        if charge == 0:
            charge = 1

        if rule_dict == None:
            elements_list = self.generate_element_list()
        else:
            elements_list = self._rule_dict_elements_list(rule_dict)

        def cartesian(ratios,
                      weights,
                      f_ratios,
                      f_weights,
                      count=1,
                      threshold=0.05):
            n_ratios = []
            n_weights = []
            normalised_ratio = [n / max(f_ratios) for n in f_ratios]

            for i in enumerate(ratios[count]):
                r = ratios[count][i[0]]
                w = weights[count][i[0]]
                for j in enumerate(normalised_ratio):
                    current_n_ratio = normalised_ratio[j[0]] * 100
                    current_s_weight = f_weights[j[0]]
                    t_w = current_n_ratio * r
                    if t_w > threshold:
                        n_ratios += [current_n_ratio * r]
                        n_weights += [current_s_weight + w]
            count = count + 1
            if count < len(ratios) and len(n_ratios) < 1000:
                n_ratios, n_weights = cartesian(ratios, weights, n_ratios,
                                                n_weights, count)
            return n_ratios, n_weights

        def isotopes():
            weights = []
            ratios = []
            for element in elements_list:
                for i in range(element.count):
                    weights.append(element.isotopic_weights)
                    ratios.append(element.ratios)
            return cartesian(ratios, weights, ratios[0], weights[0])

        def calculate_distributions(gen_iso):
            signals = dict((key, tuple(v for (k, v) in pairs))
                           for (key, pairs) in itertools.groupby(
                               sorted(gen_iso), operator.itemgetter(0)))

            distributions = {}
            peak_intensity = max(signals.values())[0]

            for mass, intensity in signals.items():
                mass = round(
                    (mass + (ELECTRON_WEIGHT * electrons)) / abs(charge), 5)
                relative_intensity = round(
                    float(sum(intensity)) * 100 / peak_intensity, 5)
                if mass not in distributions.keys():
                    distributions[mass] = relative_intensity
                else:
                    distributions[mass] += relative_intensity
            return sorted(distributions.items(), key=lambda x: x[0])

        ratios, weights = isotopes()
        distributions = calculate_distributions(
            [(weights[index], ratio) for index, ratio in enumerate(ratios)
             if ratio > 1e-6])
        return distributions
