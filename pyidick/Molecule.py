#!/usr/bin/env python

from Element import Element

import itertools, operator

from re import findall
from rdkit.Chem import rdMolDescriptors, MolFromSmiles

class Molecule(object):
    ''' Molecule class
    
    Attributes:
        smiles: A string containing SMILES notation for a given molecule.
        __mol: A rdkit.mol object, generated from SMILES.
        molecular_formula: The molecular formula that shows the total number and kinds of atoms in a molecule.
        structure_dict: A structure dict, comprised from the generated_molecula formula attribute.
    '''

    def __init__(self, smiles):
        '''Inits Molecule will smiles, and generates all other attributes.'''
        self.smiles = smiles
        self.__mol = MolFromSmiles(smiles)
        self.molecular_formula = rdMolDescriptors.CalcMolFormula(self.__mol)
        self.structure_dict = self.__split()

    def generate_element_list(self, chemical_structure=None):
        if chemical_structure == None:
            chemical_structure = self.structure_dict
        element_list = [Element(symbol, chemical_structure[symbol]) for symbol in chemical_structure.keys()]
        return element_list

    def accurate_mass(self, chemical_structure=None):
        if chemical_structure == None:
            chemical_structure = self.__split()
        element_list = self.generate_element_list(chemical_structure)
        return sum([x.calculate_weight() for x in element_list])

    def __split(self):
        structure_dict = {}
        symbols = findall("[A-Z][a-z]*", self.molecular_formula)
        symbol_count = [findall("[0-9]+", e) for e in findall("[A-Z][a-z]*[0-9]*", self.molecular_formula)]
        for index, symbol in enumerate(symbols):
            count = symbol_count[index]
            if len(count) > 0:
                count = int(count[0])
            else:
                count = 1
            structure_dict[symbol] = count
        return structure_dict

    def _rule_dict_elements_list(self, rule_dict):
        _structure_dict = self.structure_dict

        if "multiply" in rule_dict.keys():
            for element in _structure_dict.keys():
                _structure_dict[element] = _structure_dict[element] * rule_dict["multiply"]

        elif "divide" in rule_dict.keys():
            for element in _structure_dict.keys():
                _structure_dict[element] = _structure_dict[element] / rule_dict["divide"]

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
        electron_weight = 0.0005484


        if charge == 0:
            charge = 1

        if rule_dict == None:
            elements_list = self.generate_element_list()
        else:
            elements_list = self._rule_dict_elements_list(rule_dict)

        def cartesian(ratios, weights, f_ratios, f_weights, count=1, threshold=0.05):
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
                n_ratios, n_weights = cartesian(ratios, weights, n_ratios, n_weights, count)
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
            signals = dict((key, tuple(v for (k,v) in pairs))
                          for (key, pairs) in itertools.groupby(sorted(gen_iso), operator.itemgetter(0)))

            distributions = {}
            peak_intensity = max(signals.values())[0]

            for mass, intensity in signals.items():
                mass = round((mass + (electron_weight * electrons)) / abs(charge), 5)
                relative_intensity = round(float(sum(intensity)) * 100 / peak_intensity, 5)
                if mass not in distributions.keys():
                    distributions[mass] = relative_intensity
                else:
                    distributions[mass] += relative_intensity
            return sorted(distributions.items(), key=lambda x: x[0])

        ratios, weights = isotopes()
        distributions = calculate_distributions([(weights[index], ratio) for index, ratio in enumerate(ratios) if ratio > 1e-6])
        # TODO: This is just a temporary solution.
        self.__init__(self.smiles)
        return distributions

if __name__ == "__main__":
    mol = Molecule("OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O")


    print "Sucrose"
    print "Molecular Formula: ", mol.molecular_formula


    base = mol.isotopic_distribution()[0][0]

    test_one = mol.isotopic_distribution({"add" : {}, "remove" : {"H" : 1}}, electrons=1, charge=-1)[0][0]

    print "[M-H]1-", test_one, base, test_one - base
    test_two = mol.isotopic_distribution(charge=-1, electrons=0)[0][0]
    print "[M-.]1-", test_two, base, test_two - base
    test_three = mol.isotopic_distribution({"add" : {}, "remove" : {"H" : 2}}, charge=-2, electrons=2)[0][0]
    print "[M-2H]2-", test_three, base, test_three - base
    test_four = mol.isotopic_distribution({"add" : {"K" : 1}, "remove" : {}, "multiply" : 2}, charge=1, electrons=-1)[0][0]
    print "[2M+K]1+", test_four, base, test_four - base
