from pyisopach import Molecule

mol = Molecule("C12H22O11")


print(mol.molecular_weight)

print(mol.num_atoms)

rule_dict = {
    "add" : {
        "H" : 2,
        "Cl" : 1
    },
    "remove" : {
        "O" : 10
    },
    "divide" :{
        "H": 2
    },
    "multiply" : {
        "O": 10
    }
}

print(mol._apply_rule_dict(rule_dict))
