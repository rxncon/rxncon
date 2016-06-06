from rxncon.simulation.rule_based.molecule_from_string import property_ins_from_string, mol_def_from_string,\
    specification_from_string

from rxncon.semantics.fact import *

def test_one():
    mol_defs = {
        specification_from_string('A'): mol_def_from_string('A#ass/A_[x]:B_[a]~C_[a]'),
        specification_from_string('B'): mol_def_from_string('B#ass/B_[a]:A_[x]'),
        specification_from_string('C'): mol_def_from_string('C#ass/C_[a]:A_[x]')
    }

    empty_domain = OneParticleFact(mol_defs,
                                   specification_from_string('A'),
                                   property_ins_from_string(mol_defs[specification_from_string('A')],
                                                            'ass/A_[x]:'))

    print(empty_domain.complements())
