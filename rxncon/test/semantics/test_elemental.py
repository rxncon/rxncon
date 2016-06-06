from rxncon.simulation.rule_based.molecule_from_string import property_ins_from_string, mol_def_from_string,\
    specification_from_string

from rxncon.semantics.elemental import *

def test_empty_domain():
    mol_defs = {
        specification_from_string('A'): mol_def_from_string('A#ass/A_[x]:B_[a]~C_[a]'),
        specification_from_string('B'): mol_def_from_string('B#ass/B_[a]:A_[x]'),
        specification_from_string('C'): mol_def_from_string('C#ass/C_[a]:A_[x]')
    }

    empty_domain = OneParticleElemental(mol_defs,
                                        specification_from_string('A'),
                                        property_ins_from_string(mol_defs[specification_from_string('A')], 'ass/A_[x]:'))

    print(empty_domain.complements())


def test_binding():
    mol_defs = {
        specification_from_string('A'): mol_def_from_string('A#ass/A_[x]:B_[x]~C_[a]'),
        specification_from_string('B'): mol_def_from_string('B#ass/B_[x]:A_[x]~D_[b]'),
        specification_from_string('C'): mol_def_from_string('C#ass/C_[a]:A_[x]'),
        specification_from_string('D'): mol_def_from_string('D#ass/D_[b]:B_[x]')
    }

    binding = TwoParticleElemental(mol_defs, specification_from_string('A'),
                                   property_ins_from_string(mol_defs[specification_from_string('A')], 'ass/A_[x]:B_[x]'),
                                   specification_from_string('B'),
                                   property_ins_from_string(mol_defs[specification_from_string('B')], 'ass/B_[x]:A_[x]'))

    print(binding.complements())
