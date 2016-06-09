from rxncon.simulation.rule_based.rbm_from_rxncon import *
from rxncon.semantics.elemental import *
from rxncon.semantics.molecule import *
from rxncon.simulation.rule_based.molecule_from_string import *
from rxncon.syntax.rxncon_from_string import *
from rxncon.venntastic.sets import *
from rxncon.core.reaction import *

def test_complexes_from_elementals_no_bg():
    mol_defs = {
        specification_from_string('A'): mol_def_from_string('A#'),
        specification_from_string('B'): mol_def_from_string('B#mod/B_[(r)]:0~p')
    }

    elem_a = OneParticleElemental(mol_defs, specification_from_string('A'), None)
    elem_b = elemental_from_state(mol_defs, state_from_string('B_[(r)]-{0}'))

    complexes = complexes_from_elementals(mol_defs, [elem_a, elem_b], [])

    print(complexes)

def test_complexes_from_elementals_bg():
    mol_defs = {
        specification_from_string('A'): mol_def_from_string('A#mod/A_[(x)]:0~p'),
        specification_from_string('B'): mol_def_from_string('B#mod/B_[(r)]:0~p')
    }

    elem_a = OneParticleElemental(mol_defs, specification_from_string('A'), None)
    elem_b = elemental_from_state(mol_defs, state_from_string('B_[(r)]-{0}'))

    elem_a_phos = elemental_from_state(mol_defs, state_from_string('A_[(x)]-{p}'))

    complexes = complexes_from_elementals(mol_defs, [elem_a, elem_b], [elem_a_phos])

    print(complexes)


def test_complexes_from_elementals_with_binding():
    mol_defs = {
        specification_from_string('A'): mol_def_from_string('A#ass/A_[x]:B_[y]'),
        specification_from_string('B'): mol_def_from_string('B#ass/B_[y]:A_[x]')
    }

    elem_ab = elemental_from_state(mol_defs, state_from_string('A_[x]--B_[y]'))

    complexes = complexes_from_elementals(mol_defs, [elem_ab], [])

    print(complexes)


def test_complexes_binding_mutual_exclusivity_returns_none():
    mol_defs = {
        specification_from_string('A'): mol_def_from_string('A#ass/A_[x]:B_[y]~C_[z]'),
        specification_from_string('B'): mol_def_from_string('B#ass/B_[y]:A_[x]'),
        specification_from_string('C'): mol_def_from_string('C#ass/C_[z]:A_[x]')
    }

    elem_ab = elemental_from_state(mol_defs, state_from_string('A_[x]--B_[y]'))
    elem_ac = elemental_from_state(mol_defs, state_from_string('A_[x]--C_[z]'))

    complexes = complexes_from_elementals(mol_defs, [elem_ab], [elem_ac])

    print(complexes)


def test_rule_from_rxn_and_contingency_soln_single_contingency():
    mol_defs = {
        specification_from_string('A'): mol_def_from_string('A#ass/A_[x]:B_[y]'),
        specification_from_string('B'): mol_def_from_string('B#mod/B_[(r)]:0~p,ass/B_[y]:A_[x]')
    }

    contingency_soln = PropertySet(elemental_from_state(mol_defs, state_from_string('B_[(r)]-{0}')))
    reaction = reaction_from_string(REACTION_DEFINITIONS, 'A_[x]_ppi_B_[y]')

    rule = rule_from_reaction_and_contingency_soln(mol_defs, reaction, contingency_soln)

    print(rule)





