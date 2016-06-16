from rxncon.simulation.rule_based.rbm_from_rxncon import *
from rxncon.semantics.elemental import *
from rxncon.semantics.molecule import *
from rxncon.simulation.rule_based.molecule_from_string import *
from rxncon.syntax.rxncon_from_string import *
from rxncon.venntastic.sets import *
from rxncon.core.reaction import *
from rxncon.input.quick.quick import *


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

    contingency_soln = PropertySet(elemental_from_state(mol_defs, state_from_string('B_[(r)]-{p}')))
    reaction = reaction_from_string(REACTION_DEFINITIONS, 'A_[x]_ppi_B_[y]')

    rule = rule_from_reaction_and_contingency_soln(mol_defs, reaction, contingency_soln)

    print(rule)


def test_rule_from_rxn_and_contingency_soln_or_contingency():
    mol_defs = {
        specification_from_string('A'): mol_def_from_string('A#ass/A_[x]:B_[y],ass/A_[c]:C_[a]'),
        specification_from_string('B'): mol_def_from_string('B#mod/B_[(r)]:0~p,ass/B_[y]:A_[x]'),
        specification_from_string('C'): mol_def_from_string('C#ass/C_[a]:A_[c]')
    }

    contingency_soln = Union(
        PropertySet(state_from_string('A_[c]--C_[a]')),
        PropertySet(state_from_string('B_[(r)]-{p}'))
    )
    reaction = reaction_from_string(REACTION_DEFINITIONS, 'A_[x]_ppi_B_[y]')

    rules = rules_from_mol_defs_and_reaction_and_bg_state_set(mol_defs, reaction, contingency_soln)

    for rule in rules:
        print()
        print(rule)


def test_complicated_case():
    mol_defs = {
        specification_from_string('A'): mol_def_from_string('A#ass/A_[b]:B_[a],ass/A_[c]:C_[a]'),
        specification_from_string('B'): mol_def_from_string('B#ass/B_[a]:A_[b],ass/B_[e]:E_[b]'),
        specification_from_string('C'): mol_def_from_string('C#ass/C_[a]:A_[c],ass/C_[d]:D_[c]'),
        specification_from_string('D'): mol_def_from_string('D#ass/D_[c]:C_[d]'),
        specification_from_string('E'): mol_def_from_string('E#ass/E_[b]:B_[e]')
    }


    # cont_soln = Complement(Union(
    #     Intersection(PropertySet(state_from_string('A_[c]--C_[a]')),
    #                  PropertySet(state_from_string('C_[d]--D_[c]'))),
    #     Intersection(PropertySet(state_from_string('A_[c]--C_[a]')),
    #                  PropertySet(state_from_string('B_[e]--E_[b]')))))

    cont_soln = Complement(Union(
        Intersection(PropertySet(state_from_string('A_[c]--C_[a]')),
                     PropertySet(state_from_string('C_[d]--D_[c]'))),
        Intersection(PropertySet(state_from_string('A_[c]--C_[a]')),
                     PropertySet(state_from_string('B_[e]--E_[b]')))))

    reaction = reaction_from_string(REACTION_DEFINITIONS, 'A_[b]_ppi_B_[a]')

    rules = rules_from_mol_defs_and_reaction_and_bg_state_set(mol_defs, reaction, cont_soln)

    for i, rule in enumerate(rules):
        print()
        print(i)
        print()
        print(rule)


def test_rule_from_rxncon_system():
    rxncon = Quick('C_p+_A_[(r)]\nA_[x]_ppi_B_[y] ; ! A_[(r)]-{p}').rxncon_system

    rules = []

    for reaction in rxncon.reactions:
        rules += rules_from_reaction(rxncon, reaction)

    print(rules)


def test_rules_from_everything():
    pass

