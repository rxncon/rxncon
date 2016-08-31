from rxncon.core.reaction import *

def test_simple_reaction():
    rxn = reaction_from_string('A_p+_B_[(r)]')

    assert rxn.reactants_rhs == [Reactant(spec_from_string('A'), None),
                                 Reactant(spec_from_string('B'), [state_from_string('B_[(r)]-{p}')])]

def test_ppi_reaction():
    rxn = reaction_from_string('A_[x]_ppi_B_[y]')

    assert rxn.reactants_lhs == [Reactant(spec_from_string('A'), [state_from_string('A_[x]--0')]),
                                 Reactant(spec_from_string('B'), [state_from_string('B_[y]--0')])]

    assert rxn.reactants_rhs == [Reactant(spec_from_string('A'), spec_from_string('A_[x]~B_[y]')),
                                 Reactant(spec_from_string('B'), spec_from_string('A_[x]~B_[y]')),
                                 Reactant(spec_from_string('A_[x]~B_[y]'), [state_from_string('A_[x]--B_[y]')])]


def test_ipi_reaction():
    rxn = reaction_from_string('A_[n]_ipi_A_[m]')

    assert rxn.reactants_lhs == [Reactant(spec_from_string('A'), [state_from_string('A_[n]--0'),
                                                                  state_from_string('A_[m]--0')])]

    assert rxn.reactants_rhs == [Reactant(spec_from_string('A'), spec_from_string('A_[n]~A_[m]')),
                                 Reactant(spec_from_string('A_[n]~A_[m]'), [state_from_string('A_[n]--[m]')])]

def test_wrong():
    rxn = reaction_from_string('A_ppi_B')
    assert rxn.produced_states == [state_from_string('A_[B]--B_[A]')]

    rxn = reaction_from_string('A_trsl_B')
    assert Reactant(spec_from_string('BmRNA'), None) in rxn.reactants_rhs
