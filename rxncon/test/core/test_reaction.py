from rxncon.core.reaction import *


def test_simple_reaction():
    rxn = reaction_from_str('A_p+_B_[(r)]')

    assert rxn.terms_rhs == [ReactionTerm(spec_from_str('A'), [], []),
                             ReactionTerm(spec_from_str('B'), [state_from_str('B_[(r)]-{p}')], [])]


def test_ppi_reaction():
    rxn = reaction_from_str('A_[x]_ppi+_B_[y]')

    assert rxn.terms_lhs == [ReactionTerm(spec_from_str('A'), [state_from_str('A_[x]--0')], []),
                             ReactionTerm(spec_from_str('B'), [state_from_str('B_[y]--0')], [])]

    assert rxn.terms_rhs == [ReactionTerm(spec_from_str('A'), [], [spec_from_str('A_[x]~B_[y]')]),
                             ReactionTerm(spec_from_str('B'), [], [spec_from_str('A_[x]~B_[y]')]),
                             BondReactionTerm(spec_from_str('A_[x]~B_[y]'), [state_from_str('A_[x]--B_[y]')])]


def test_ipi_reaction():
    rxn = reaction_from_str('A_[n]_ipi_A_[m]')

    assert rxn.terms_lhs == [ReactionTerm(spec_from_str('A'), [state_from_str('A_[n]--0'),
                                                               state_from_str('A_[m]--0')], [])]

    assert rxn.terms_rhs == [ReactionTerm(spec_from_str('A'), [], [spec_from_str('A_[n]~A_[m]')]),
                             BondReactionTerm(spec_from_str('A_[n]~A_[m]'), [state_from_str('A_[n]--[m]')])]


def test_wrong():
    rxn = reaction_from_str('A_ppi+_B')
    assert rxn.produced_states == [state_from_str('A_[B]--B_[A]')]

    rxn = reaction_from_str('A_trsl_B')
    assert ReactionTerm(spec_from_str('BmRNA'), [], []) in rxn.terms_rhs
