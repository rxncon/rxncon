from rxncon.core.reaction import *


def test_simple_reaction() -> None:
    rxn = reaction_from_str('A_p+_B_[(r)]')

    assert rxn.terms_rhs == [ReactionTerm([spec_from_str('A')], []),
                             ReactionTerm([spec_from_str('B')], [state_from_str('B_[(r)]-{p}')])]


def test_ppi_reaction() -> None:
    rxn = reaction_from_str('A_[x]_ppi+_B_[y]')

    assert rxn.terms_lhs == [ReactionTerm([spec_from_str('A')], [state_from_str('A_[x]--0')]),
                             ReactionTerm([spec_from_str('B')], [state_from_str('B_[y]--0')])]

    assert rxn.terms_rhs == [ReactionTerm([spec_from_str('A'), spec_from_str('B')], [state_from_str('A_[x]--B_[y]')])]


def test_ipi_reaction() -> None:
    rxn = reaction_from_str('A_[n]_ipi+_A_[m]')

    assert rxn.terms_lhs == [ReactionTerm([spec_from_str('A')], [state_from_str('A_[n]--0'), state_from_str('A_[m]--0')])]
    assert rxn.terms_rhs == [ReactionTerm([spec_from_str('A')], [state_from_str('A_[n]--[m]')])]

