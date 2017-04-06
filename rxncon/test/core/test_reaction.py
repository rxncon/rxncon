import pytest
from rxncon.core.reaction import reaction_from_str, ReactionTerm, initialize_reaction_defs
from rxncon.core.state import state_from_str
from rxncon.core.spec import spec_from_str


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


def test_output_reaction() -> None:
    rxn = reaction_from_str('[Output]')

    assert rxn.components_rhs == []
    assert rxn.components_lhs == []
    assert rxn.degraded_components == []
    assert rxn.synthesised_components == []

    assert rxn.consumed_states == []
    assert rxn.produced_states == []
    assert rxn.synthesised_states == []
    assert rxn.degraded_states == []


def test_equality_output_reaction() -> None:
    assert not reaction_from_str('A_p+_B_[(x)]') == reaction_from_str('[output]')
    assert not reaction_from_str('[output]') == reaction_from_str('A_p+_B_[(x)]')
    assert reaction_from_str('[output]') == reaction_from_str('[output]')
    assert not reaction_from_str('[output]') == reaction_from_str('[output2]')


def test_modifier() -> None:
    rxn = reaction_from_str('A_trsl_BmRNA')
    assert rxn.modifier_components == [spec_from_str('A'), spec_from_str('BmRNA')]
    assert rxn.modifier_states == []


def test_dynamical_reactions() -> None:
    with pytest.raises(SyntaxError):
        rxn = reaction_from_str('A_agex_A_[(r)]')

    initialize_reaction_defs([
        {
            '!UID:Reaction': 'auto-GuanineNucleotideExchange',
            '!UID:ReactionKey': 'agex',
            '!BidirectionalVerb': 'no',
            '!MolTypeX': 'Protein',
            '!ResolutionX': 'component',
            '!MolTypeY': 'Protein',
            '!ResolutionY': 'residue',
            '!SkeletonRule': '$y%#$y%-{0} -> $y%#$y%-{GTP}'
        }
    ])

    rxn = reaction_from_str('A_agex_A_[(r)]')
    assert rxn.produced_states == [state_from_str('A_[(r)]-{gtp}')]
    assert rxn.consumed_states == [state_from_str('A_[(r)]-{0}')]

    initialize_reaction_defs()

    with pytest.raises(SyntaxError):
        rxn = reaction_from_str('A_agex_A_[(r)]')
