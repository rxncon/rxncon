import pytest

from rxncon.input.quick.quick import Quick
from rxncon.core.effector import StateEffector
from rxncon.core.state import state_from_str
from rxncon.core.reaction import reaction_from_str
from collections import namedtuple
from rxncon.util.utils import elems_eq
from rxncon.core.spec import spec_from_str


ReactionTestCase = namedtuple('ReactionTestCase', ['synthesised_states', 'produced_states', 'consumed_states'])

def test_translation() -> None:
    rxncon_sys = Quick('''A_trsl_BmRNA
                          C_p+_B_[(r1)]
                          D_p+_B_[(r2)] ; ! B@1_[(r1)]-{p}
                          D_[x]_ppi+_B_[y]
                          E_ub+_B_[(r1)]''').rxncon_system

    expected_rxns = {
        'A_trsl_BmRNA': ReactionTestCase(
            [
                state_from_str('B_[(r1)]-{0}'),
                state_from_str('B_[(r2)]-{0}'),
                state_from_str('B_[y]--0')
            ],
            [],
            []
        ),
        'C_p+_B_[(r1)]': ReactionTestCase(
            [],
            [state_from_str('B_[(r1)]-{p}')],
            [state_from_str('B_[(r1)]-{0}')]
        ),
        'D_p+_B_[(r2)]': ReactionTestCase(
            [],
            [state_from_str('B_[(r2)]-{p}')],
            [state_from_str('B_[(r2)]-{0}')]
        ),
        'D_[x]_ppi+_B_[y]': ReactionTestCase(
            [],
            [state_from_str('D_[x]--B_[y]')],
            [state_from_str('D_[x]--0'), state_from_str('B_[y]--0')]
        ),
        'E_ub+_B_[(r1)]': ReactionTestCase(
            [],
            [state_from_str('B_[(r1)]-{ub}')],
            [state_from_str('B_[(r1)]-{0}')]
        )
    }

    for rxn in rxncon_sys.reactions:
        assert elems_eq(rxn.synthesised_states, expected_rxns[str(rxn)].synthesised_states)
        assert elems_eq(rxn.produced_states, expected_rxns[str(rxn)].produced_states)
        assert elems_eq(rxn.consumed_states, expected_rxns[str(rxn)].consumed_states)

    assert rxncon_sys.states_for_component_grouped(spec_from_str('A')) == {}
    assert rxncon_sys.states_for_component_grouped(spec_from_str('C')) == {}
    assert rxncon_sys.states_for_component_grouped(spec_from_str('E')) == {}

    assert elems_eq(list(rxncon_sys.states_for_component_grouped(spec_from_str('B')).values()), [
        [state_from_str('B_[y]--0'), state_from_str('B_[y]--D_[x]')],
        [state_from_str('B_[(r1)]-{0}'), state_from_str('B_[(r1)]-{p}'), state_from_str('B_[(r1)]-{ub}')],
        [state_from_str('B_[(r2)]-{0}'), state_from_str('B_[(r2)]-{p}')]
    ])
    assert elems_eq(list(rxncon_sys.states_for_component_grouped(spec_from_str('D')).values()), [
        [state_from_str('D_[x]--0'), state_from_str('B_[y]--D_[x]')]
    ])


def test_non_elemental_contingency() -> None:
    rxncon_sys = Quick('''A_trsl_BmRNA
                       C_p+_B_[(r1)]
                       D_p+_B_[(r2)]
                       D_[x]_ppi_B_[y] ; ! B-{p}''').rxncon_system

    contingencies = rxncon_sys.contingencies_for_reaction(reaction_from_str('D_[x]_ppi+_B_[y]'))

    assert len(contingencies) == 1
    assert state_from_str('B@1_[(r1)]-{p}') in contingencies[0].effector.states
    assert state_from_str('B@1_[(r2)]-{p}') in contingencies[0].effector.states
    assert contingencies[0].effector.name == 'B-{p}'


def test_non_elemental_contingency_single_state() -> None:
    rxncon_sys = Quick('''A_trsl_BmRNA
                       C_p+_B_[(r1)]
                       D_[x]_ppi_B_[y] ; ! B-{p}''').rxncon_system

    contingencies = rxncon_sys.contingencies_for_reaction(reaction_from_str('D_[x]_ppi+_B_[y]'))

    assert len(contingencies) == 1
    print(contingencies[0].effector)
    assert contingencies[0].effector.name == 'B-{p}'


def test_inconsistent_system() -> None:
    with pytest.raises(AssertionError):
        Quick('''A_ppi_B ; ! <X>
              C_p+_B_[(r1)]
              D_p+_B_[(r2)]
              <X> ; AND B_[(r1)]-{p}
              <X> ; AND B_[(r1)]-{0}''').rxncon_system

    with pytest.raises(AssertionError):
        Quick('''A_ppi_B ; ! B_[(r1)]-{p}
              A_ppi_B ; ! B_[(r1)]-{0}
              C_p+_B_[(r1)]
              D_p+_B_[(r2)]''').rxncon_system


def test_output_reactions() -> None:
    rxncon_sys = Quick("""A_p+_B_[(x)]
                    [output]; ! B_[(x)]-{p}""").rxncon_system

    contingencies = rxncon_sys.contingencies_for_reaction(reaction_from_str('[output]'))

    assert len(contingencies) == 1
    assert isinstance(contingencies[0].effector, StateEffector)
    assert [state_from_str('B@2_[(x)]-{p}')] == contingencies[0].effector.states
