from rxncon.core.spec import spec_from_str
from rxncon.core.state import state_from_str
from rxncon.core.reaction import reaction_from_str
from rxncon.input.quick.quick import Quick
from rxncon.core.contingency import ContingencyType


def test_single_reaction() -> None:
    rxncon_sys = Quick('A_p+_B_[(r)]').rxncon_system
    assert state_from_str('B_[(r)]-{0}') in rxncon_sys.consumed_states
    assert state_from_str('B_[(r)]-{p}') in rxncon_sys.produced_states

    assert state_from_str('B_[(r)]-{0}') in rxncon_sys.states_for_component(spec_from_str('B'))
    assert state_from_str('B_[(r)]-{p}') in rxncon_sys.states_for_component(spec_from_str('B'))


def test_single_contingency() -> None:
    rxncon_sys = Quick('''A_p+_B_[(r)]; ! A_[(x)]-{p}
                       C_p+_A_[(x)]''').rxncon_system
    assert state_from_str('A_[(x)]-{0}') in rxncon_sys.consumed_states
    assert state_from_str('B_[(r)]-{0}') in rxncon_sys.consumed_states
    assert state_from_str('A_[(x)]-{p}') in rxncon_sys.produced_states
    assert state_from_str('B_[(r)]-{p}') in rxncon_sys.produced_states

    assert state_from_str('A_[(x)]-{0}') in rxncon_sys.states_for_component(spec_from_str('A'))
    assert state_from_str('B_[(r)]-{0}') in rxncon_sys.states_for_component(spec_from_str('B'))
    assert state_from_str('A_[(x)]-{p}') in rxncon_sys.states_for_component(spec_from_str('A'))
    assert state_from_str('B_[(r)]-{p}') in rxncon_sys.states_for_component(spec_from_str('B'))

    contingencies = rxncon_sys.contingencies_for_reaction(reaction_from_str('A_p+_B_[(r)]'))
    assert len(contingencies) == 1
    assert contingencies[0].effector.states == [state_from_str('A@0_[(x)]-{p}')]
    assert contingencies[0].contingency_type == ContingencyType.requirement


def test_bidirectional_reactions() -> None:
    rxncon_sys = Quick('''A_[x]_ppi_B_[y]; ! A_[(x)]-{p}
                       C_p+_A_[(x)]''').rxncon_system

    #  Forward direction
    contingencies = rxncon_sys.contingencies_for_reaction(reaction_from_str('A_[x]_ppi+_B_[y]'))
    assert len(contingencies) == 1
    assert contingencies[0].effector.states == [state_from_str('A@0_[(x)]-{p}')]
    assert contingencies[0].contingency_type == ContingencyType.requirement

    #  Reverse direction
    contingencies = rxncon_sys.contingencies_for_reaction(reaction_from_str('A_[x]_ppi-_B_[y]'))
    assert len(contingencies) == 0
