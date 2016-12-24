from rxncon.core.spec import spec_from_str
from rxncon.core.state import state_from_str
from rxncon.input.quick.quick import *


def test_single_reaction():
    rxncon_sys = Quick('A_p+_B_[(r)]').rxncon_system
    assert state_from_str('B_[(r)]-{0}') in rxncon_sys.consumed_states
    assert state_from_str('B_[(r)]-{p}') in rxncon_sys.produced_states

    assert state_from_str('B_[(r)]-{0}') in rxncon_sys.states_for_component(spec_from_str('B'))
    assert state_from_str('B_[(r)]-{p}') in rxncon_sys.states_for_component(spec_from_str('B'))


def test_single_contingency():
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

    assert len(rxncon_sys.contingencies_for_reaction(reaction_from_str('A_p+_B_[(r)]'))) == 1

