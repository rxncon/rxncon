from rxncon.core.spec import spec_from_string
from rxncon.core.state import state_from_string
from rxncon.input.quick.quick import *
from rxncon.util.utils import elems_eq


def test_simple_quick():
    rxncon_sys = Quick('A_p+_B_[(r)]').rxncon_system
    assert state_from_string('B_[(r)]-{0}') in rxncon_sys.consumed_states
    assert state_from_string('B_[(r)]-{p}') in rxncon_sys.produced_states

    assert state_from_string('B_[(r)]-{0}') in rxncon_sys.states_for_component(spec_from_string('B'))
    assert state_from_string('B_[(r)]-{p}') in rxncon_sys.states_for_component(spec_from_string('B'))

    assert len(rxncon_sys.states_for_component_grouped(spec_from_string('B'))) == 1
    assert len(rxncon_sys.states_for_component_grouped(spec_from_string('B'))[0]) == 2
    assert state_from_string('B_[(r)]-{0}') in rxncon_sys.states_for_component_grouped(spec_from_string('B'))[0]
    assert state_from_string('B_[(r)]-{p}') in rxncon_sys.states_for_component_grouped(spec_from_string('B'))[0]

    return elems_eq(rxncon_sys.states_for_component_grouped(spec_from_string('B')),
                    [[state_from_string('B_[(r)]-{0}'), state_from_string('B_[(r)]-{p}')]])


def test_with_cont():
    rxncon_sys = Quick('''A_p+_B_[(r)]; ! A_[(x)]-{p}
                       C_p+_A_[(x)]''').rxncon_system
    assert state_from_string('A_[(x)]-{0}') in rxncon_sys.consumed_states
    assert state_from_string('B_[(r)]-{0}') in rxncon_sys.consumed_states
    assert state_from_string('A_[(x)]-{p}') in rxncon_sys.produced_states
    assert state_from_string('B_[(r)]-{p}') in rxncon_sys.produced_states

    assert state_from_string('A_[(x)]-{0}') in rxncon_sys.states_for_component(spec_from_string('A'))
    assert state_from_string('B_[(r)]-{0}') in rxncon_sys.states_for_component(spec_from_string('B'))
    assert state_from_string('A_[(x)]-{p}') in rxncon_sys.states_for_component(spec_from_string('A'))
    assert state_from_string('B_[(r)]-{p}') in rxncon_sys.states_for_component(spec_from_string('B'))

    assert rxncon_sys.contingencies_for_reaction(reaction_from_string('A_p+_B_[(r)]'))


def test_ppi():
    rxncon_sys = Quick('''A_[x]_ppi_B_[y]''').rxncon_system
    assert state_from_string('A_[x]--0') in rxncon_sys.consumed_states
    assert state_from_string('B_[y]--0') in rxncon_sys.consumed_states
    assert state_from_string('A_[x]--B_[y]') in rxncon_sys.produced_states

    assert state_from_string('A_[x]--0') in rxncon_sys.states_for_component(spec_from_string('A'))
    assert state_from_string('B_[y]--0') in rxncon_sys.states_for_component(spec_from_string('B'))
    assert state_from_string('A_[x]--B_[y]') in rxncon_sys.states_for_component(spec_from_string('A'))
    assert state_from_string('A_[x]--B_[y]') in rxncon_sys.states_for_component(spec_from_string('B'))


def test_ppi_with_cont():
    rxncon_sys = Quick('''A_[n]_ppi_B_[m]; ! A_[(x)]-{p}
                       C_p+_A_[(x)]''').rxncon_system

    assert elems_eq(rxncon_sys.states_for_component_grouped(spec_from_string('A')),
                    [[state_from_string('A_[n]--0'), state_from_string('A_[n]--B_[m]')],
                     [state_from_string('A_[(x)]-{0}'), state_from_string('A_[(x)]-{p}')]])

