from rxncon.input.quick.quick import *
from rxncon.core.state import state_from_string

def test_simple_quick():
    rxncon_sys = Quick('A_p+_B_[(r)]').rxncon_system
    assert state_from_string('B_[(r)]-{0}') in rxncon_sys.consumed_states
    assert state_from_string('B_[(r)]-{p}') in rxncon_sys.produced_states


def test_with_cont():
    rxncon_sys = Quick('''A_p+_B_[(r)]; ! A_[(x)]-{p}
                       C_p+_A_[(x)]''').rxncon_system
    assert state_from_string('A_[(x)]-{0}') in rxncon_sys.consumed_states
    assert state_from_string('B_[(r)]-{0}') in rxncon_sys.consumed_states
    assert state_from_string('A_[(x)]-{p}') in rxncon_sys.produced_states
    assert state_from_string('B_[(r)]-{p}') in rxncon_sys.produced_states
    assert rxncon_sys.contingencies_for_reaction(reaction_from_string('A_p+_B_[(r)]'))


def test_ppi():
    rxncon_sys = Quick('''A_[x]_ppi_B_[y]''').rxncon_system
    assert state_from_string('A_[x]--0') in rxncon_sys.consumed_states
    assert state_from_string('B_[y]--0') in rxncon_sys.consumed_states
    assert state_from_string('A_[x]--B_[y]') in rxncon_sys.produced_states
