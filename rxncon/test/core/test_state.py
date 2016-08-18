from rxncon.core.state import state_from_string
from rxncon.core.spec import bond_spec_from_string, mol_spec_from_string

def test_ppi_states():
    state = state_from_string('A--B')
    assert state.target == bond_spec_from_string('A~B')
    assert not state.is_elemental

    elem_state = state_from_string('A_[d1]--B_[d2]')
    assert elem_state.target == bond_spec_from_string('A_[d1]~B_[d2]')
    assert elem_state.is_elemental

    assert state.is_superset_of(elem_state)
    assert elem_state.is_subset_of(state)


def test_ppi_equality():
    print(state_from_string('A--B'))
    print(state_from_string('B--A'))



def test_empty_ppi_states():
    state = state_from_string('A_[x]--0')
    assert state.target == mol_spec_from_string('A_[x]')
    assert state.is_elemental

    state = state_from_string('0--B_[x]')
    assert state.target == mol_spec_from_string('B_[x]')
    assert state.is_elemental


def test_ipi_states():
    state = state_from_string('A_[n]--[m]')
    assert state.target == bond_spec_from_string('A_[m]~A_[n]')
    assert state.is_elemental


def test_mod_states():
    state = state_from_string('A_[(r)]-{0}')
    assert state.target == mol_spec_from_string('A_[(r)]')
    assert state.is_elemental
    assert state.is_subset_of(state_from_string('A-{0}'))


def test_neutral_states():
    state = state_from_string('A_[(r)]-{p}')
    assert not state.is_neutral

    neutral = state_from_string('A_[(r)]-{0}')
    assert neutral.is_neutral
    assert state.neutral_states == [neutral]

    state = state_from_string('A_[x]--B_[y]')
    assert not state.is_neutral

    first_neutral = state_from_string('A_[x]--0')
    second_neutral = state_from_string('B_[y]--0')
    assert first_neutral.is_neutral
    assert second_neutral.is_neutral

    assert first_neutral in state.neutral_states
    assert second_neutral in state.neutral_states

    assert first_neutral.neutral_states == [first_neutral]
    assert second_neutral.neutral_states == [second_neutral]
