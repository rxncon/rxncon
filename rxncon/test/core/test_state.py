import pytest

from rxncon.core.state import State, StateDef, StateModifier, STATE_DEFS, state_from_string, mol_spec_from_string


def test_ppi_states():
    state = state_from_string(STATE_DEFS, 'A--B')
    assert state.components == [mol_spec_from_string('A'), mol_spec_from_string('B')]
    assert not state.is_elemental

    elem_state = state_from_string(STATE_DEFS, 'A_[d1]--B_[d2]')
    assert elem_state.components == [mol_spec_from_string('A_[d1]'), mol_spec_from_string('B_[d2]')]
    assert elem_state.is_elemental

    assert state.is_superset_of(elem_state)
    assert elem_state.is_subset_of(state)


def test_ipi_states():
    state = state_from_string(STATE_DEFS, 'A_[n]--[m]')
    print(state.components)
