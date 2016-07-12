import pytest

from rxncon.core.state import State, StateDef, StateModifier, STATE_DEFS, state_from_string, spec_from_string


def test_states():
    state = state_from_string(STATE_DEFS, 'A--B')
    assert state.components == [spec_from_string('A'), spec_from_string('B')]
    assert not state.is_elemental

    elem_state = state_from_string(STATE_DEFS, 'A_[d1]--B_[d2]')
    assert elem_state.components == [spec_from_string('A_[d1]'), spec_from_string('B_[d2]')]
    assert elem_state.is_elemental

    assert state.is_superset_of(elem_state)
    assert elem_state.is_subset_of(state)
