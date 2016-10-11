from rxncon.core.state import state_from_str, FullyNeutralState
from rxncon.core.spec import bond_spec_from_str, spec_from_str
from rxncon.util.utils import elems_eq
import pytest

###                      ###
# Test modification states #
###                      ###

def test_modification_props():
    # Elemental state, neutral.
    state = state_from_str('A_[(res)]-{0}')
    assert state.is_elemental
    assert elems_eq(state.components, [spec_from_str('A')])
    assert state.target == spec_from_str('A_[(res)]')
    assert state.is_neutral
    assert elems_eq(state.neutral_states, [state])

    # Elemental state, non-neutral.
    state = state_from_str('A_[(res)]-{p}')
    assert state.is_elemental
    assert elems_eq(state.components, [spec_from_str('A')])
    assert state.target == spec_from_str('A_[(res)]')
    assert not state.is_neutral
    assert elems_eq(state.neutral_states, [state_from_str('A_[(res)]-{0}')])

    # Non-elemental state, neutral.
    state = state_from_str('A_[dom]-{0}')
    assert not state.is_elemental
    assert elems_eq(state.components, [spec_from_str('A')])
    assert state.target == spec_from_str('A_[dom]')
    assert state.is_neutral
    assert elems_eq(state.neutral_states, [state])

    # Non-elemental state, non-neutral.
    state = state_from_str('A_[dom]-{p}')
    assert not state.is_elemental
    assert elems_eq(state.components, [spec_from_str('A')])
    assert state.target == spec_from_str('A_[dom]')
    assert not state.is_neutral
    assert elems_eq(state.neutral_states, [state_from_str('A_[dom]-{0}')])


def test_modification_parsing():
    # Upper/lower case should not matter for StateModifier.
    assert state_from_str('A-{P}') == state_from_str('A-{p}')
    assert state_from_str('A_[(res)]-{UB}') == state_from_str('A_[(res)]-{Ub}')

    # Unknown StateModifier should raise.
    with pytest.raises(ValueError):
        state_from_str('A_[(res)]-{kryptonite}')


def test_modification_superset_subset():
    # Happy path, same modification.
    assert state_from_str('A-{p}').is_superset_of(state_from_str('A_[dom]-{p}'))
    assert state_from_str('A_[dom]-{p}').is_superset_of(state_from_str('A_[dom(res)]-{p}'))
    assert state_from_str('A_[dom(res)]-{p}').is_subset_of(state_from_str('A_[dom]-{p}'))
    assert state_from_str('A_[dom]-{p}').is_subset_of(state_from_str('A-{p}'))

    # Sad path, same modification.
    assert not state_from_str('A-{p}').is_subset_of(state_from_str('A_[dom]-{p}'))
    assert not state_from_str('A_[dom]-{p}').is_subset_of(state_from_str('A_[dom(res)]-{p}'))
    assert not state_from_str('A_[dom(res)]-{p}').is_superset_of(state_from_str('A_[dom]-{p}'))
    assert not state_from_str('A_[dom]-{p}').is_superset_of(state_from_str('A-{p}'))

    # Different modifications.
    assert not state_from_str('A-{p}').is_superset_of(state_from_str('A_[dom]-{0}'))
    assert not state_from_str('A_[dom(res)]-{ub}').is_subset_of(state_from_str('A_[dom]-{p}'))


def test_modification_structured_index():
    state = state_from_str('A@3-{p}')
    assert state.components[0].struct_index == 3

###                                     ###
# Test protein-protein-interaction states #
###                                     ###

def test_ppi_props():
    # Elemental state, free binding domain.
    state = state_from_str('A_[m]--0')
    assert state.is_elemental
    assert elems_eq(state.components, [spec_from_str('A')])
    assert state.target == spec_from_str('A_[m]')
    assert state.is_neutral
    assert elems_eq(state.neutral_states, [state])

    # Elemental state, bond.
    state = state_from_str('A_[m]--B_[n]')
    assert state.is_elemental
    assert elems_eq(state.components, [spec_from_str('A'), spec_from_str('B')])
    assert state.target == bond_spec_from_str('A_[m]~B_[n]')
    assert not state.is_neutral
    assert elems_eq(state.neutral_states, [state_from_str('A_[m]--0'), state_from_str('B_[n]--0')])

    # Non-elemental state, free binding domain.
    state = state_from_str('A--0')
    assert not state.is_elemental
    assert elems_eq(state.components, [spec_from_str('A')])
    assert state.target == spec_from_str('A')
    assert state.is_neutral
    assert elems_eq(state.neutral_states, [state])

    # Non-elemental state, bond.
    state = state_from_str('A--B_[n]')
    assert not state.is_elemental
    assert elems_eq(state.components, [spec_from_str('A'), spec_from_str('B')])
    assert state.target == bond_spec_from_str('A~B_[n]')
    assert not state.is_neutral
    # @todo JCR 20160920 This is not correct, the state A--0 should be replaced by the list of domains of A that
    # @todo              can bind to B. Perhaps we should invent a new notation for this, since this property cannot
    # @todo              be known by just looking at the state itself, but requires knowledge of the rest of the system.
    assert elems_eq(state.neutral_states, [state_from_str('A--0'), state_from_str('B_[n]--0')])


def test_ppi_parsing():
    # @todo JCR 20160920 Bonds are symmetric. We need extra structure to allow this is the StateDef.
    # assert state_from_string('A_[x]--B_[y]') == state_from_string('B_[y]--A_[x]')

    # Too fine resolution (higher than elemental) raises.
    with pytest.raises(SyntaxError):
        state_from_str('A_[(x)]--B_[(y)]')


def test_ppi_superset_subset():
    # Happy path, free binding domain.
    assert state_from_str('A--0').is_superset_of(state_from_str('A_[m]--0'))
    assert state_from_str('A_[m]--0').is_subset_of(state_from_str('A--0'))

    # Sad path, free binding domain.
    assert not state_from_str('A--0').is_subset_of(state_from_str('A_[m]--0'))
    assert not state_from_str('A_[m]--0').is_superset_of(state_from_str('A--0'))

    # Happy path, superset, bond.
    assert state_from_str('A--B').is_superset_of(state_from_str('A_[m]--B_[n]'))
    assert state_from_str('A--B').is_superset_of(state_from_str('A_[m]--B'))
    assert state_from_str('A--B').is_superset_of(state_from_str('A--B_[n]'))
    assert state_from_str('A_[m]--B').is_superset_of(state_from_str('A_[m]--B_[n]'))
    assert state_from_str('A--B_[n]').is_superset_of(state_from_str('A_[m]--B_[n]'))

    # Happy path, subset, bond.
    assert state_from_str('A_[m]--B_[n]').is_subset_of(state_from_str('A--B'))
    assert state_from_str('A_[m]--B').is_subset_of(state_from_str('A--B'))
    assert state_from_str('A--B_[n]').is_subset_of(state_from_str('A--B'))
    assert state_from_str('A_[m]--B_[n]').is_subset_of(state_from_str('A_[m]--B'))
    assert state_from_str('A_[m]--B_[n]').is_subset_of(state_from_str('A--B_[n]'))

    # Sad path, superset, bond.
    assert not state_from_str('A--B').is_subset_of(state_from_str('A_[m]--B_[n]'))
    assert not state_from_str('A--B').is_subset_of(state_from_str('A_[m]--B'))
    assert not state_from_str('A--B').is_subset_of(state_from_str('A--B_[n]'))
    assert not state_from_str('A_[m]--B').is_subset_of(state_from_str('A_[m]--B_[n]'))
    assert not state_from_str('A--B_[n]').is_subset_of(state_from_str('A_[m]--B_[n]'))

    # Sad path, subset, bond.
    assert not state_from_str('A_[m]--B_[n]').is_superset_of(state_from_str('A--B'))
    assert not state_from_str('A_[m]--B').is_superset_of(state_from_str('A--B'))
    assert not state_from_str('A--B_[n]').is_superset_of(state_from_str('A--B'))
    assert not state_from_str('A_[m]--B_[n]').is_superset_of(state_from_str('A_[m]--B'))
    assert not state_from_str('A_[m]--B_[n]').is_superset_of(state_from_str('A--B_[n]'))


def test_ppi_structured_index():
    state = state_from_str('X@4--Z@3')
    assert elems_eq([(spec.component_name, spec.struct_index) for spec in state.mol_specs], [('X', 4), ('Z', 3)])

###                          ###
# Test self-interaction states #
###                          ###

def test_ipi_props():
    # NOTE : States with free binding domain are the same states as tested in 'test_ppi_props', tested there.

    # Elemental state, bond.
    state = state_from_str('A_[m]--[n]')
    assert state.is_elemental
    assert elems_eq(state.components, [spec_from_str('A')])
    assert state.target == bond_spec_from_str('A_[m]~A_[n]')
    assert not state.is_neutral
    assert elems_eq(state.neutral_states, [state_from_str('A_[m]--0'), state_from_str('A_[n]--0')])

    # Non-elemental state, bond.
    state = state_from_str('A--A_[n]')
    assert not state.is_elemental
    assert elems_eq(state.components, [spec_from_str('A')])
    assert state.target == bond_spec_from_str('A~A_[n]')
    assert not state.is_neutral
    # @todo JCR 20160920 This is not correct, the state A--0 should be replaced by the list of domains of A that
    # @todo              can bind to B. Perhaps we should invent a new notation for this, since this property cannot
    # @todo              be known by just looking at the state itself, but requires knowledge of the rest of the system.
    # @todo              Furthermore, the state A--0 is a superset of the other one, so it should not appear independently.
    assert elems_eq(state.neutral_states, [state_from_str('A--0'), state_from_str('A_[n]--0')])


def test_ipi_parsing():
    # @todo JCR 20160920 Bonds are symmetric. We need extra structure to capture this is the StateDef.
    # assert state_from_string('A_[x]--B_[y]') == state_from_string('B_[y]--A_[x]')

    # Too fine resolution (higher than elemental) raises.
    with pytest.raises(SyntaxError):
        state_from_str('A_[(x)]--[(y)]')


def test_ipi_superset_subset():
    # NOTE : States with free binding domain are the same states as tested in 'test_ppi_superset_subset', tested there.

    # @todo JCR 20160920 There is a problem here: the StateDef requires the second variable in the State string to be
    # @todo              a Locus (i.e. '[n]') and not a full Spec (i.e. 'A_[n]').
    pass


###                      ###
# Test Fully Neutral state #
###                      ###

def test_fully_neutral():
    assert state_from_str('0') == FullyNeutralState()


def test_global():
    x = state_from_str('[IN]')
    print(x)
