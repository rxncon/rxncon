from rxncon.core.state import state_from_str, FullyNeutralState, GlobalState, EmptyBindingState
from rxncon.core.spec import spec_from_str
from rxncon.util.utils import elems_eq
import pytest


###                      ###
# Test modification states #
###                      ###

def test_modification_props() -> None:
    # Elemental state, neutral.
    state = state_from_str('A_[(res)]-{0}')
    assert state.is_elemental
    assert elems_eq(state.components, [spec_from_str('A')])
    assert state.is_neutral
    assert elems_eq(state.neutral_states, [state])

    # Elemental state, non-neutral.
    state = state_from_str('A_[(res)]-{p}')
    assert state.is_elemental
    assert elems_eq(state.components, [spec_from_str('A')])
    assert not state.is_neutral
    assert elems_eq(state.neutral_states, [state_from_str('A_[(res)]-{0}')])

    # Non-elemental state, neutral.
    state = state_from_str('A_[dom]-{0}')
    assert not state.is_elemental
    assert elems_eq(state.components, [spec_from_str('A')])
    assert state.is_neutral
    assert elems_eq(state.neutral_states, [state])

    # Non-elemental state, non-neutral.
    state = state_from_str('A_[dom]-{p}')
    assert not state.is_elemental
    assert elems_eq(state.components, [spec_from_str('A')])
    assert not state.is_neutral
    assert elems_eq(state.neutral_states, [state_from_str('A_[dom]-{0}')])


def test_modification_parsing() -> None:
    # Upper/lower case should not matter for StateModifier.
    assert state_from_str('A-{P}') == state_from_str('A-{p}')
    assert state_from_str('A_[(res)]-{UB}') == state_from_str('A_[(res)]-{Ub}')

    # Unknown StateModifier should raise.
    with pytest.raises(ValueError):
        state_from_str('A_[(res)]-{kryptonite}')


def test_modification_superset_subset() -> None:
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


def test_modification_structured_index() -> None:
    state = state_from_str('A@3-{p}')
    assert state.components[0].struct_index == 3

###                                     ###
# Test protein-protein-interaction states #
###                                     ###

def test_ppi_props() -> None:
    # Elemental state, free binding domain.
    state = state_from_str('A_[m]--0')
    assert state.is_elemental
    assert elems_eq(state.components, [spec_from_str('A')])
    assert state.is_neutral
    assert elems_eq(state.neutral_states, [state])

    # Elemental state, bond.
    state = state_from_str('A_[m]--B_[n]')
    assert state.is_elemental
    assert elems_eq(state.components, [spec_from_str('A'), spec_from_str('B')])
    assert not state.is_neutral
    assert elems_eq(state.neutral_states, [state_from_str('A_[m]--0'), state_from_str('B_[n]--0')])

    # Non-elemental state, free binding domain.
    state = state_from_str('A--0')
    assert not state.is_elemental
    assert elems_eq(state.components, [spec_from_str('A')])
    assert state.is_neutral
    assert elems_eq(state.neutral_states, [state])

    # Non-elemental state, bond.
    state = state_from_str('A--B_[n]')
    assert not state.is_elemental
    assert elems_eq(state.components, [spec_from_str('A'), spec_from_str('B')])
    assert not state.is_neutral
    # @todo JCR 20160920 This is not correct, the state A--0 should be replaced by the list of domains of A that
    # @todo              can bind to B. Perhaps we should invent a new notation for this, since this property cannot
    # @todo              be known by just looking at the state itself, but requires knowledge of the rest of the system.
    assert elems_eq(state.neutral_states, [state_from_str('A--0'), state_from_str('B_[n]--0')])


def test_ppi_parsing() -> None:
    # @todo JCR 20160920 Bonds are symmetric. We need extra structure to allow this is the StateDef.
    # assert state_from_string('A_[x]--B_[y]') == state_from_string('B_[y]--A_[x]')

    # Too fine resolution (higher than elemental) raises.
    with pytest.raises(SyntaxError):
        state_from_str('A_[(x)]--B_[(y)]')


def test_ppi_superset_subset() -> None:
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


def test_ppi_structured_index() -> None:
    state = state_from_str('X@4--Z@3')
    assert spec_from_str('X@4') in state.specs
    assert spec_from_str('Z@3') in state.specs


def test_ppi_mutual_exclusivity() -> None:
    assert state_from_str('A_[x]--B_[y]').is_mutually_exclusive_with(state_from_str('A_[x]--C_[z]'))

def test_ppi_is_global() -> None:
    assert not state_from_str('A_[x]--B_[y]').is_global
###                          ###
# Test self-interaction states #
###                          ###

def test_ipi_props() -> None:
    # NOTE : States with free binding domain are the same states as tested in 'test_ppi_props', tested there.

    # Elemental state, bond.
    state = state_from_str('A_[m]--[n]')
    assert state.is_elemental
    assert state.components == [spec_from_str('A')]
    assert not state.is_neutral
    assert elems_eq(state.neutral_states, [state_from_str('A_[m]--0'), state_from_str('A_[n]--0')])

    # # Non-elemental state, bond.
    # state = state_from_str('A--A_[n]')
    # assert not state.is_elemental
    # assert elems_eq(state.components, [spec_from_str('A')])
    # assert not state.is_neutral
    # # @todo JCR 20160920 This is not correct, the state A--0 should be replaced by the list of domains of A that
    # # @todo              can bind to B. Perhaps we should invent a new notation for this, since this property cannot
    # # @todo              be known by just looking at the state itself, but requires knowledge of the rest of the system.
    # # @todo              Furthermore, the state A--0 is a superset of the other one, so it should not appear independently.
    # assert elems_eq(state.neutral_states, [state_from_str('A--0'), state_from_str('A_[n]--0')])


def test_ipi_parsing() -> None:
    # @todo JCR 20160920 Bonds are symmetric. We need extra structure to capture this is the StateDef.
    # assert state_from_string('A_[x]--B_[y]') == state_from_string('B_[y]--A_[x]')

    # Too fine resolution (higher than elemental) raises.
    with pytest.raises(SyntaxError):
        state_from_str('A_[(x)]--[(y)]')


def test_ipi_superset_subset() -> None:
    # NOTE : States with free binding domain are the same states as tested in 'test_ppi_superset_subset', tested there.

    # @todo JCR 20160920 There is a problem here: the StateDef requires the second variable in the State string to be
    # @todo              a Locus (i.e. '[n]') and not a full Spec (i.e. 'A_[n]').
    pass


def test_ipi_sorting() -> None:
    assert state_from_str('A_[x]--[y]') == state_from_str('A_[x]--[y]')
    assert state_from_str('A_[x]--[y]') == state_from_str('A_[y]--[x]')
    assert state_from_str('A_[x]--[y]') != state_from_str('B_[x]--[y]')
    assert state_from_str('A_[x]--[y]') < state_from_str('B_[x]--[y]')
    assert state_from_str('A_[a]--[b]') < state_from_str('A_[x]--[y]')

    #with pytest.raises(NotImplemented):
    #    state_from_str('A_[a]--[b]') == spec_from_str('A_[a]')


def test_ipi_get_item() -> None:
    assert spec_from_str('A_[x]') == state_from_str('A_[x]--[y]')['$x']
    assert spec_from_str('A_[y]') == state_from_str('A_[x]--[y]')['$y']

    with pytest.raises(AssertionError):
        assert spec_from_str('A_[y]') == state_from_str('A_[y]--[x]')['$x']

    with pytest.raises(KeyError):
        state_from_str('A_[y]--[x]')['$z']


def test_ipi_is_mutually_exclusive_with() -> None:
    assert not state_from_str('A_[x]--[y]').is_mutually_exclusive_with(state_from_str('A_[x]--[y]'))
    assert state_from_str('A_[x]--[y]').is_mutually_exclusive_with(state_from_str('A_[x]--B_[y]'))
    assert state_from_str('A_[x]--[y]').is_mutually_exclusive_with(state_from_str('A_[y]--0'))


def test_ipi_is_homodimer() -> None:
    assert not state_from_str('A_[x]--[y]').is_homodimer


def test_ipi_update_spec() -> None:
    state = state_from_str('A_[x]--[y]')
    state.update_specs({spec_from_str('A_[x]'):spec_from_str('A_[z]')})
    assert state == state_from_str('A_[y]--[z]')
    state.update_specs({spec_from_str('A_[y]'):spec_from_str('A_[z]')})
    assert state == state_from_str('A_[z]--[z]')

def test_ipi_is_global() -> None:
    assert not state_from_str('A_[x]--[y]').is_global

def test_ipi_is_structured() -> None:
    assert not state_from_str('A_[x]--[y]').is_structured
    assert state_from_str('A@0_[x]--[y]').is_structured

def test_ipi_to_structured_from_spec() -> None:
    homodimer = state_from_str('A_[x]--[y]')
    structured_homodimer = homodimer.to_structured_from_spec(spec_from_str('A@0'))
    assert structured_homodimer == state_from_str('A@0_[x]--[y]')

def test_ipi_to_structured_from_state() -> None:
    structured_homodimer = state_from_str('A_[x]--[y]').to_structured_from_state(state_from_str('A@0_[w]--[z]'))
    assert structured_homodimer == state_from_str('A@0_[x]--[y]')

def test_ipi_is_superset_of():

    # Happy path, subset, bond.
    assert state_from_str('A_[m]--[n]').is_subset_of(state_from_str('A_[m]--[n]'))
    assert state_from_str('A_[m]--[n]').is_superset_of(state_from_str('A_[m]--[n]'))

    # Sad path, superset, bond.
    assert not state_from_str('A_[m]--B').is_subset_of(state_from_str('A_[m]--[n]'))
    assert not state_from_str('A--A').is_subset_of(state_from_str('A_[m]--[n]'))
    assert not state_from_str('A_[m]--[n]').is_subset_of(state_from_str('A--A'))
    assert not state_from_str('A_[m]--B_[n]').is_subset_of(state_from_str('A_[m]--[n]'))


###                      ###
# Test Fully Neutral state #
###                      ###

def test_fully_neutral() -> None:
    assert state_from_str('0') == FullyNeutralState()


def test_global() -> None:
    x = state_from_str('[IN]')
    assert x.is_global


def test_properties_fully_neutral() -> None:
    fully_neutral_state = state_from_str('0')
    assert not fully_neutral_state.is_structured

    with pytest.raises(AssertionError):
        fully_neutral_state.is_subset_of(state_from_str('0'))

    with pytest.raises(AssertionError):
        fully_neutral_state.is_superset_of(state_from_str('0'))

    with pytest.raises(AssertionError):
        fully_neutral_state.is_elemental

    with pytest.raises(AssertionError):
        fully_neutral_state.is_neutral

    with pytest.raises(AssertionError):
        fully_neutral_state.neutral_states

    with pytest.raises(AssertionError):
        fully_neutral_state.to_structured_from_spec(spec_from_str('A_[m]'))

    with pytest.raises(AssertionError):
        fully_neutral_state.is_global

    with pytest.raises(AssertionError):
        fully_neutral_state.is_mutually_exclusive_with(state_from_str('0'))

    with pytest.raises(AssertionError):
        fully_neutral_state.to_structured_from_state(state_from_str('A@0_[x]--B@1_[y]'))

    with pytest.raises(AssertionError):
        fully_neutral_state.specs

    with pytest.raises(AssertionError):
        fully_neutral_state.is_homodimer


def test_sort_fully_neutral():
    assert state_from_str('0') < state_from_str('A_[x]--0')
    assert state_from_str('0') < state_from_str('A_[x]--A_[y]')
    assert state_from_str('0') < state_from_str('A_[x]--B_[y]')
    assert state_from_str('0') < state_from_str('A_[(x)]-{p}')
    assert state_from_str('0') < state_from_str('A_[(x)]-{0}')
    assert state_from_str('0') < state_from_str('[OUT]')

###                      ###
#   Test homodimer state   #
###                      ###
def test_homodi_to_structured_from_spec_non_sturcured() -> None:
    homodimer = state_from_str('A_[x]--A_[y]')
    structured_spec0 = spec_from_str('A@0')
    structured_spec1 = spec_from_str('A@1')
    structured_homodimer = homodimer.to_structured_from_spec(structured_spec0)
    structured_homodimer = structured_homodimer.to_structured_from_spec(structured_spec1)
    assert structured_homodimer == state_from_str('A@0_[x]--A@1_[y]')


def test_homodi_to_structured_from_spec_structured() -> None:
    homodimer = state_from_str('A@0_[x]--A@1_[y]')
    structured_spec0 = spec_from_str('A@2')
    structured_spec1 = spec_from_str('A@3')
    structured_homodimer = homodimer.to_structured_from_spec(structured_spec0)
    structured_homodimer = structured_homodimer.to_structured_from_spec(structured_spec1)
    assert structured_homodimer == state_from_str('A@0_[x]--A@1_[y]')


def test_homodi_to_structured_from_state() -> None:
    structured_homodimer = state_from_str('A@0_[x]--A@1_[y]')
    non_structured_diff_homodimer = state_from_str('B_[z]--B_[w]')
    non_structured_homodimer = state_from_str('A_[z]--A_[w]')

    with pytest.raises(AssertionError):
        non_structured_diff_homodimer.to_structured_from_state(structured_homodimer)

    non_structured_homodimer = non_structured_homodimer.to_structured_from_state(structured_homodimer)
    new_idx = set([spec.struct_index for spec in non_structured_homodimer.specs])
    assert len(new_idx) == 2
    assert all(idx in [0, 1] for idx in new_idx)


##### EmptyBindingState #####
def test_resolution_EmptyBindingState() -> None:
    with pytest.raises(SyntaxError):
        EmptyBindingState(spec_from_str('A_[(x)]'))


def test_ordering_EmptyBindingState() -> None:
    assert state_from_str('A_[x]--0') == state_from_str('A_[x]--0')
    assert state_from_str('A_[x]--0') < state_from_str('B_[x]--0')
    assert state_from_str('A_[x]--0') < state_from_str('A_[y]--0')

    #with pytest.raises(NotImplemented):
    #    state_from_str('A_[x]--0') == spec_from_str('A_[x]')

###                      ###
#   Test global state      #
###                      ###
def test_input_state() -> None:
    state = state_from_str('[BLA]')
    assert isinstance(state, GlobalState)

def test_global_state_propterties() -> None:
    state = state_from_str('[BLA]')
    assert state.is_global
    assert state.is_elemental
    assert not all(state.is_mutually_exclusive_with(state_to_test)for state_to_test in [state_from_str('[BLUB]'),
                                                                                        state_from_str('A_[x]--B_[y]'),
                                                                                        state_from_str('A_[x]--0'),
                                                                                        state_from_str('A_[(r)]-{p}'),
                                                                                        state_from_str('A_[(r)]-{0}')] )
    assert not state.is_homodimer
    assert state.is_structured
    assert state.neutral_states == []

def test_global_get_item():
    with pytest.raises(KeyError):
        state_from_str('[BLA]')['$y']

    assert state_from_str('[BLA]') > state_from_str('X_[x]--0')
    assert state_from_str('[BLA]') > state_from_str('X_[x]--A_[y]')
    assert state_from_str('[BLA]') > state_from_str('X_[x]--B_[y]')
    assert state_from_str('[BLA]') > state_from_str('X_[(x)]-{p}')
    assert state_from_str('[BLA]') > state_from_str('X_[(x)]-{0}')
    assert state_from_str('[BLA]') < state_from_str('[OUT]')