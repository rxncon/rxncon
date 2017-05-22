from rxncon.core.state import state_from_str, FullyNeutralState, GlobalState, EmptyBindingState, initialize_state_modifiers
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
    assert elems_eq(state.neutral_states, [state_from_str('A--0'), state_from_str('B_[n]--0')])


def test_ppi_parsing() -> None:
    assert state_from_str('A_[x]--B_[y]') == state_from_str('B_[y]--A_[x]')

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
    assert not state_from_str('A_[x]--B_[y]').is_mutually_exclusive_with(state_from_str('A_[y]--C_[z]'))


def test_ppi_is_global() -> None:
    assert not state_from_str('A_[x]--B_[y]').is_global


###                          ###
# Test self-interaction states #
###                          ###

def test_ipi_props() -> None:
    # Elemental state, bond.
    state = state_from_str('A_[m]--[n]')
    assert state.is_elemental
    assert state.components == [spec_from_str('A')]
    assert not state.is_neutral
    assert elems_eq(state.neutral_states, [state_from_str('A_[m]--0'), state_from_str('A_[n]--0')])


def test_ipi_parsing() -> None:
    assert state_from_str('A_[x]--B_[y]') == state_from_str('B_[y]--A_[x]')

    # Too fine resolution (higher than elemental) raises.
    with pytest.raises(SyntaxError):
        state_from_str('A_[(x)]--[(y)]')

    with pytest.raises(AssertionError):
        state = state_from_str('A@1_[x]--[y]')
        state.update_specs({
            spec_from_str('A@1_[x]'): spec_from_str('A@5_[x]'),
        })

def test_ipi_superset_subset() -> None:
    # Happy path, subset, bond.
    assert state_from_str('A_[m]--[n]').is_subset_of(state_from_str('A_[m]--[n]'))
    assert state_from_str('A_[m]--[n]').is_superset_of(state_from_str('A_[m]--[n]'))

    # Sad path, superset, bond.
    assert not state_from_str('A_[m]--B').is_subset_of(state_from_str('A_[m]--[n]'))
    assert not state_from_str('A--A').is_subset_of(state_from_str('A_[m]--[n]'))
    assert not state_from_str('A_[m]--[n]').is_subset_of(state_from_str('A--A'))
    assert not state_from_str('A_[m]--B_[n]').is_subset_of(state_from_str('A_[m]--[n]'))


def test_ipi_sorting() -> None:
    assert state_from_str('A_[x]--[y]') == state_from_str('A_[x]--[y]')
    assert state_from_str('A_[x]--[y]') == state_from_str('A_[y]--[x]')
    assert state_from_str('A_[x]--[y]') != state_from_str('B_[x]--[y]')
    assert state_from_str('A_[x]--[y]') < state_from_str('B_[x]--[y]')
    assert state_from_str('A_[a]--[b]') < state_from_str('A_[x]--[y]')


def test_ipi_is_mutually_exclusive_with() -> None:
    assert not state_from_str('A_[x]--[y]').is_mutually_exclusive_with(state_from_str('A_[x]--[y]'))
    assert state_from_str('A_[x]--[y]').is_mutually_exclusive_with(state_from_str('A_[x]--B_[y]'))
    assert state_from_str('A_[x]--[y]').is_mutually_exclusive_with(state_from_str('A_[y]--0'))


def test_ipi_is_homodimer() -> None:
    assert not state_from_str('A_[x]--[y]').is_homodimer


def test_ipi_update_spec() -> None:
    state = state_from_str('A_[x]--[y]')
    state.update_specs({spec_from_str('A_[x]'): spec_from_str('A_[z]')})
    assert state == state_from_str('A_[y]--[z]')
    state.update_specs({spec_from_str('A_[y]'): spec_from_str('A_[z]')})
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

    with pytest.raises(NotImplementedError):
        fully_neutral_state.is_subset_of(state_from_str('0'))

    with pytest.raises(NotImplementedError):
        fully_neutral_state.is_superset_of(state_from_str('0'))

    with pytest.raises(NotImplementedError):
        fully_neutral_state.is_elemental

    with pytest.raises(NotImplementedError):
        fully_neutral_state.is_neutral

    with pytest.raises(NotImplementedError):
        fully_neutral_state.neutral_states

    with pytest.raises(NotImplementedError):
        fully_neutral_state.to_structured_from_spec(spec_from_str('A_[m]'))

    with pytest.raises(NotImplementedError):
        fully_neutral_state.is_global

    with pytest.raises(NotImplementedError):
        fully_neutral_state.is_mutually_exclusive_with(state_from_str('0'))

    with pytest.raises(NotImplementedError):
        fully_neutral_state.to_structured_from_state(state_from_str('A@0_[x]--B@1_[y]'))

    with pytest.raises(NotImplementedError):
        fully_neutral_state.specs

    with pytest.raises(NotImplementedError):
        fully_neutral_state.is_homodimer


def test_sort_fully_neutral() -> None:
    assert state_from_str('0') < state_from_str('A_[x]--0')
    assert state_from_str('0') < state_from_str('A_[x]--A_[y]')
    assert state_from_str('0') < state_from_str('A_[x]--B_[y]')
    assert state_from_str('0') < state_from_str('A_[(x)]-{p}')
    assert state_from_str('0') < state_from_str('A_[(x)]-{0}')
    assert state_from_str('0') < state_from_str('[OUT]')


###                      ###
#   Test homodimer state   #
###                      ###
def test_homodimer_same_struct_index() -> None:
    with pytest.raises(AssertionError):
        state_from_str('A@4_[x]--A@4_[x]')


def test_homodimer_to_structured_from_spec_non_structured() -> None:
    structured_homodimer = state_from_str('A_[x]--A_[y]')\
        .to_structured_from_spec(spec_from_str('A@0'))\
        .to_structured_from_spec(spec_from_str('A@1'))

    assert structured_homodimer == state_from_str('A@0_[x]--A@1_[y]')


def test_homodimer_to_structured_from_spec_previously_structured() -> None:
    # A previously structured homodimer should remain invariant upon updating its specs.
    structured_homodimer = state_from_str('A@0_[x]--A@1_[y]')\
        .to_structured_from_spec(spec_from_str('A@2'))\
        .to_structured_from_spec(spec_from_str('A@3'))

    assert structured_homodimer == state_from_str('A@0_[x]--A@1_[y]')


def test_homodimer_to_structured_from_state_overlapping() -> None:
    structured_homodimer = state_from_str('A_[z]--A_[w]').to_structured_from_state(state_from_str('A@0_[x]--A@1_[y]'))
    assert structured_homodimer == state_from_str('A@0_[w]--A@1_[z]')


def test_homodimer_to_structured_from_state_non_overlapping() -> None:
    with pytest.raises(AssertionError):
        state_from_str('B_[z]--B_[w]').to_structured_from_state(state_from_str('A@0_[x]--A@1_[y]'))


def test_homodimer_update_specs() -> None:
    structured_homodimer = state_from_str('A@0_[x]--A@1_[y]')

    structured_homodimer.update_specs({
        spec_from_str('A@0_[x]'): spec_from_str('A@5_[x]'),
        spec_from_str('A@1_[y]'): spec_from_str('A@7_[y]')
    })

    assert structured_homodimer == state_from_str('A@5_[x]--A@7_[y]')


def test_homodimer_update_specs_flip() -> None:
    structured_homodimer = state_from_str('A@0_[x]--A@1_[y]')

    structured_homodimer.update_specs({
        spec_from_str('A@0_[x]'): spec_from_str('A@1_[x]'),
        spec_from_str('A@1_[y]'): spec_from_str('A@0_[y]')
    })

    assert structured_homodimer == state_from_str('A@1_[x]--A@0_[y]')


###                       ###
#     EmptyBindingState     #
###                       ###
def test_resolution_EmptyBindingState() -> None:
    with pytest.raises(SyntaxError):
        state_from_str('A_[(x)]--0')

    with pytest.raises(SyntaxError):
        EmptyBindingState(spec_from_str('A_[(x)]'))


def test_ordering_EmptyBindingState() -> None:
    assert state_from_str('A_[x]--0') == state_from_str('A_[x]--0')
    assert state_from_str('A_[x]--0') < state_from_str('B_[x]--0')
    assert state_from_str('A_[x]--0') < state_from_str('A_[y]--0')


###                      ###
#   Test global state      #
###                      ###
def test_input_state() -> None:
    assert isinstance(state_from_str('[BLA]'), GlobalState)
    assert isinstance(state_from_str('[BLA123]'), GlobalState)


def test_global_state_properties() -> None:
    state = state_from_str('[BLA]')
    assert state.is_global
    assert state.is_elemental
    assert not any(state.is_mutually_exclusive_with(state_to_test) for state_to_test in [
        state_from_str('[BLUB]'),
        state_from_str('A_[x]--B_[y]'),
        state_from_str('A_[x]--0'),
        state_from_str('A_[(r)]-{p}'),
        state_from_str('A_[(r)]-{0}')
    ])
    assert not state.is_homodimer
    assert state.is_structured
    assert state.neutral_states == []


def test_global_ordering() -> None:
    assert state_from_str('[BLA]') > state_from_str('X_[x]--0')
    assert state_from_str('[BLA]') > state_from_str('X_[x]--A_[y]')
    assert state_from_str('[BLA]') > state_from_str('X_[x]--B_[y]')
    assert state_from_str('[BLA]') > state_from_str('X_[(x)]-{p}')
    assert state_from_str('[BLA]') > state_from_str('X_[(x)]-{0}')
    assert state_from_str('[BLA]') < state_from_str('[OUT]')

###                       ###
#   Test dynamical states   #
###                       ###

def test_dynamical_modification_state() -> None:
    with pytest.raises(ValueError):
        state_from_str('A_[(r)]-{bla}')

    initialize_state_modifiers({'bla': 'bla'})

    assert str(state_from_str('A_[(r)]-{bla}')) == 'A_[(r)]-{bla}'

    initialize_state_modifiers()

    with pytest.raises(ValueError):
        state_from_str('A_[(r)]-{bla}')
