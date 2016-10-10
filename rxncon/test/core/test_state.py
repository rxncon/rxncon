from rxncon.core.state import state_from_str, FullyNeutralState
from rxncon.core.spec import bond_spec_from_str, mol_spec_from_str
from rxncon.util.utils import elems_eq
import pytest

###                      ###
# Test modification states #
###                      ###

def test_modification_props():
    # Elemental state, neutral.
    state = state_from_str('A_[(res)]-{0}')
    assert state.is_elemental
    assert elems_eq(state.components, [mol_spec_from_str('A')])
    assert state.target == mol_spec_from_str('A_[(res)]')
    assert state.is_neutral
    assert elems_eq(state.neutral_states, [state])

    # Elemental state, non-neutral.
    state = state_from_str('A_[(res)]-{p}')
    assert state.is_elemental
    assert elems_eq(state.components, [mol_spec_from_str('A')])
    assert state.target == mol_spec_from_str('A_[(res)]')
    assert not state.is_neutral
    assert elems_eq(state.neutral_states, [state_from_str('A_[(res)]-{0}')])

    # Non-elemental state, neutral.
    state = state_from_str('A_[dom]-{0}')
    assert not state.is_elemental
    assert elems_eq(state.components, [mol_spec_from_str('A')])
    assert state.target == mol_spec_from_str('A_[dom]')
    assert state.is_neutral
    assert elems_eq(state.neutral_states, [state])

    # Non-elemental state, non-neutral.
    state = state_from_str('A_[dom]-{p}')
    assert not state.is_elemental
    assert elems_eq(state.components, [mol_spec_from_str('A')])
    assert state.target == mol_spec_from_str('A_[dom]')
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
    assert elems_eq(state.components, [mol_spec_from_str('A')])
    assert state.target == mol_spec_from_str('A_[m]')
    assert state.is_neutral
    assert elems_eq(state.neutral_states, [state])

    # Elemental state, bond.
    state = state_from_str('A_[m]--B_[n]')
    assert state.is_elemental
    assert elems_eq(state.components, [mol_spec_from_str('A'), mol_spec_from_str('B')])
    assert state.target == bond_spec_from_str('A_[m]~B_[n]')
    assert not state.is_neutral
    assert elems_eq(state.neutral_states, [state_from_str('A_[m]--0'), state_from_str('B_[n]--0')])

    # Non-elemental state, free binding domain.
    state = state_from_str('A--0')
    assert not state.is_elemental
    assert elems_eq(state.components, [mol_spec_from_str('A')])
    assert state.target == mol_spec_from_str('A')
    assert state.is_neutral
    assert elems_eq(state.neutral_states, [state])

    # Non-elemental state, bond.
    state = state_from_str('A--B_[n]')
    assert not state.is_elemental
    assert elems_eq(state.components, [mol_spec_from_str('A'), mol_spec_from_str('B')])
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
    assert elems_eq(state.components, [mol_spec_from_str('A')])
    assert state.target == bond_spec_from_str('A_[m]~A_[n]')
    assert not state.is_neutral
    assert elems_eq(state.neutral_states, [state_from_str('A_[m]--0'), state_from_str('A_[n]--0')])

    # Non-elemental state, bond.
    state = state_from_str('A--A_[n]')
    assert not state.is_elemental
    assert elems_eq(state.components, [mol_spec_from_str('A')])
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



# @pytest.fixture
# def the_case_hierarchy():
#     return [
#         # HierarchyTestCase(sta.state_from_string('B'),
#         #                   [sta.state_from_string('B'),
#         #                    sta.state_from_string('A--B'), sta.state_from_string('A--B_[d]'),
#         #                    sta.state_from_string('A--B_[d/s]'),
#         #                    sta.state_from_string('A--B_[d/s(r)]')
#         #                    ]),
#         HierarchyTestCase(sta.state_from_string('A-{p}'),
#                           [sta.state_from_string('A_[n]-{p}')]),
# >>>>>>> regulatory-graph
#
# def test_empty_ppi_states():
#     state = state_from_string('A_[x]--0')
#     assert state.target == mol_spec_from_string('A_[x]')
#     assert state.is_elemental
#
#     state = state_from_string('0--B_[x]')
#     assert state.target == mol_spec_from_string('B_[x]')
#     assert state.is_elemental
#
#
# def test_ipi_states():
#     state = state_from_string('A_[n]--[m]')
#     assert state.target == bond_spec_from_string('A_[m]~A_[n]')
#     assert state.is_elemental
#
#
#
# def test_mod_states():
#     state = state_from_string('A_[(r)]-{0}')
#     assert state.target == mol_spec_from_string('A_[(r)]')
#     assert state.is_elemental
#     assert state.is_subset_of(state_from_string('A-{0}'))
#
#
# def test_neutral_states():
#     state = state_from_string('A_[(r)]-{p}')
#     assert not state.is_neutral
#
#     neutral = state_from_string('A_[(r)]-{0}')
#     assert neutral.is_neutral
#     assert state.neutral_states == [neutral]
#
#     state = state_from_string('A_[x]--B_[y]')
#     assert not state.is_neutral
#     first_neutral = state_from_string('A_[x]--0')
#     second_neutral = state_from_string('B_[y]--0')
#     assert first_neutral.is_neutral
#     assert second_neutral.is_neutral
#             # HierarchyTestCase(sta.state_from_string('A@0'),
#             #                   [sta.state_from_string('A@0'),
#             #                    sta.state_from_string('A@0--B@1'), sta.state_from_string('A@0_[d]--B@1'), sta.state_from_string('A@0_[d/s]--B@1'),
#             #                    sta.state_from_string('A@0_[d/s(r)]--B@1'),
#             #                    sta.state_from_string('A@0-{p}'), sta.state_from_string('A@0_[d]-{p}'), sta.state_from_string('A@0_[d/s]-{p}'),
#             #                    sta.state_from_string('A@0_[d/s(r)]-{p}'), sta.state_from_string('A@0_[(r)]-{p}')
#             #                    ]),
#             #todo: how do we handle Input states with @
#             #HierarchyTestCase(sta.state_from_string('[INPUT]'),
#             #                  [sta.state_from_string('[INPUT]')])
#     ]
#
#
# def test_not_super_or_subspecification(the_case_no_hierarchy):
#     for test_case in the_case_no_hierarchy:
#         for state in test_case.superset_of:
#             assert not test_case.state.is_superset_of(state)
#             assert not test_case.state.is_subset_of(state)
#             assert not state.is_superset_of(test_case.state)
#             assert not state.is_subset_of(test_case.state)
#
#
# @pytest.fixture
# def the_case_no_hierarchy():
#     return [
#         HierarchyTestCase(sta.state_from_string('A@0_[n]--B@1'),
#                           [sta.state_from_string('A@0--B@1_[m]')]),
#
#         HierarchyTestCase(sta.state_from_string('A_[m]--B'),
#                           [sta.state_from_string('A_[m]--0')]),
#
#         HierarchyTestCase(sta.state_from_string('A_[m]--A_[m]'),
#                           [sta.state_from_string('A_[m]--0')]),
#
#         HierarchyTestCase(sta.state_from_string('A@0_[n]--B@1'),
#                           [sta.state_from_string('A@0_[n]-{P}')]),
#
#         HierarchyTestCase(sta.state_from_string('A@0--B@1'),
#                           [sta.state_from_string('A@0-{P}')]),
#
#         HierarchyTestCase(sta.state_from_string('A@0_[d1]-{P}'),
#                           [sta.state_from_string('A@0_[d2]-{P}')]),
#
#         HierarchyTestCase(sta.state_from_string('A@0_[d/s1]-{P}'),
#                           [sta.state_from_string('A@0_[d/s2]-{P}')]),
#
#         HierarchyTestCase(sta.state_from_string('A@0_[d/s(r1)]-{P}'),
#                           [sta.state_from_string('A@0_[d/s(r2)]-{P}')]),
#
#         HierarchyTestCase(sta.state_from_string('A@0_[d1]--B@1'),
#                           [sta.state_from_string('A@0_[d2]--B@1')]),
#
#         HierarchyTestCase(sta.state_from_string('A@0_[d/s1]--B@1'),
#                           [sta.state_from_string('A@0_[d/s2]--B@1')]),
#
#         HierarchyTestCase(sta.state_from_string('A@0_[d/s(r1)]--B@1'),
#                           [sta.state_from_string('A@0_[d/s(r2)]--B@1')]),
#
#         HierarchyTestCase(sta.state_from_string('A@0-{p}'),
#                           [sta.state_from_string('A@0_[(n)]-{ub}')])
#
#         # HierarchyTestCase(sta.state_from_string('A@0'),
#         #                   [sta.state_from_string('B@1')]),
#
#         #HierarchyTestCase(sta.state_from_string('[INPUT]'),
#         #                  [sta.state_from_string('B@1'), sta.state_from_string('A@0--B@1'),
#         #                   sta.state_from_string('A@0-{P}')])
#     ]
#
#
# # StateTestCase = namedtuple('StateTestCase', ['state', 'expected_specifications', 'expected_string'])
# #
# # @pytest.fixture
# # def the_case_state():
# #     return [
# #
# #         StateTestCase(sta.state_from_string('Agene'),
# #                       [spec.DnaSpecification('A', None, spec.Domain(None, None, None))],
# #                       'Agene'),
# #         StateTestCase(sta.state_from_string('Agene_[d/s(r)]'),
# #                       [spec.DnaSpecification('A', None, spec.Domain('d', 's', 'r'))],
# #                       'Agene_[d/s(r)]'),
# #         StateTestCase(sta.state_from_string('AmRNA'),
# #                       [spec.RnaSpecification('A', None, spec.Domain(None, None, None))],
# #                       'AmRNA'),
# #         StateTestCase(sta.state_from_string('AmRNA_[d/s(r)]'),
# #                       [spec.RnaSpecification('A', None, spec.Domain('d', 's', 'r'))],
# #                       'AmRNA_[d/s(r)]'),
# #         StateTestCase(sta.state_from_string('A'),
# #                       [spec.ProteinSpecification('A', None, spec.Domain(None, None, None))],
# #                       'A'),
# #         StateTestCase(sta.state_from_string('A_[d/s(r)]'),
# #                       [spec.ProteinSpecification('A', None, spec.Domain('d', 's', 'r'))],
# #                       'A_[d/s(r)]')
# #     ]
# #
# #
# # @pytest.fixture
# # def the_case_state_structured():
# #     return [
# #         StateTestCase(sta.state_from_string('Agene@0'),
# #                       [spec.DnaSpecification('A', 0, spec.Domain(None, None, None))],
# #                       'Agene@0'),
# #         StateTestCase(sta.state_from_string('Agene@0_[d/s(r)]'),
# #                       [spec.DnaSpecification('A', 0, spec.Domain('d', 's', 'r'))],
# #                       'Agene@0_[d/s(r)]'),
# #         StateTestCase(sta.state_from_string('AmRNA@0'),
# #                       [spec.RnaSpecification('A', 0, spec.Domain(None, None, None))],
# #                       'AmRNA@0'),
# #         StateTestCase(sta.state_from_string('AmRNA@0_[d/s(r)]'),
# #                       [spec.RnaSpecification('A', 0, spec.Domain('d', 's', 'r'))],
# #                       'AmRNA@0_[d/s(r)]'),
# #         StateTestCase(sta.state_from_string('A@0'),
# #                       [spec.ProteinSpecification('A', 0, spec.Domain(None, None, None))],
# #                       'A@0'),
# #         StateTestCase(sta.state_from_string('A@0_[d/s(r)]'),
# #                       [spec.ProteinSpecification('A', 0, spec.Domain('d', 's', 'r'))],
# #                       'A@0_[d/s(r)]'),
# #
# #     ]
#
#
# # def test_state_building(the_case_state, the_case_state_structured):
# #     for the_case in the_case_state + the_case_state_structured:
# #         is_state_correct(the_case)
#
#
# # def is_state_correct(the_case):
# #     assert all(the_case.state.variables[variable] in the_case.expected_specifications for variable in the_case.state.variables)
# #     assert str(the_case.state) == the_case.expected_string
#
# ResolutionTestCase = namedtuple('ResolutionTestCase', ['state', 'expected_resolution', 'expected_neutral_modifier', 'expected_is_elemental'])
#
# @pytest.fixture
# def test_case_resolution_and_neutral_modifier():
#     return [
#         # ResolutionTestCase(sta.state_from_string('A'),
#         #                    [(spec.ProteinSpecification('A', None, spec.Domain(None, None, None)),
#         #                      spec.SpecificationResolution.component)],
#         #                    None,
#         #                    True),
#         ResolutionTestCase(sta.state_from_string('A_[m]--[n]'),
#                            [(spec.ProteinSpecification('A', None, spec.Domain('m', None, None)),
#                              spec.SpecificationResolution.domain)],
#                            None,
#                            True),
#
#         ResolutionTestCase(sta.state_from_string('A_[(r)]-{p}'),
#                            [(spec.ProteinSpecification('A', None, spec.Domain(None, None, 'r')),
#                              spec.SpecificationResolution.residue)],
#                            sta.StateModifier.neutral,
#                            True),
#
#         ResolutionTestCase(sta.state_from_string('A_[n]-{p}'),
#                            [(spec.ProteinSpecification('A', None, spec.Domain('n', None, None)),
#                              spec.SpecificationResolution.residue)],
#                            sta.StateModifier.neutral,
#                            False),
#     ]
# @pytest.fixture
# def test_case_resolution_and_neutral_modifier_structured():
#     return [
#
#
#         # ResolutionTestCase(sta.state_from_string('A@0'),
#         #                    [(spec.ProteinSpecification('A', 0, spec.Domain(None, None, None)),
#         #                     spec.SpecificationResolution.component)],
#         #                    None,
#         #                    True),
#
#         # ResolutionTestCase(sta.state_from_string('Agene@0'),
#         #                    [(spec.DnaSpecification('A', 0, spec.Domain(None, None, None)),
#         #                      spec.SpecificationResolution.component)],
#         #                    None,
#         #                    True),
#         #
#         # ResolutionTestCase(sta.state_from_string('AmRNA@0'),
#         #                    [(spec.RnaSpecification('A', 0, spec.Domain(None, None, None)),
#         #                      spec.SpecificationResolution.component)],
#         #                    None,
#         #                    True),
#
#         #todo: discuss this What is with the DomainResolution ?
#         ResolutionTestCase(sta.state_from_string('A@0_[m]--[n]'),
#                            [(spec.ProteinSpecification('A', 0, spec.Domain('m', None, None)),
#                              spec.SpecificationResolution.domain)],
#                            None,
#                            True),
#
#         ResolutionTestCase(sta.state_from_string('A@0_[(r)]-{p}'),
#                            [(spec.ProteinSpecification('A', 0, spec.Domain(None, None, 'r')),
#                              spec.SpecificationResolution.residue)],
#                            sta.StateModifier.neutral,
#                            True),
#
#         ResolutionTestCase(sta.state_from_string('A@0_[d(r)]-{p}'),
#                            [(spec.ProteinSpecification('A', 0, spec.Domain('d', None, 'r')),
#                              spec.SpecificationResolution.residue)],
#                            sta.StateModifier.neutral,
#                            True),
#
#         ResolutionTestCase(sta.state_from_string('A@0_[d/s(r)]-{p}'),
#                            [(spec.ProteinSpecification('A', 0, spec.Domain('d', 's', 'r')),
#                              spec.SpecificationResolution.residue)],
#                            sta.StateModifier.neutral,
#                            True),
#
#         ResolutionTestCase(sta.state_from_string('Agene@0_[d/s(r)]-{p}'),
#                            [(spec.DnaSpecification('A', 0, spec.Domain('d', 's', 'r')),
#                              spec.SpecificationResolution.residue)],
#                            sta.StateModifier.neutral,
#                            True),
#
#         ResolutionTestCase(sta.state_from_string('AmRNA@0_[d/s(r)]-{p}'),
#                            [(spec.RnaSpecification('A', 0, spec.Domain('d', 's', 'r')),
#                              spec.SpecificationResolution.residue)],
#                            sta.StateModifier.neutral,
#                            True),
#
#         ResolutionTestCase(sta.state_from_string('A@0_[n]-{p}'),
#                            [(spec.ProteinSpecification('A', 0, spec.Domain('n', None, None)),
#                              spec.SpecificationResolution.residue)],
#                            sta.StateModifier.neutral,
#                            False),
#
#         ResolutionTestCase(sta.state_from_string('A@0_[d]--B@1_[d]'),
#                            [(spec.ProteinSpecification('A', 0, spec.Domain('d', None, None)),
#                              spec.SpecificationResolution.domain),
#                             (spec.ProteinSpecification('B', 1, spec.Domain('d', None, None)),
#                              spec.SpecificationResolution.domain)
#                             ],
#                            None,
#                            True),
#
#         ResolutionTestCase(sta.state_from_string('A@0_[d/s]--B@1_[d/s]'),
#                            [(spec.ProteinSpecification('A', 0, spec.Domain('d', 's', None)),
#                              spec.SpecificationResolution.domain),
#                             (spec.ProteinSpecification('B', 1, spec.Domain('d', 's', None)),
#                              spec.SpecificationResolution.domain)
#                             ],
#                            None,
#                            False),
#
#         ResolutionTestCase(sta.state_from_string('A@0--B@1'),
#                            [(spec.ProteinSpecification('A', 0, spec.Domain(None, None, None)),
#                              spec.SpecificationResolution.domain),
#                             (spec.ProteinSpecification('B', 1, spec.Domain(None, None, None)),
#                              spec.SpecificationResolution.domain)
#                             ],
#                            None,
#                            False)
#     ]
#
# def test_resolution_and_default_modifier(test_case_resolution_and_neutral_modifier, test_case_resolution_and_neutral_modifier_structured):
#     for the_case in test_case_resolution_and_neutral_modifier + test_case_resolution_and_neutral_modifier_structured:
#         is_resolution_and_neutral_modifier_correct(the_case)
#
# def is_resolution_and_neutral_modifier_correct(the_case):
#     for component in the_case.state.components():
#         assert (component, the_case.state.elemental_resolution(component)) in the_case.expected_resolution
#     assert the_case.state.neutral_modifier == the_case.expected_neutral_modifier
#     assert the_case.state.is_elemental() == the_case.expected_is_elemental
#     assert first_neutral in state.neutral_states
#     assert second_neutral in state.neutral_states
#
#     assert first_neutral.neutral_states == [first_neutral]
#     assert second_neutral.neutral_states == [second_neutral]
#
# def test_test():
#     a = sta.state_from_string('A-{p}')
#     str(a)
#     a = sta.state_from_string('A_[m]--B')
#     b = sta.state_from_string('A_[m]--0')
#     a.is_subset_of(b)
#     b.is_subset_of(a)
#     #sta.state_from_string('A_[m]--B').is_subset_of(sta.state_from_string('A_[m]--0'))
#     #sta.state_from_string('A_[m]--0').is_subset_of(sta.state_from_string('A_[m]--B'))
