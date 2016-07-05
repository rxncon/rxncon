import pytest
import rxncon.core.state as sta
import rxncon.core.specification as spec
from collections import namedtuple


HierarchyTestCase = namedtuple('HierarchyTestCase', ['state', 'superset_of'])


def test_hierarchy(the_case_hierarchy):
    for test_case in the_case_hierarchy:
        is_hierarchy_correct(test_case)


def is_hierarchy_correct(test_case):
    assert test_case.superset_of
    for state in test_case.superset_of:
        assert test_case.state.is_superset_of(state)
        if test_case.state != state:
            assert not test_case.state.is_subset_of(state) and test_case.state.is_subset_of(state) is not None
        assert state.is_subset_of(test_case.state)


@pytest.fixture
def the_case_hierarchy():
    return [
            HierarchyTestCase(sta.state_from_string('A@0-{p}'),
                             [sta.state_from_string('A@0_[n]-{p}')]),

            HierarchyTestCase(sta.state_from_string('A@0--B@1'),
                              [sta.state_from_string('A@0_[n]--B@1'), sta.state_from_string('A@0--B@1_[m]'), sta.state_from_string('A@0_[n]--B@1_[m]')],
                              ),

            HierarchyTestCase(sta.state_from_string('A@0_[n]--B@1'),
                              [sta.state_from_string('A@0_[n]--B@1_[m]')]),

            HierarchyTestCase(sta.state_from_string('A@0--B@1_[m]'),
                              [sta.state_from_string('A@0_[n]--B@1_[m]')]),

            HierarchyTestCase(sta.state_from_string('A@0'),
                              [sta.state_from_string('A@0'),
                               sta.state_from_string('A@0--B@1'), sta.state_from_string('A@0_[d]--B@1'), sta.state_from_string('A@0_[d/s]--B@1'),
                               sta.state_from_string('A@0_[d/s(r)]--B@1'),
                               sta.state_from_string('A@0-{p}'), sta.state_from_string('A@0_[d]-{p}'), sta.state_from_string('A@0_[d/s]-{p}'),
                               sta.state_from_string('A@0_[d/s(r)]-{p}'), sta.state_from_string('A@0_[(r)]-{p}')
                               ]),

            HierarchyTestCase(sta.state_from_string('B@1'),
                              [sta.state_from_string('B@1'),
                               sta.state_from_string('A@0--B@1'), sta.state_from_string('A@0--B@1_[d]'), sta.state_from_string('A@0--B@1_[d/s]'),
                               sta.state_from_string('A@0--B@1_[d/s(r)]')
                               ]),
            #todo: how do we handle Input states with @
            #HierarchyTestCase(sta.state_from_string('[INPUT]'),
            #                  [sta.state_from_string('[INPUT]')])
    ]


def test_not_super_or_subspecification(the_case_no_hierarchy):
    for test_case in the_case_no_hierarchy:
        for state in test_case.superset_of:
            assert not test_case.state.is_superset_of(state)
            assert not test_case.state.is_subset_of(state)
            assert not state.is_superset_of(test_case.state)
            assert not state.is_subset_of(test_case.state)


@pytest.fixture
def the_case_no_hierarchy():
    return [
        HierarchyTestCase(sta.state_from_string('A@0_[n]--B@1'),
                          [sta.state_from_string('A@0--B@1_[m]')]),

        HierarchyTestCase(sta.state_from_string('A@0_[n]--B@1'),
                          [sta.state_from_string('A@0_[n]-{P}')]),

        HierarchyTestCase(sta.state_from_string('A@0--B@1'),
                          [sta.state_from_string('A@0-{P}')]),

        HierarchyTestCase(sta.state_from_string('A@0_[d1]-{P}'),
                          [sta.state_from_string('A@0_[d2]-{P}')]),

        HierarchyTestCase(sta.state_from_string('A@0_[d/s1]-{P}'),
                          [sta.state_from_string('A@0_[d/s2]-{P}')]),

        HierarchyTestCase(sta.state_from_string('A@0_[d/s(r1)]-{P}'),
                          [sta.state_from_string('A@0_[d/s(r2)]-{P}')]),

        HierarchyTestCase(sta.state_from_string('A@0_[d1]--B@1'),
                          [sta.state_from_string('A@0_[d2]--B@1')]),

        HierarchyTestCase(sta.state_from_string('A@0_[d/s1]--B@1'),
                          [sta.state_from_string('A@0_[d/s2]--B@1')]),

        HierarchyTestCase(sta.state_from_string('A@0_[d/s(r1)]--B@1'),
                          [sta.state_from_string('A@0_[d/s(r2)]--B@1')]),

        HierarchyTestCase(sta.state_from_string('A@0'),
                          [sta.state_from_string('B@1')]),

        #HierarchyTestCase(sta.state_from_string('[INPUT]'),
        #                  [sta.state_from_string('B@1'), sta.state_from_string('A@0--B@1'),
        #                   sta.state_from_string('A@0-{P}')])
    ]


StateTestCase = namedtuple('StateTestCase', ['state', 'expected_specifications', 'expected_string'])


@pytest.fixture
def the_case_state():
    return [
        StateTestCase(sta.state_from_string('Agene@0'),
                      [spec.DnaSpecification('A', 0, spec.DomainDefinition(None, None, None))],
                      'Agene@0'),
        StateTestCase(sta.state_from_string('Agene@0_[d/s(r)]'),
                      [spec.DnaSpecification('A', 0, spec.DomainDefinition('d', 's', 'r'))],
                      'Agene@0_[d/s(r)]'),
        StateTestCase(sta.state_from_string('AmRNA@0'),
                      [spec.RnaSpecification('A', 0, spec.DomainDefinition(None, None, None))],
                      'AmRNA@0'),
        StateTestCase(sta.state_from_string('AmRNA@0_[d/s(r)]'),
                      [spec.RnaSpecification('A', 0, spec.DomainDefinition('d', 's', 'r'))],
                      'AmRNA@0_[d/s(r)]'),
        StateTestCase(sta.state_from_string('A@0'),
                      [spec.ProteinSpecification('A', 0, spec.DomainDefinition(None, None, None))],
                      'A@0'),
        StateTestCase(sta.state_from_string('A@0_[d/s(r)]'),
                      [spec.ProteinSpecification('A', 0, spec.DomainDefinition('d', 's', 'r'))],
                      'A@0_[d/s(r)]')
    ]


def test_state_building(the_case_state):
    for the_case in the_case_state:
        is_state_correct(the_case)


def is_state_correct(the_case):
    assert all(the_case.state.variables[variable] in the_case.expected_specifications for variable in the_case.state.variables)
    assert str(the_case.state) == the_case.expected_string

ResolutionTestCase = namedtuple('ResolutionTestCase', ['state', 'expected_resolution', 'expected_neutral_modifier', 'expected_is_elemental'])


@pytest.fixture
def test_case_resolution_and_neutral_modifier():
    return [
        ResolutionTestCase(sta.state_from_string('A@0'),
                           [(spec.ProteinSpecification('A', 0, spec.DomainDefinition(None, None, None)),
                            spec.SpecificationResolution.component)],
                           None,
                           True),

        ResolutionTestCase(sta.state_from_string('Agene@0'),
                           [(spec.DnaSpecification('A', 0, spec.DomainDefinition(None, None, None)),
                             spec.SpecificationResolution.component)],
                           None,
                           True),

        ResolutionTestCase(sta.state_from_string('AmRNA@0'),
                           [(spec.RnaSpecification('A', 0, spec.DomainDefinition(None, None, None)),
                             spec.SpecificationResolution.component)],
                           None,
                           True),

        #todo: discuss this What is with the DomainResolution ?
        ResolutionTestCase(sta.state_from_string('A@0_[m]--[n]'),
                           [(spec.ProteinSpecification('A', 0, spec.DomainDefinition('m', None, None)),
                             spec.SpecificationResolution.domain)],
                           None,
                           True),
        ResolutionTestCase(sta.state_from_string('A@0_[(r)]-{p}'),
                           [(spec.ProteinSpecification('A', 0, spec.DomainDefinition(None, None, 'r')),
                             spec.SpecificationResolution.residue)],
                           sta.StateModifier.neutral,
                           True),

        ResolutionTestCase(sta.state_from_string('A@0_[d(r)]-{p}'),
                           [(spec.ProteinSpecification('A', 0, spec.DomainDefinition('d', None, 'r')),
                             spec.SpecificationResolution.residue)],
                           sta.StateModifier.neutral,
                           True),

        ResolutionTestCase(sta.state_from_string('A@0_[d/s(r)]-{p}'),
                           [(spec.ProteinSpecification('A', 0, spec.DomainDefinition('d', 's', 'r')),
                             spec.SpecificationResolution.residue)],
                           sta.StateModifier.neutral,
                           True),

        ResolutionTestCase(sta.state_from_string('Agene@0_[d/s(r)]-{p}'),
                           [(spec.DnaSpecification('A', 0, spec.DomainDefinition('d', 's', 'r')),
                             spec.SpecificationResolution.residue)],
                           sta.StateModifier.neutral,
                           True),

        ResolutionTestCase(sta.state_from_string('AmRNA@0_[d/s(r)]-{p}'),
                           [(spec.RnaSpecification('A', 0, spec.DomainDefinition('d', 's', 'r')),
                             spec.SpecificationResolution.residue)],
                           sta.StateModifier.neutral,
                           True),

        ResolutionTestCase(sta.state_from_string('A@0_[n]-{p}'),
                           [(spec.ProteinSpecification('A', 0, spec.DomainDefinition('n', None, None)),
                             spec.SpecificationResolution.residue)],
                           sta.StateModifier.neutral,
                           False),

        ResolutionTestCase(sta.state_from_string('A@0_[d]--B@1_[d]'),
                           [(spec.ProteinSpecification('A', 0, spec.DomainDefinition('d', None, None)),
                             spec.SpecificationResolution.domain),
                            (spec.ProteinSpecification('B', 1, spec.DomainDefinition('d', None, None)),
                             spec.SpecificationResolution.domain)
                            ],
                           None,
                           True),

        ResolutionTestCase(sta.state_from_string('A@0_[d/s]--B@1_[d/s]'),
                           [(spec.ProteinSpecification('A', 0, spec.DomainDefinition('d', 's', None)),
                             spec.SpecificationResolution.domain),
                            (spec.ProteinSpecification('B', 1, spec.DomainDefinition('d', 's', None)),
                             spec.SpecificationResolution.domain)
                            ],
                           None,
                           False),

        ResolutionTestCase(sta.state_from_string('A@0--B@1'),
                           [(spec.ProteinSpecification('A', 0, spec.DomainDefinition(None, None, None)),
                             spec.SpecificationResolution.domain),
                            (spec.ProteinSpecification('B', 1, spec.DomainDefinition(None, None, None)),
                             spec.SpecificationResolution.domain)
                            ],
                           None,
                           False)
    ]

def test_resolution_and_default_modifier(test_case_resolution_and_neutral_modifier):
    for the_case in test_case_resolution_and_neutral_modifier:
        is_resolution_and_neutral_modifier_correct(the_case)

def is_resolution_and_neutral_modifier_correct(the_case):
    for component in the_case.state.components():
        assert (component, the_case.state.elemental_resolution(component)) in the_case.expected_resolution
    assert the_case.state.neutral_modifier == the_case.expected_neutral_modifier
    assert the_case.state.is_elemental() == the_case.expected_is_elemental

def test_unboundstate():
    sta.state_from_string('A@0_[m]--0')
@pytest.fixture
def test_case_indirect_subset():
    return [HierarchyTestCase(sta.state_from_string('A@0'),
                      [sta.state_from_string('A@0_[m]--[n]')]),
    ]

# the STATE_DEFINITION gets changed here!!!
def test_indirect_subset(test_case_indirect_subset):
    sta.STATE_DEFINITION = [

        sta.StateDefinition('interaction-state',
                            '$x--$y',
                            { '$x': (spec.Specification, spec.SpecificationResolution.domain),
                              '$y': (spec.Specification, spec.SpecificationResolution.domain) },
                            ['component-state']
                            ),

        sta.StateDefinition('self-interaction-state',
                            '$x--[$y]',
                            { '$x': (spec.Specification, spec.SpecificationResolution.domain),
                              '$y': (spec.DomainDefinition, spec.SpecificationResolution.domain) },
                            ['interaction-state']),

        sta.StateDefinition('component-state',
                            '$x',
                            { '$x': (spec.Specification, spec.SpecificationResolution.component) },
                            [],
                            ),
    ]

    for the_case in test_case_indirect_subset:
        is_hierarchy_correct(the_case)


