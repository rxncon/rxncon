import pytest
import rxncon.core.state as sta
from rxncon.core.spec import RnaSpec, SpecificationResolution, ProteinSpec, DNASpec, Domain, Spec
from collections import namedtuple


HierarchyTestCase = namedtuple('HierarchyTestCase', ['state', 'superset_of'])


def test_hierarchy(the_case_hierarchy, the_case_hierarchy_structured):
    for test_case in the_case_hierarchy +the_case_hierarchy_structured:
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
        HierarchyTestCase(state_from_string('B'),
                          [state_from_string('B'),
                           state_from_string('A--B'), state_from_string('A--B_[d]'),
                           state_from_string('A--B_[d/s]'),
                           state_from_string('A--B_[d/s(r)]')
                           ]),
        HierarchyTestCase(state_from_string('A-{p}'),
                          [state_from_string('A_[n]-{p}')]),

        HierarchyTestCase(state_from_string('A--B'),
                          [state_from_string('A_[n]--B'), state_from_string('A--B_[m]'),
                           state_from_string('A_[n]--B_[m]')],
                          ),

        HierarchyTestCase(state_from_string('A_[n]--B'),
                          [state_from_string('A_[n]--B_[m]')]),

        HierarchyTestCase(state_from_string('A--B_[m]'),
                          [state_from_string('A_[n]--B_[m]')]),

        HierarchyTestCase(state_from_string('A'),
                          [state_from_string('A'),
                           state_from_string('A--B'), state_from_string('A_[d]--B'),
                           state_from_string('A_[d/s]--B'),
                           state_from_string('A_[d/s(r)]--B'),
                           state_from_string('A-{p}'), state_from_string('A_[d]-{p}'),
                           state_from_string('A_[d/s]-{p}'),
                           state_from_string('A_[d/s(r)]-{p}'), state_from_string('A_[(r)]-{p}')
                           ]),

        HierarchyTestCase(state_from_string('B'),
                          [state_from_string('B'),
                           state_from_string('A--B'), state_from_string('A--B_[d]'),
                           state_from_string('A--B_[d/s]'),
                           state_from_string('A--B_[d/s(r)]')
                           ]),
    ]


@pytest.fixture
def the_case_hierarchy_structured():
    return [
            HierarchyTestCase(state_from_string('A@0-{p}'),
                              [state_from_string('A@0_[n]-{p}')]),

            HierarchyTestCase(state_from_string('A@0--B@1'),
                              [state_from_string('A@0_[n]--B@1'), state_from_string('A@0--B@1_[m]'), state_from_string('A@0_[n]--B@1_[m]')],
                              ),

            HierarchyTestCase(state_from_string('A@0_[n]--B@1'),
                              [state_from_string('A@0_[n]--B@1_[m]')]),

            HierarchyTestCase(state_from_string('A@0--B@1_[m]'),
                              [state_from_string('A@0_[n]--B@1_[m]')]),

            HierarchyTestCase(state_from_string('A@0'),
                              [state_from_string('A@0'),
                               state_from_string('A@0--B@1'), state_from_string('A@0_[d]--B@1'), state_from_string('A@0_[d/s]--B@1'),
                               state_from_string('A@0_[d/s(r)]--B@1'),
                               state_from_string('A@0-{p}'), state_from_string('A@0_[d]-{p}'), state_from_string('A@0_[d/s]-{p}'),
                               state_from_string('A@0_[d/s(r)]-{p}'), state_from_string('A@0_[(r)]-{p}')
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
        HierarchyTestCase(state_from_string('A@0_[n]--B@1'),
                          [state_from_string('A@0--B@1_[m]')]),

        HierarchyTestCase(state_from_string('A@0_[n]--B@1'),
                          [state_from_string('A@0_[n]-{P}')]),

        HierarchyTestCase(state_from_string('A@0--B@1'),
                          [state_from_string('A@0-{P}')]),

        HierarchyTestCase(state_from_string('A@0_[d1]-{P}'),
                          [state_from_string('A@0_[d2]-{P}')]),

        HierarchyTestCase(state_from_string('A@0_[d/s1]-{P}'),
                          [state_from_string('A@0_[d/s2]-{P}')]),

        HierarchyTestCase(state_from_string('A@0_[d/s(r1)]-{P}'),
                          [state_from_string('A@0_[d/s(r2)]-{P}')]),

        HierarchyTestCase(state_from_string('A@0_[d1]--B@1'),
                          [state_from_string('A@0_[d2]--B@1')]),

        HierarchyTestCase(state_from_string('A@0_[d/s1]--B@1'),
                          [state_from_string('A@0_[d/s2]--B@1')]),

        HierarchyTestCase(state_from_string('A@0_[d/s(r1)]--B@1'),
                          [state_from_string('A@0_[d/s(r2)]--B@1')]),

        HierarchyTestCase(state_from_string('A@0'),
                          [state_from_string('B@1')]),

        #HierarchyTestCase(sta.state_from_string('[INPUT]'),
        #                  [sta.state_from_string('B@1'), sta.state_from_string('A@0--B@1'),
        #                   sta.state_from_string('A@0-{P}')])
    ]


StateTestCase = namedtuple('StateTestCase', ['state', 'expected_specifications', 'expected_string'])

@pytest.fixture
def the_case_state():
    return [

        StateTestCase(state_from_string('Agene'),
                      [DNASpec('A', None, Domain(None, None, None))],
                      'Agene'),
        StateTestCase(state_from_string('Agene_[d/s(r)]'),
                      [DNASpec('A', None, Domain('d', 's', 'r'))],
                      'Agene_[d/s(r)]'),
        StateTestCase(state_from_string('AmRNA'),
                      [RnaSpec('A', None, Domain(None, None, None))],
                      'AmRNA'),
        StateTestCase(state_from_string('AmRNA_[d/s(r)]'),
                      [RnaSpec('A', None, Domain('d', 's', 'r'))],
                      'AmRNA_[d/s(r)]'),
        StateTestCase(state_from_string('A'),
                      [ProteinSpec('A', None, Domain(None, None, None))],
                      'A'),
        StateTestCase(state_from_string('A_[d/s(r)]'),
                      [ProteinSpec('A', None, Domain('d', 's', 'r'))],
                      'A_[d/s(r)]')
    ]


@pytest.fixture
def the_case_state_structured():
    return [
        StateTestCase(state_from_string('Agene@0'),
                      [DNASpec('A', 0, Domain(None, None, None))],
                      'Agene@0'),
        StateTestCase(state_from_string('Agene@0_[d/s(r)]'),
                      [DNASpec('A', 0, Domain('d', 's', 'r'))],
                      'Agene@0_[d/s(r)]'),
        StateTestCase(state_from_string('AmRNA@0'),
                      [RnaSpec('A', 0, Domain(None, None, None))],
                      'AmRNA@0'),
        StateTestCase(state_from_string('AmRNA@0_[d/s(r)]'),
                      [RnaSpec('A', 0, Domain('d', 's', 'r'))],
                      'AmRNA@0_[d/s(r)]'),
        StateTestCase(state_from_string('A@0'),
                      [ProteinSpec('A', 0, Domain(None, None, None))],
                      'A@0'),
        StateTestCase(state_from_string('A@0_[d/s(r)]'),
                      [ProteinSpec('A', 0, Domain('d', 's', 'r'))],
                      'A@0_[d/s(r)]'),

    ]


def test_state_building(the_case_state, the_case_state_structured):
    for the_case in the_case_state + the_case_state_structured:
        is_state_correct(the_case)


def is_state_correct(the_case):
    assert all(the_case.state.variables[variable] in the_case.expected_specifications for variable in the_case.state.variables)
    assert str(the_case.state) == the_case.expected_string

ResolutionTestCase = namedtuple('ResolutionTestCase', ['state', 'expected_resolution', 'expected_neutral_modifier', 'expected_is_elemental'])

@pytest.fixture
def test_case_resolution_and_neutral_modifier():
    return [
        ResolutionTestCase(state_from_string('A'),
                           [(ProteinSpec('A', None, Domain(None, None, None)),
                             SpecificationResolution.component)],
                           None,
                           True),
        ResolutionTestCase(state_from_string('A_[m]--[n]'),
                           [(ProteinSpec('A', None, Domain('m', None, None)),
                             SpecificationResolution.domain)],
                           None,
                           True),

        ResolutionTestCase(state_from_string('A_[(r)]-{p}'),
                           [(ProteinSpec('A', None, Domain(None, None, 'r')),
                             SpecificationResolution.residue)],
                           StateModifier.neutral,
                           True),

        ResolutionTestCase(state_from_string('A_[n]-{p}'),
                           [(ProteinSpec('A', None, Domain('n', None, None)),
                             SpecificationResolution.residue)],
                           StateModifier.neutral,
                           False),
    ]
@pytest.fixture
def test_case_resolution_and_neutral_modifier_structured():
    return [


        ResolutionTestCase(state_from_string('A@0'),
                           [(ProteinSpec('A', 0, Domain(None, None, None)),
                             SpecificationResolution.component)],
                           None,
                           True),

        ResolutionTestCase(state_from_string('Agene@0'),
                           [(DNASpec('A', 0, Domain(None, None, None)),
                             SpecificationResolution.component)],
                           None,
                           True),

        ResolutionTestCase(state_from_string('AmRNA@0'),
                           [(RnaSpec('A', 0, Domain(None, None, None)),
                             SpecificationResolution.component)],
                           None,
                           True),

        #todo: discuss this What is with the DomainResolution ?
        ResolutionTestCase(state_from_string('A@0_[m]--[n]'),
                           [(ProteinSpec('A', 0, Domain('m', None, None)),
                             SpecificationResolution.domain)],
                           None,
                           True),

        ResolutionTestCase(state_from_string('A@0_[(r)]-{p}'),
                           [(ProteinSpec('A', 0, Domain(None, None, 'r')),
                             SpecificationResolution.residue)],
                           StateModifier.neutral,
                           True),

        ResolutionTestCase(state_from_string('A@0_[d(r)]-{p}'),
                           [(ProteinSpec('A', 0, Domain('d', None, 'r')),
                             SpecificationResolution.residue)],
                           StateModifier.neutral,
                           True),

        ResolutionTestCase(state_from_string('A@0_[d/s(r)]-{p}'),
                           [(ProteinSpec('A', 0, Domain('d', 's', 'r')),
                             SpecificationResolution.residue)],
                           StateModifier.neutral,
                           True),

        ResolutionTestCase(state_from_string('Agene@0_[d/s(r)]-{p}'),
                           [(DNASpec('A', 0, Domain('d', 's', 'r')),
                             SpecificationResolution.residue)],
                           StateModifier.neutral,
                           True),

        ResolutionTestCase(state_from_string('AmRNA@0_[d/s(r)]-{p}'),
                           [(RnaSpec('A', 0, Domain('d', 's', 'r')),
                             SpecificationResolution.residue)],
                           StateModifier.neutral,
                           True),

        ResolutionTestCase(state_from_string('A@0_[n]-{p}'),
                           [(ProteinSpec('A', 0, Domain('n', None, None)),
                             SpecificationResolution.residue)],
                           StateModifier.neutral,
                           False),

        ResolutionTestCase(state_from_string('A@0_[d]--B@1_[d]'),
                           [(ProteinSpec('A', 0, Domain('d', None, None)),
                             SpecificationResolution.domain),
                            (ProteinSpec('B', 1, Domain('d', None, None)),
                             SpecificationResolution.domain)
                            ],
                           None,
                           True),

        ResolutionTestCase(state_from_string('A@0_[d/s]--B@1_[d/s]'),
                           [(ProteinSpec('A', 0, Domain('d', 's', None)),
                             SpecificationResolution.domain),
                            (ProteinSpec('B', 1, Domain('d', 's', None)),
                             SpecificationResolution.domain)
                            ],
                           None,
                           False),

        ResolutionTestCase(state_from_string('A@0--B@1'),
                           [(ProteinSpec('A', 0, Domain(None, None, None)),
                             SpecificationResolution.domain),
                            (ProteinSpec('B', 1, Domain(None, None, None)),
                             SpecificationResolution.domain)
                            ],
                           None,
                           False)
    ]

def test_resolution_and_default_modifier(test_case_resolution_and_neutral_modifier, test_case_resolution_and_neutral_modifier_structured):
    for the_case in test_case_resolution_and_neutral_modifier + test_case_resolution_and_neutral_modifier_structured:
        is_resolution_and_neutral_modifier_correct(the_case)

def is_resolution_and_neutral_modifier_correct(the_case):
    for component in the_case.state.components():
        assert (component, the_case.state._elemental_resolutions(component)) in the_case.expected_resolution
    assert the_case.state.neutral_modifier == the_case.expected_neutral_modifier
    assert the_case.state.is_elemental() == the_case.expected_is_elemental

def test_unboundstate():
    state_from_string('A@0_[m]--0')
@pytest.fixture
def test_case_indirect_subset():
    return [HierarchyTestCase(state_from_string('A@0'),
                              [state_from_string('A@0_[m]--[n]')]),
    ]

# the STATE_DEFINITION gets changed here!!!
def test_indirect_subset(test_case_indirect_subset):
    sta.STATE_DEFS = [

        sta.StateDef('interaction-state',
                            '$x--$y',
                     { '$x': (Spec, SpecificationResolution.domain),
                              '$y': (Spec, SpecificationResolution.domain)},
                     ['component-state']
                     ),

        sta.StateDef('self-interaction-state',
                            '$x--[$y]',
                     { '$x': (Spec, SpecificationResolution.domain),
                              '$y': (Domain, SpecificationResolution.domain)},
                     ['interaction-state']),

        sta.StateDef('component-state',
                            '$x',
                     { '$x': (Spec, SpecificationResolution.component)},
                     [],
                     ),
    ]

    for the_case in test_case_indirect_subset:
        is_hierarchy_correct(the_case)


