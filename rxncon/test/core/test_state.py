import pytest
import rxncon.syntax.rxncon_from_string as rfs
import rxncon.core.state as sta
import rxncon.core.specification as spec
from collections import namedtuple


HierarchyTestCase = namedtuple('HierarchyTestCase', ['state', 'superset_of'])


def test_hierarchy(the_case_hierarchy):
    for test_case in the_case_hierarchy:
        is_hierarchy_correct(test_case)


def is_hierarchy_correct(test_case):
    assert test_case.superset_of
    for super_state in test_case.superset_of:
        assert test_case.state.is_superset_of(super_state)
        if test_case.state != super_state:
            assert not test_case.state.is_subset_of(super_state)
        assert super_state.is_subspecification_of(test_case.state)


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
        for spec in test_case.superspecification_of:
            assert not test_case.state.is_superset_of(spec)
            assert not test_case.state.is_subset_of(spec)
            assert not spec.is_superspecification_of(test_case.state)
            assert not spec.is_subspecification_of(test_case.state)


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

        HierarchyTestCase(sta.state_from_string('[INPUT]'),
                          [sta.state_from_string('B@1'), sta.state_from_string('A@0--B@1'),
                           sta.state_from_string('A@0-{P}')])
    ]

StateTestCase = namedtuple('StateTestCase', ['state', 'expected_specifications', 'expected_string'])

@pytest.fixture
def the_case_state():
    return [
        StateTestCase(sta.state_from_string('Agene@0'),
                      [spec.DnaSpecification('A', 0, spec.DomainResolution(None, None, None))],
                      'DnaSpecification: Agene@0'),
        StateTestCase(sta.state_from_string('Agene@0_[d/s(r)]'),
                      [spec.DnaSpecification('A', 0, spec.DomainResolution('d', 's', 'r'))],
                      'DnaSpecification: Agene@0_[d/s(r)]'),
        StateTestCase(sta.state_from_string('AmRNA@0'),
                      [spec.RnaSpecification('A', 0, spec.DomainResolution(None, None, None))],
                      'RnaSpecification: AmRNA@0'),
        StateTestCase(sta.state_from_string('AmRNA@0_[d/s(r)]'),
                      [spec.RnaSpecification('A', 0, spec.DomainResolution('d', 's', 'r'))],
                      'RnaSpecification: AmRNA@0_[d/s(r)]'),
        StateTestCase(sta.state_from_string('A@0'),
                      [spec.ProteinSpecification('A', 0, spec.DomainResolution(None, None, None))],
                      'ProteinSpecification: A@0'),
        StateTestCase(sta.state_from_string('A@0_[d/s(r)]'),
                      [spec.ProteinSpecification('A', 0, spec.DomainResolution('d', 's', 'r'))],
                      'ProteinSpecification: A@0_[d/s(r)]')
    ]

def test_state_building(the_case_state):
    for the_case in the_case_state:
        is_state_correct(the_case)

def is_state_correct(the_case):
    assert all(the_case.state.variables[variable] in the_case.expected_specifications for variable in the_case.state.variables)
    assert str(the_case.state) == the_case.expected_string


