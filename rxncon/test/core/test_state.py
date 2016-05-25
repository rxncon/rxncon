import pytest
import rxncon.syntax.rxncon_from_string as rfs
import rxncon.core.state as sta
from collections import namedtuple


HierarchyTestCase = namedtuple('HierarchyTestCase', ['state', 'superspecification_of'])


def test_hierarchy(the_case_hierarchy):
    for test_case in the_case_hierarchy:
        is_hierarchy_correct(test_case)


def is_hierarchy_correct(test_case):
    assert test_case.superspecification_of
    for super_spec in test_case.superspecification_of:
        assert test_case.state.is_superspecification_of(super_spec)
        if test_case.state != super_spec:
            assert not test_case.state.is_subspecification_of(super_spec)
        assert super_spec.is_subspecification_of(test_case.state)


@pytest.fixture
def the_case_hierarchy():
    return [
            HierarchyTestCase(sta.state_from_string('A-{p}'),
                              [sta.state_from_string('A_[n]-{p}')]),

            HierarchyTestCase(sta.state_from_string('A--B'),
                              [sta.state_from_string('A_[n]--B'), sta.state_from_string('A--B_[m]'), sta.state_from_string('A_[n]--B_[m]')],
                              ),

            HierarchyTestCase(sta.state_from_string('A_[n]--B'),
                              [sta.state_from_string('A_[n]--B_[m]')]),

            HierarchyTestCase(sta.state_from_string('A--B_[m]'),
                              [sta.state_from_string('A_[n]--B_[m]')]),

            HierarchyTestCase(sta.state_from_string('A'),
                              [sta.state_from_string('A'),
                               sta.state_from_string('A--B'), sta.state_from_string('A_[d]--B'), sta.state_from_string('A_[d/s]--B'),
                               sta.state_from_string('A_[d/s(r)]--B'),
                               sta.state_from_string('A-{p}'), sta.state_from_string('A_[d]-{p}'), sta.state_from_string('A_[d/s]-{p}'),
                               sta.state_from_string('A_[d/s(r)]-{p}'), sta.state_from_string('A_[(r)]-{p}')
                               ]),

            HierarchyTestCase(sta.state_from_string('B'),
                              [sta.state_from_string('B'),
                               sta.state_from_string('A--B'), sta.state_from_string('A--B_[d]'), sta.state_from_string('A--B_[d/s]'),
                               sta.state_from_string('A--B_[d/s(r)]')
                               ]),

            HierarchyTestCase(sta.state_from_string('[INPUT]'),
                              [sta.state_from_string('[INPUT]')])
    ]


def test_not_super_or_subspecification(the_case_no_hierarchy):
    for test_case in the_case_no_hierarchy:
        for spec in test_case.superspecification_of:
            assert not test_case.state.is_superspecification_of(spec)
            assert not test_case.state.is_subspecification_of(spec)
            assert not spec.is_superspecification_of(test_case.state)
            assert not spec.is_subspecification_of(test_case.state)


@pytest.fixture
def the_case_no_hierarchy():
    return [
        HierarchyTestCase(sta.state_from_string('A_[n]--B'),
                          [sta.state_from_string('A--B_[m]')]),

        HierarchyTestCase(sta.state_from_string('A_[n]--B'),
                          [sta.state_from_string('A_[n]-{P}')]),

        HierarchyTestCase(sta.state_from_string('A--B'),
                          [sta.state_from_string('A-{P}')]),

        HierarchyTestCase(sta.state_from_string('A_[d1]-{P}'),
                          [sta.state_from_string('A_[d2]-{P}')]),

        HierarchyTestCase(sta.state_from_string('A_[d/s1]-{P}'),
                          [sta.state_from_string('A_[d/s2]-{P}')]),

        HierarchyTestCase(sta.state_from_string('A_[d/s(r1)]-{P}'),
                          [sta.state_from_string('A_[d/s(r2)]-{P}')]),

        HierarchyTestCase(sta.state_from_string('A_[d1]--B'),
                          [sta.state_from_string('A_[d2]--B')]),

        HierarchyTestCase(sta.state_from_string('A_[d/s1]--B'),
                          [sta.state_from_string('A_[d/s2]--B')]),

        HierarchyTestCase(sta.state_from_string('A_[d/s(r1)]--B'),
                          [sta.state_from_string('A_[d/s(r2)]--B')]),

        HierarchyTestCase(sta.state_from_string('A'),
                          [sta.state_from_string('B')]),

        HierarchyTestCase(sta.state_from_string('[INPUT]'),
                          [sta.state_from_string('B'), sta.state_from_string('A--B'),
                           sta.state_from_string('A-{P}')])
    ]