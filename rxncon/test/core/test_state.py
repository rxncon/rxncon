import pytest
import rxncon.syntax.rxncon_from_string as rfs
from collections import namedtuple


HierarchyTestCase = namedtuple('HierarchyTestCase', ['state', 'superspecification_of'])


#def test_subspecification(the_case_subspecification):
#    for test_case in the_case_subspecification:
#        for sub_spec in test_case.subspecification_of:
#            assert test_case.state.is_subspecification_of(sub_spec)
#            assert not test_case.state.is_superspecification_of(sub_spec)


def test_hierarchy(the_case_hierarchy):
    for test_case in the_case_hierarchy:
        is_hierarchy_correct(test_case)


def is_hierarchy_correct(test_case):
    assert test_case.superspecification_of
    for super_spec in test_case.superspecification_of:
        assert test_case.state.is_superspecification_of(super_spec)
        if test_case.state != super_spec:
            assert not test_case.state.is_subspecification_of(super_spec)
        # if not super_spec.is_subspecification_of(test_case.state):
        #     print()
        assert super_spec.is_subspecification_of(test_case.state)


@pytest.fixture
def the_case_hierarchy():
    return [
            HierarchyTestCase(rfs.state_from_string('A-{p}'),
                              [rfs.state_from_string('A_[n]-{p}')]),

            HierarchyTestCase(rfs.state_from_string('A--B'),
                              [rfs.state_from_string('A_[n]--B'), rfs.state_from_string('A--B_[m]'), rfs.state_from_string('A_[n]--B_[m]')],
                              ),

            HierarchyTestCase(rfs.state_from_string('A_[n]--B'),
                              [rfs.state_from_string('A_[n]--B_[m]')]),

            HierarchyTestCase(rfs.state_from_string('A--B_[m]'),
                              [rfs.state_from_string('A_[n]--B_[m]')]),

            HierarchyTestCase(rfs.state_from_string('A'),
                              [rfs.state_from_string('A'),
                               rfs.state_from_string('A--B'), rfs.state_from_string('A_[d]--B'), rfs.state_from_string('A_[d/s]--B'),
                               rfs.state_from_string('A_[d/s(r)]--B'),
                               rfs.state_from_string('A-{p}'), rfs.state_from_string('A_[d]-{p}'), rfs.state_from_string('A_[d/s]-{p}'),
                               rfs.state_from_string('A_[d/s(r)]-{p}'), rfs.state_from_string('A_[(r)]-{p}')
                               ]),
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
        HierarchyTestCase(rfs.state_from_string('A_[n]--B'),
                          [rfs.state_from_string('A--B_[m]')]),
    ]

#def test_superspec_subspec_ppi():
    # Superspecification
    #assert rfs.state_from_string('A--B').is_superspecification_of(rfs.state_from_string('A_[n]--B'))
    #assert rfs.state_from_string('A--B').is_superspecification_of(rfs.state_from_string('A--B_[m]'))
    #assert rfs.state_from_string('A--B').is_superspecification_of(rfs.state_from_string('A_[n]--B_[m]'))

    #assert not rfs.state_from_string('A_[n]--B').is_superspecification_of(rfs.state_from_string('A--B_[m]'))
    #assert not rfs.state_from_string('A--B_[m]').is_superspecification_of(rfs.state_from_string('A_[n]--B'))
    #assert rfs.state_from_string('A_[n]--B').is_superspecification_of(rfs.state_from_string('A_[n]--B_[m]'))
    #assert rfs.state_from_string('A--B_[m]').is_superspecification_of(rfs.state_from_string('A_[n]--B_[m]'))

    # Subspecification
    #assert rfs.state_from_string('A_[n]--B').is_subspecification_of(rfs.state_from_string('A--B'))
    #assert rfs.state_from_string('A--B_[m]').is_subspecification_of(rfs.state_from_string('A--B'))
    #assert rfs.state_from_string('A_[n]--B_[m]').is_subspecification_of(rfs.state_from_string('A--B'))

    #assert not rfs.state_from_string('A--B_[m]').is_subspecification_of(rfs.state_from_string('A_[n]--B'))
    #assert not rfs.state_from_string('A_[n]--B').is_subspecification_of(rfs.state_from_string('A--B_[m]'))
    #assert rfs.state_from_string('A_[n]--B_[m]').is_subspecification_of(rfs.state_from_string('A_[n]--B'))
    #assert rfs.state_from_string('A_[n]--B_[m]').is_subspecification_of(rfs.state_from_string('A--B_[m]'))

# def test_superspec_subspec_component_state():
#     assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A'))
#
#     assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A--B'))
#     assert rfs.state_from_string('B').is_superspecification_of(rfs.state_from_string('A--B'))
#     assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A_[d]--B'))
#     assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A_[d/s]--B'))
#     assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A_[d/s(r)]--B'))
#
#     assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A-{p}'))
#     assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A_[d]-{p}'))
#     assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A_[d/s]-{p}'))
#     assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A_[d/s(r)]-{p}'))
#     assert rfs.state_from_string('A').is_superspecification_of(rfs.state_from_string('A_[(r)]-{p}'))
#
#     assert rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A'))
#     assert not rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A--B'))
#     assert not rfs.state_from_string('B').is_subspecification_of(rfs.state_from_string('A--B'))
#     assert not rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A_[d]--B'))
#     assert not rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A_[d/s]--B'))
#     assert not rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A_[d/s(r)]--B'))
#
#     assert not rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A-{p}'))
#     assert not rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A_[d]-{p}'))
#     assert not rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A_[d/s]-{p}'))
#     assert not rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A_[d/s(r)]-{p}'))
#     assert not rfs.state_from_string('A').is_subspecification_of(rfs.state_from_string('A_[(r)]-{p}'))
