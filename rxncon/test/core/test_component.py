import pytest
import typecheck as tc

from collections import namedtuple

import rxncon.core.specification as com
import rxncon.syntax.rxncon_from_string as cfs

ComponentSubspecTestCase = namedtuple("ComponentSubspecTestCase", ['state', 'is_subspecification', 'is_no_subspecification'])
ComponentInitTestCase = namedtuple('ComponentInitTestCase', ['seeded', 'empty'])


def test_component_without_name_raises():
    with pytest.raises(tc.InputParameterError):
        component = com.ProteinSpecification(None, None, None, None)

    with pytest.raises(tc.InputParameterError):
        component = com.RnaSpecification(None, None, None, None)


def test_component_subspecification(the_case_subspecification):
    for the_case in the_case_subspecification:
        is_hierarchy_correct(the_case)


def test_component_initialized(the_case_initialized):
    for the_case in the_case_initialized:
        is_correct_initialized(the_case)



### SPECIFICATION, SUBSPECIFICATION TESTS ###
# All of these components describe the same protein at different specification resolutions.
no_domain = cfs.specification_from_string('A')
domain = cfs.specification_from_string('A_[b]')
domain_subdomain = cfs.specification_from_string('A_[b/c]')
domain_residue = cfs.specification_from_string('A_[b(d)]')
domain_subdomain_residue = cfs.specification_from_string('A_[b/c(d)]')
residue = cfs.specification_from_string('A_[(d)]')


@pytest.fixture
def the_case_subspecification():
    return [
        ComponentSubspecTestCase(no_domain, [no_domain, domain, domain_subdomain, domain_subdomain_residue, domain_residue, residue],  []),
        ComponentSubspecTestCase(domain, [domain, domain_subdomain, domain_subdomain_residue, domain_residue],  [no_domain, residue]),
        ComponentSubspecTestCase(domain_subdomain, [domain_subdomain, domain_subdomain_residue],  [no_domain, domain, domain_residue, residue]),
        ComponentSubspecTestCase(domain_residue, [domain_residue, domain_subdomain_residue, residue], [no_domain, domain, domain_subdomain]),
        ComponentSubspecTestCase(domain_subdomain_residue, [domain_subdomain_residue, domain_residue, residue],  [no_domain, domain, domain_subdomain]),
        ComponentSubspecTestCase(residue, [domain_residue, domain_subdomain_residue, residue], [no_domain, domain, domain_subdomain])
    ]


@pytest.fixture
def the_case_initialized():
    return [
        ComponentInitTestCase([no_domain.name], [no_domain.domain, no_domain.subdomain, no_domain.residue]),
        ComponentInitTestCase([domain.name, domain.domain], [domain.subdomain, domain.residue]),
        ComponentInitTestCase([domain_subdomain.name, domain_subdomain.domain, domain_subdomain.subdomain], [domain_subdomain.residue]),
        ComponentInitTestCase([domain_subdomain_residue.name, domain_subdomain_residue.domain,
                               domain_subdomain_residue.subdomain, domain_subdomain_residue.residue], []),
        ComponentInitTestCase([domain_residue.name, domain_residue.domain, domain_residue.residue], [domain_residue.subdomain]),
        ComponentInitTestCase([residue.name, residue.residue], [residue.domain, residue.subdomain])

    ]


def is_correct_initialized(the_case):
    for initializated in the_case.seeded:
        assert initializated
    for empty in the_case.empty:
        assert not empty


def is_hierarchy_correct(the_case):
    for subspec in the_case.is_subspecification:
        assert subspec.is_subspecification_of(the_case.state)
        if subspec != the_case.state:
            if (subspec.residue is None and the_case.state.residue is None) or (subspec.residue != the_case.state.residue):
                assert not subspec.is_superspecification_of(the_case.state)
            else:
                assert subspec.is_superspecification_of(the_case.state)
        assert the_case.state.is_superspecification_of(subspec)

    for no_subspec in the_case.is_no_subspecification:
        assert not no_subspec.is_subspecification_of(the_case.state)


# Tests for a.is_subspecification_of(b), i.e. a is more precisely specified than b.
#def test_subspecification_no_domain():
    # assert no_domain.is_subspecification_of(no_domain)
    # assert domain.is_subspecification_of(no_domain)
    # assert domain_subdomain.is_subspecification_of(no_domain)
    # assert domain_residue.is_subspecification_of(no_domain)
    # assert domain_subdomain_residue.is_subspecification_of(no_domain)
    # assert residue.is_subspecification_of(no_domain)


# def test_subspecification_domain():
#     assert not no_domain.is_subspecification_of(domain)
#     assert domain.is_subspecification_of(domain)
#     assert domain_subdomain.is_subspecification_of(domain)
#     assert domain_residue.is_subspecification_of(domain)
#     assert domain_subdomain_residue.is_subspecification_of(domain)
#     assert not residue.is_subspecification_of(domain)


# def test_subspecification_domain_subdomain():
#     assert not no_domain.is_subspecification_of(domain_subdomain)
#     assert not domain.is_subspecification_of(domain_subdomain)
#     assert domain_subdomain.is_subspecification_of(domain_subdomain)
#     assert not domain_residue.is_subspecification_of(domain_subdomain)
#     assert domain_subdomain_residue.is_subspecification_of(domain_subdomain)
#     assert not residue.is_subspecification_of(domain_subdomain)

#
# def test_subspecification_domain_residue():
#     assert not no_domain.is_subspecification_of(domain_residue)
#     assert not domain.is_subspecification_of(domain_residue)
#     assert not domain_subdomain.is_subspecification_of(domain_residue)
#     assert domain_residue.is_subspecification_of(domain_residue)
#     assert domain_subdomain_residue.is_subspecification_of(domain_residue)
#     assert residue.is_subspecification_of(domain_residue)


# def test_subspecification_domain_subdomain_residue():
#     assert not no_domain.is_subspecification_of(domain_subdomain_residue)
#     assert not domain.is_subspecification_of(domain_subdomain_residue)
#     assert not domain_subdomain.is_subspecification_of(domain_subdomain_residue)
#     assert domain_residue.is_subspecification_of(domain_subdomain_residue)
#     assert domain_subdomain_residue.is_subspecification_of(domain_subdomain_residue)
#     assert residue.is_subspecification_of(domain_subdomain_residue)


# def test_subspecification_residue():
#     assert not no_domain.is_subspecification_of(residue)
#     assert not domain.is_subspecification_of(residue)
#     assert not domain_subdomain.is_subspecification_of(residue)
#     assert domain_residue.is_subspecification_of(residue)
#     assert domain_subdomain_residue.is_subspecification_of(residue)
#     assert residue.is_subspecification_of(residue)


# Tests for a.is_superspecification_of(b), i.e. a is less precisely specified than b.
# def test_superspecification_no_domain():
#     assert no_domain.is_superspecification_of(no_domain)
#     assert not domain.is_superspecification_of(no_domain)
#     assert not domain_subdomain.is_superspecification_of(no_domain)
#     assert not domain_residue.is_superspecification_of(no_domain)
#     assert not domain_subdomain_residue.is_superspecification_of(no_domain)
#     assert not residue.is_superspecification_of(no_domain)


def test_superspecification_domain():
    assert no_domain.is_superspecification_of(domain)
    assert domain.is_superspecification_of(domain)
    assert not domain_subdomain.is_superspecification_of(domain)
    assert not domain_residue.is_superspecification_of(domain)
    assert not domain_subdomain_residue.is_superspecification_of(domain)
    assert not residue.is_superspecification_of(domain)


def test_superspecification_domain_subdomain():
    assert no_domain.is_superspecification_of(domain_subdomain)
    assert domain.is_superspecification_of(domain_subdomain)
    assert domain_subdomain.is_superspecification_of(domain_subdomain)
    assert not domain_residue.is_superspecification_of(domain_subdomain)
    assert not domain_subdomain_residue.is_superspecification_of(domain_subdomain)
    assert not residue.is_superspecification_of(domain_subdomain)


def test_superspecification_domain_residue():
    assert no_domain.is_superspecification_of(domain_residue)
    assert domain.is_superspecification_of(domain_residue)
    assert not domain_subdomain.is_superspecification_of(domain_residue)
    assert domain_residue.is_superspecification_of(domain_residue)
    assert domain_subdomain_residue.is_superspecification_of(domain_residue)
    assert residue.is_superspecification_of(domain_residue)


def test_superspecification_domain_subdomain_residue():
    assert no_domain.is_superspecification_of(domain_subdomain_residue)
    assert domain.is_superspecification_of(domain_subdomain_residue)
    assert domain_subdomain.is_superspecification_of(domain_subdomain_residue)
    assert domain_residue.is_superspecification_of(domain_subdomain_residue)
    assert domain_subdomain_residue.is_superspecification_of(domain_subdomain_residue)
    assert residue.is_superspecification_of(domain_subdomain_residue)


def test_superspecification_residue():
    assert no_domain.is_superspecification_of(residue)
    assert not domain.is_superspecification_of(residue)
    assert not domain_subdomain.is_superspecification_of(residue)
    assert domain_residue.is_superspecification_of(residue)
    assert domain_subdomain_residue.is_superspecification_of(residue)
    assert residue.is_superspecification_of(residue)


# Tests for different proteins with identical domains etc.
def test_different_proteins():
    first_component = cfs.specification_from_string('A_[b/d]')
    second_component = cfs.specification_from_string('B_[b]')

    assert not first_component.is_equivalent_to(second_component)
    assert not first_component.is_subspecification_of(second_component)
    assert not first_component.is_superspecification_of(second_component)

    assert not second_component.is_equivalent_to(first_component)
    assert not second_component.is_subspecification_of(first_component)
    assert not second_component.is_superspecification_of(first_component)

