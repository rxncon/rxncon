import pytest
import typecheck as tc

import rxncon.core.specification as com
import rxncon.syntax.rxncon_from_string as cfs


def test_component_without_name_raises():
    with pytest.raises(tc.InputParameterError):
        component = com.Specification(None, None, None, None)


### SPECIFICATION, SUBSPECIFICATION TESTS ###
# All of these components describe the same protein at different specification resolutions.
no_domain = cfs.specification_from_string('A')
domain = cfs.specification_from_string('A_[b]')
domain_subdomain = cfs.specification_from_string('A_[b/c]')
domain_residue = cfs.specification_from_string('A_[b(d)]')
domain_subdomain_residue = cfs.specification_from_string('A_[b/c(d)]')
residue = cfs.specification_from_string('A_[(d)]')


def test_component_initialization():
    assert no_domain.name
    assert not no_domain.domain
    assert not no_domain.subdomain
    assert not no_domain.residue

    assert domain.name
    assert domain.domain
    assert not domain.subdomain
    assert not domain.residue

    assert domain_subdomain.name
    assert domain_subdomain.domain
    assert domain_subdomain.subdomain
    assert not domain_subdomain.residue

    assert domain_subdomain_residue.name
    assert domain_subdomain_residue.domain
    assert domain_subdomain_residue.subdomain
    assert domain_subdomain_residue.residue

    assert domain_residue.name
    assert domain_residue.domain
    assert not domain_residue.subdomain
    assert domain_residue.residue

    assert residue.name
    assert not residue.domain
    assert not residue.subdomain
    assert residue.residue


# Tests for a.is_subspecification_of(b), i.e. a is more precisely specified than b.
def test_subspecification_no_domain():
    assert no_domain.is_subspecification_of(no_domain)
    assert domain.is_subspecification_of(no_domain)
    assert domain_subdomain.is_subspecification_of(no_domain)
    assert domain_residue.is_subspecification_of(no_domain)
    assert domain_subdomain_residue.is_subspecification_of(no_domain)
    assert residue.is_subspecification_of(no_domain)


def test_subspecification_domain():
    assert not no_domain.is_subspecification_of(domain)
    assert domain.is_subspecification_of(domain)
    assert domain_subdomain.is_subspecification_of(domain)
    assert domain_residue.is_subspecification_of(domain)
    assert domain_subdomain_residue.is_subspecification_of(domain)
    assert not residue.is_subspecification_of(domain)


def test_subspecification_domain_subdomain():
    assert not no_domain.is_subspecification_of(domain_subdomain)
    assert not domain.is_subspecification_of(domain_subdomain)
    assert domain_subdomain.is_subspecification_of(domain_subdomain)
    assert not domain_residue.is_subspecification_of(domain_subdomain)
    assert domain_subdomain_residue.is_subspecification_of(domain_subdomain)
    assert not residue.is_subspecification_of(domain_subdomain)


def test_subspecification_domain_residue():
    assert not no_domain.is_subspecification_of(domain_residue)
    assert not domain.is_subspecification_of(domain_residue)
    assert not domain_subdomain.is_subspecification_of(domain_residue)
    assert domain_residue.is_subspecification_of(domain_residue)
    assert domain_subdomain_residue.is_subspecification_of(domain_residue)
    assert residue.is_subspecification_of(domain_residue)


def test_subspecification_domain_subdomain_residue():
    assert not no_domain.is_subspecification_of(domain_subdomain_residue)
    assert not domain.is_subspecification_of(domain_subdomain_residue)
    assert not domain_subdomain.is_subspecification_of(domain_subdomain_residue)
    assert domain_residue.is_subspecification_of(domain_subdomain_residue)
    assert domain_subdomain_residue.is_subspecification_of(domain_subdomain_residue)
    assert residue.is_subspecification_of(domain_subdomain_residue)


def test_subspecification_residue():
    assert not no_domain.is_subspecification_of(residue)
    assert not domain.is_subspecification_of(residue)
    assert not domain_subdomain.is_subspecification_of(residue)
    assert domain_residue.is_subspecification_of(residue)
    assert domain_subdomain_residue.is_subspecification_of(residue)
    assert residue.is_subspecification_of(residue)


# Tests for a.is_superspecification_of(b), i.e. a is less precisely specified than b.
def test_superspecification_no_domain():
    assert no_domain.is_superspecification_of(no_domain)
    assert not domain.is_superspecification_of(no_domain)
    assert not domain_subdomain.is_superspecification_of(no_domain)
    assert not domain_residue.is_superspecification_of(no_domain)
    assert not domain_subdomain_residue.is_superspecification_of(no_domain)
    assert not residue.is_superspecification_of(no_domain)


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

