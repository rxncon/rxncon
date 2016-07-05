import pytest
import typecheck as tc

from collections import namedtuple

import rxncon.core.specification as com
import rxncon.syntax.specification_from_string as cfs

ComponentInitTestCase = namedtuple('ComponentInitTestCase', ['seeded', 'empty'])
ComponentHierarchicalTestCase = namedtuple("ComponentHierarchicalTestCase", ['state', 'is_superspecification'])
SpecificationEquivalenceTestCase = namedtuple('SpecificationEquivalenceTestCase', ['state', 'equivalent', 'notequivalent'])

### SPECIFICATION, SUBSPECIFICATION TESTS ###
# All of these components describe the same protein at different specification resolutions.
no_domain = cfs.specification_from_string('A@0')
domain = cfs.specification_from_string('A@0_[b]')
domain_subdomain = cfs.specification_from_string('A@0_[b/c]')
domain_residue = cfs.specification_from_string('A@0_[b(d)]')
domain_subdomain_residue = cfs.specification_from_string('A@0_[b/c(d)]')
residue = cfs.specification_from_string('A@0_[(d)]')


def test_component_without_name_raises():
    with pytest.raises(tc.InputParameterError):
        com.ProteinSpecification(None, 0, com.DomainDefinition(None, None, None))

    with pytest.raises(tc.InputParameterError):
        com.RnaSpecification(None,  0, com.DomainDefinition(None, None, None))

    with pytest.raises(tc.InputParameterError):
        com.DnaSpecification(None,  0, com.DomainDefinition(None, None, None))


def test_component_initialized(the_case_initialized):
    for the_case in the_case_initialized:
        is_correct_initialized(the_case)


def test_component_subspecification_superspecification(the_case_superspecification):
    for the_case in the_case_superspecification:
        is_hierarchy_correct(the_case)


def test_components_both_sub_and_superspecification(the_case_both_superspecification_and_subspecification):
    for the_case in the_case_both_superspecification_and_subspecification:
        is_both_sub_and_superspecification(the_case)


def test_component_equivalence(the_case_is_equivalent_to):
    for the_case in the_case_is_equivalent_to:
        is_equivalent(the_case)


@pytest.fixture
def the_case_initialized():
    return [
        ComponentInitTestCase([no_domain.name], [no_domain.spec_resolution.domain, no_domain.spec_resolution.subdomain,
                                                 no_domain.spec_resolution.residue]),
        ComponentInitTestCase([domain.name, domain.spec_resolution.domain], [domain.spec_resolution.subdomain,
                                                                             domain.spec_resolution.residue]),
        ComponentInitTestCase([domain_subdomain.name, domain_subdomain.spec_resolution.domain,
                               domain_subdomain.spec_resolution.subdomain], [domain_subdomain.spec_resolution.residue]),
        ComponentInitTestCase([domain_subdomain_residue.name, domain_subdomain_residue.spec_resolution.domain,
                               domain_subdomain_residue.spec_resolution.subdomain, domain_subdomain_residue.spec_resolution.residue], []),
        ComponentInitTestCase([domain_residue.name, domain_residue.spec_resolution.domain, domain_residue.spec_resolution.residue],
                              [domain_residue.spec_resolution.subdomain]),
        ComponentInitTestCase([residue.name, residue.spec_resolution.residue], [residue.spec_resolution.domain, residue.spec_resolution.subdomain])

    ]


def is_correct_initialized(the_case):
    for initializated in the_case.seeded:
        assert initializated
    for empty in the_case.empty:
        assert not empty


def is_hierarchy_correct(the_case):
    for spec in the_case.is_superspecification:
        assert the_case.state.is_superspecification_of(spec)
        assert spec.is_subspecification_of(the_case.state)
        if spec != the_case.state:
            assert not spec.is_superspecification_of(the_case.state)
        else:
            assert spec.is_superspecification_of(the_case.state)


def is_both_sub_and_superspecification(the_case):
    for spec in the_case.is_superspecification:
        assert the_case.state.is_superspecification_of(spec)
        assert spec.is_superspecification_of(the_case.state)
        assert the_case.state.is_subspecification_of(spec)
        assert spec.is_subspecification_of(the_case.state)


def is_equivalent(the_case):
    for spec in the_case.equivalent:
        assert spec.is_equivalent_to(the_case.state)
        assert the_case.state.is_equivalent_to(spec)

    for spec in the_case.notequivalent:
        assert not spec.is_equivalent_to(the_case.state)
        assert not the_case.state.is_equivalent_to(spec)


@pytest.fixture
def the_case_superspecification():
    return [
        ComponentHierarchicalTestCase(no_domain, [no_domain, domain, domain_subdomain, domain_residue, domain_subdomain_residue]),
        ComponentHierarchicalTestCase(domain, [domain, domain_subdomain, domain_residue, domain_subdomain_residue]),
        ComponentHierarchicalTestCase(domain_subdomain, [domain_subdomain, domain_subdomain_residue]),
    ]


@pytest.fixture
def the_case_both_superspecification_and_subspecification():
    return [
        ComponentHierarchicalTestCase(domain_residue, [domain_residue, domain_subdomain_residue, residue])
    ]


@pytest.fixture
def the_case_is_equivalent_to():
    return [
        SpecificationEquivalenceTestCase(no_domain, [no_domain], [domain, domain_subdomain, domain_residue, domain_subdomain_residue]),
        SpecificationEquivalenceTestCase(domain, [domain], [domain_subdomain, domain_residue, domain_subdomain_residue, residue]),
        SpecificationEquivalenceTestCase(domain_subdomain, [domain_subdomain], [domain_residue, domain_subdomain_residue, residue]),
        SpecificationEquivalenceTestCase(domain_residue, [domain_residue, domain_subdomain_residue, residue], []),
        SpecificationEquivalenceTestCase(cfs.specification_from_string('A@0'), [], [cfs.specification_from_string('B@1'),
                                                                                    cfs.specification_from_string('B@1_[d]'),
                                                                                    cfs.specification_from_string('B@1_[d/s]'),
                                                                                    cfs.specification_from_string('B@1_[d/s(r)]'),
                                                                                    cfs.specification_from_string('B@1_[d(r)]'),
                                                                                    cfs.specification_from_string('B@1_[(r)]')]),
        SpecificationEquivalenceTestCase(cfs.specification_from_string('A@0_[d]'), [], [cfs.specification_from_string('B@1'),
                                                                                        cfs.specification_from_string('B@1_[d]'),
                                                                                        cfs.specification_from_string('B@1_[d/s]'),
                                                                                        cfs.specification_from_string('B@1_[d/s(r)]'),
                                                                                        cfs.specification_from_string('B@1_[d(r)]'),
                                                                                        cfs.specification_from_string('B@1_[(r)]')]),
        SpecificationEquivalenceTestCase(cfs.specification_from_string('A@0_[d/s]'), [], [cfs.specification_from_string('B@1'),
                                                                                          cfs.specification_from_string('B@1_[d]'),
                                                                                          cfs.specification_from_string('B@1_[d/s]'),
                                                                                          cfs.specification_from_string('B@1_[d/s(r)]'),
                                                                                          cfs.specification_from_string('B@1_[d(r)]'),
                                                                                          cfs.specification_from_string('B@1_[(r)]')]),
        SpecificationEquivalenceTestCase(cfs.specification_from_string('A@0_[d/s(r)]'), [], [cfs.specification_from_string('B@1'),
                                                                                             cfs.specification_from_string('B@1_[d]'),
                                                                                             cfs.specification_from_string('B@1_[d/s]'),
                                                                                             cfs.specification_from_string('B@1_[d/s(r)]'),
                                                                                             cfs.specification_from_string('B@1_[d(r)]'),
                                                                                             cfs.specification_from_string('B@1_[(r)]')]),
        SpecificationEquivalenceTestCase(cfs.specification_from_string('A@0_[d(r)]'), [], [cfs.specification_from_string('B@1'),
                                                                                           cfs.specification_from_string('B@1_[d]'),
                                                                                           cfs.specification_from_string('B@1_[d/s]'),
                                                                                           cfs.specification_from_string('B@1_[d/s(r)]'),
                                                                                           cfs.specification_from_string('B@1_[d(r)]'),
                                                                                           cfs.specification_from_string('B@1_[(r)]')]),
        SpecificationEquivalenceTestCase(cfs.specification_from_string('A@0_[(r)]'), [], [cfs.specification_from_string('B@1'),
                                                                                          cfs.specification_from_string('B@1_[d]'),
                                                                                          cfs.specification_from_string('B@1_[d/s]'),
                                                                                          cfs.specification_from_string('B@1_[d/s(r)]'),
                                                                                          cfs.specification_from_string('B@1_[d(r)]'),
                                                                                          cfs.specification_from_string('B@1_[(r)]')])
    ]


