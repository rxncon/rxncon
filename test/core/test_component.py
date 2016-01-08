import core.component as com
import input.shared.from_string as cfs


no_domain = cfs.component_from_string('A')
domain = cfs.component_from_string('A_[b]')
domain_subdomain = cfs.component_from_string('A_[b/c]')
domain_subdomain_residue = cfs.component_from_string('A_[b/c(d)]')
domain_residue = cfs.component_from_string('A_[b/(d)]')
residue = cfs.component_from_string('A_[(d)]')


# def test_component_initialization():
#     assert no_domain.name
#     assert not no_domain.domain
#     assert not no_domain.subdomain
#     assert not no_domain.residue
#
#     assert domain.name
#     assert domain.domain
#     assert not domain.subdomain
#     assert not domain.residue
#
#     assert domain_subdomain.name
#     assert domain_subdomain.domain
#     assert domain_subdomain.subdomain
#     assert not domain_subdomain.residue
#
#     assert domain_subdomain_residue.name
#     assert domain_subdomain_residue.domain
#     assert domain_subdomain_residue.subdomain
#     assert domain_subdomain_residue.residue
#
#     assert domain_residue.name
#     assert domain_residue.domain
#     assert not domain_residue.subdomain
#     assert domain_residue.residue
#
#     assert residue.name
#     assert not residue.domain
#     assert not residue.subdomain
#     assert residue.residue
