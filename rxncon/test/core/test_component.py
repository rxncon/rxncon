import rxncon.syntax.rxncon_from_string as cfs


no_domain = cfs.component_from_string('A')
domain = cfs.component_from_string('A_[b]')
domain_subdomain = cfs.component_from_string('A_[b/c]')
domain_residue = cfs.component_from_string('A_[b(d)]')
domain_subdomain_residue = cfs.component_from_string('A_[b/c(d)]')
residue = cfs.component_from_string('A_[(d)]')


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



