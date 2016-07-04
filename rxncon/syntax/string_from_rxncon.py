from enum import unique
from rxncon.util.utils import OrderedEnum
import typing as tp
import typecheck as tc


# Note: Here be no type annotations, because this would require importing the state, component and reaction modules,
# which would induce a circular dependency.

@unique
class SpecificationSuffix(OrderedEnum):
    mrna = "mRNA"
    gene = "gene"
    protein = ""


@tc.typecheck
def string_from_specification(specification, prefix: OrderedEnum) -> str:
    if str(specification.spec_resolution):
        return '{0}: {1}{2}_[{3}]'.format(type(specification).__name__, create_structured_name(specification), prefix.value, str(specification.spec_resolution), )
    else:
        return '{0}: {1}{2}'.format(type(specification).__name__, create_structured_name(specification), prefix.value)


def string_from_domain_information(domain_resolution):

    if domain_resolution.domain and domain_resolution.subdomain and domain_resolution.residue:
        return '{0}/{1}({2})'.format(domain_resolution.domain, domain_resolution.subdomain, domain_resolution.residue)

    elif domain_resolution.domain and not domain_resolution.subdomain and domain_resolution.residue:
        return '{0}({1})'.format(domain_resolution.domain, domain_resolution.residue)

    elif domain_resolution.domain and domain_resolution.subdomain and not domain_resolution.residue:
        return '{0}/{1}'.format(domain_resolution.domain, domain_resolution.subdomain)

    elif not domain_resolution.domain and not domain_resolution.subdomain and domain_resolution.residue:
        return '({0})'.format(domain_resolution.residue)

    elif domain_resolution.domain and not domain_resolution.subdomain and not domain_resolution.residue:
        return '{0}'.format(domain_resolution.domain)

    elif not domain_resolution.domain and not domain_resolution.subdomain and not domain_resolution.residue:
        return ''

    else:
        raise AssertionError


def string_from_rna_specification(specification):
    return string_from_specification(specification, SpecificationSuffix.mrna)


def string_from_gene_specification(specification):
    return string_from_specification(specification, SpecificationSuffix.gene)


def string_from_protein_specification(specification):
    return string_from_specification(specification, SpecificationSuffix.protein)


@tc.typecheck
def create_structured_name(specification):
    if specification.structure_index is not None:
        return "{0}@{1}".format(specification.name, specification.structure_index)
    else:
        return "{0}".format(specification.name)