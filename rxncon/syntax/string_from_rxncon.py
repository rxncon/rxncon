from enum import unique
from rxncon.util.utils import OrderedEnum
import typing as tp
import typecheck as tc


# Note: Here be no type annotations, because this would require importing the state, component and reaction modules,
# which would induce a circular dependency.

@unique
class SpecificationSuffix(OrderedEnum):
    mrna = "mRNA"
    gene = "Gene"
    protein = ""


@tc.typecheck
def string_from_specification(specification, prefix: OrderedEnum) -> str:
    if str(specification.spec_resolution):
        return '{0}_[{1}]'.format(create_name(specification, prefix), str(specification.spec_resolution))
    else:
        return '{0}'.format(create_name(specification, prefix))


def string_from_domain_resolution(domain_resolution):

    if domain_resolution.domain and domain_resolution.subdomain and domain_resolution.residue:
        return '{0}/{1}({2})'.format(domain_resolution.domain, domain_resolution.subdomain, domain_resolution.residue)

    elif domain_resolution.domain and not domain_resolution.subdomain and domain_resolution.residue:
        return '{0}({1})'.format(domain_resolution.domain, domain_resolution.residue)

    elif domain_resolution.domain and domain_resolution.subdomain and not domain_resolution.residue:
        return '{0}/{1}'.format(domain_resolution.domain, domain_resolution.subdomain)

    elif not domain_resolution.domain and not domain_resolution.subdomain and domain_resolution.residue:
        return '({0})'.format(domain_resolution, domain_resolution.residue)

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
def create_name(specification, prefix: tp.Optional[OrderedEnum]):
    return "{0}{1}".format(specification.name, prefix.value)


def string_from_reaction(reaction) -> str:
    return '{0}_{1}_{2}'.format(reaction.subject, reaction.verb.value, reaction.object)


def string_from_inter_protein_interaction_state(state) -> str:
    return '{0}--{1}'.format(state.first_component, state.second_component)


def string_from_intra_protein_interaction_state(state) -> str:
    return '{0}--[{1}]'.format(state.first_component, state.second_component.domain)


def string_from_covalent_modification_state(state) -> str:
    return '{0}-{{{1}}}'.format(state.substrate, state.modifier.value)


def string_from_translocation_state(state) -> str:
    return '{0}-{{{1}}}'.format(state.substrate, state.compartment.value)


def string_from_synthesis_degradation_state(state) -> str:
    return '{}'.format(state.component)


def string_from_input_state(state) -> str:
    return '{}'.format(state.name)


def string_from_component_state(state) -> str:
    return '{}'.format(state.component)
