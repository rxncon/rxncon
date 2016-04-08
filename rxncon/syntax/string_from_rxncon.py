from enum import Enum, unique
import typing as tp
import typecheck as tc

# Note: Here be no type annotations, because this would require importing the state, component and reaction modules,
# which would induce a circular dependency.

@unique
class SpecificationSuffix(Enum):
    mrna = "mRNA"
    protein = ""

@tc.typecheck
def string_from_specification(specification, prefix: Enum) -> str:
    if specification.domain and specification.subdomain and specification.residue:
        return '{0}_[{1}/{2}({3})]'.format(create_name(specification, prefix), specification.domain, specification.subdomain, specification.residue)

    elif specification.domain and not specification.subdomain and specification.residue:
        return '{0}_[{1}({2})]'.format(create_name(specification, prefix), specification.domain, specification.residue)

    elif specification.domain and specification.subdomain and not specification.residue:
        return '{0}_[{1}/{2}]'.format(create_name(specification, prefix), specification.domain, specification.subdomain)

    elif not specification.domain and not specification.subdomain and specification.residue:
        return '{0}_[({1})]'.format(create_name(specification, prefix), specification.residue)

    elif specification.domain and not specification.subdomain and not specification.residue:
        return '{0}_[{1}]'.format(create_name(specification, prefix), specification.domain)

    elif not specification.domain and not specification.subdomain and not specification.residue:
        return '{0}'.format(create_name(specification, prefix))

    else:
        raise AssertionError

def string_from_rna_specification(specification):
    return string_from_specification(specification, SpecificationSuffix.mrna)

def string_from_protein_specification(specification):
    return string_from_specification(specification, SpecificationSuffix.protein)

@tc.typecheck
def create_name(specification, prefix: tp.Optional[Enum]):
    if prefix == SpecificationSuffix.mrna:
        return "{0}{1}".format(specification.name, SpecificationSuffix.mrna.value)
    else:
        return specification.name


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
