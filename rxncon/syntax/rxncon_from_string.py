import re
from typing import Union
from collections import OrderedDict
from typecheck import typecheck

from rxncon.core.specification import DnaSpec, Spec, RnaSpec, ProteinSpec
from rxncon.core.state import State, InteractionState, TranslocationState, StateModifier, CovalentModificationState, \
    GlobalQuantityState, SelfInteractionState, ComponentState
from enum import unique
from rxncon.util.utils import OrderedEnum


@unique
class SpecificationSuffix(OrderedEnum):
    mrna = "mRNA"
    gene = "Gene"
    protein = ""


mapping_suffix_to_specification = OrderedDict([(SpecificationSuffix.mrna, RnaSpec),
                                               (SpecificationSuffix.protein, ProteinSpec),
                                               (SpecificationSuffix.gene, DnaSpec)])


def create_specification_from_name_suffix(name, domain, subdomain, residue, struct_index=None):
    for suffix in mapping_suffix_to_specification:
        if name.endswith(suffix.value):
            return mapping_suffix_to_specification[suffix](name, domain, subdomain, residue, struct_index)


@typecheck
def specification_from_string(specification_string: str) -> Spec:
    DOMAIN_DELIMITER = '_'
    INDEX_DELIMITER = '@'

    DOMAIN_SUBDOMAIN_RESIDUE_REGEX = '^[\w:-]+\/[\w:-]+\([\w:-]+\)$'
    DOMAIN_RESIDUE_REGEX = '^[\w:-]+\([\w:-]+\)$'
    DOMAIN_SUBDOMAIN_REGEX = '^[\w:-]+\/[\w:-]+$'
    RESIDUE_REGEX = '^\([\w:-]+\)$'
    DOMAIN_REGEX = '^[\w:-]+$'

    items = specification_string.split(DOMAIN_DELIMITER, maxsplit=1)

    name_items = items[0].split(INDEX_DELIMITER)
    if len(name_items) == 1:
        name = name_items[0]
        index = None
    elif len(name_items) == 2:
        name = name_items[0]
        index = name_items[1]
    else:
        raise SyntaxError

    if len(items) == 1:
        return create_specification_from_name_suffix(name, None, None, None, struct_index=index)

    elif len(items) == 2:
        full_domain_string = items[1].strip('[]')

        if re.match(DOMAIN_SUBDOMAIN_RESIDUE_REGEX, full_domain_string):
            domain = full_domain_string.split('/')[0]
            subdomain = full_domain_string.split('/')[1].split('(')[0]
            residue = full_domain_string.split('/')[1].split('(')[1].strip(')')

        elif re.match(DOMAIN_RESIDUE_REGEX, full_domain_string):
            domain = full_domain_string.split('(')[0]
            subdomain = None
            residue = full_domain_string.split('(')[1].strip(')')

        elif re.match(DOMAIN_SUBDOMAIN_REGEX, full_domain_string):
            domain = full_domain_string.split('/')[0]
            subdomain = full_domain_string.split('/')[1]
            residue = None

        elif re.match(RESIDUE_REGEX, full_domain_string):
            domain = None
            subdomain = None
            residue = full_domain_string.strip('()')

        elif re.match(DOMAIN_REGEX, full_domain_string):
            domain = full_domain_string
            subdomain = None
            residue = None

        else:
            raise SyntaxError('Could not parse specification string {}'.format(specification_string))

        return create_specification_from_name_suffix(name, domain, subdomain, residue, struct_index=index)

    else:
        raise SyntaxError('Could not parse specification string {}'.format(specification_string))


@typecheck
def state_from_string(state_string: str) -> State:
    OUTPUT_REGEX = '^\[.+?\]$'
    INTERACTION_REGEX = '^.+?--.+?$'
    MODIFIER_REGEX = '.+?-{.+?}'
    MODIFIER_VALUE_REGEX = '{.+?}'
    COMPONENT_REGEX = '\w'

    assert isinstance(state_string, str)

    if re.match(INTERACTION_REGEX, state_string):
        return _interaction_state_from_string(state_string)

    elif re.match(MODIFIER_REGEX, state_string):
        modifier_string = re.search(MODIFIER_VALUE_REGEX, state_string).group(0).strip('{}').lower()
        try:
            modifier = StateModifier(modifier_string)
            return _covalent_modification_state_from_string(state_string, modifier)

        except ValueError:
            return _translocation_state_from_string(state_string)

    elif re.match(OUTPUT_REGEX, state_string):
        return GlobalQuantityState(state_string)
    elif re.match(COMPONENT_REGEX, state_string):
        return _component_state_from_string(state_string)
    else:
        raise SyntaxError('Could not parse state string {} into State'.format(state_string))


@typecheck
def _interaction_state_from_string(state_string: str) -> Union[InteractionState,
                                                               SelfInteractionState]:
    component_strings = state_string.split('--')

    if component_strings[1].startswith('['):
        first_component = specification_from_string(component_strings[0])
        second_component = specification_from_string("{0}_{1}".format(first_component.name, component_strings[1]))

        return SelfInteractionState(first_component, second_component)

    else:
        first_component = specification_from_string(component_strings[0])
        second_component = specification_from_string(component_strings[1])

        return InteractionState(first_component, second_component)


@typecheck
def _covalent_modification_state_from_string(state_string: str, modifier: StateModifier) -> CovalentModificationState:
    substrate_string = state_string.split('-{')[0]
    substrate = specification_from_string(substrate_string)

    return CovalentModificationState(substrate, modifier)


def _component_state_from_string(state_string: str):
    return ComponentState(specification_from_string(state_string))


@typecheck
def _translocation_state_from_string(state_string) -> TranslocationState:
    # @todo Implement this.
    pass

