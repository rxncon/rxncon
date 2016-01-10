import re
from typing import Tuple

import rxncon.core.component as com
import rxncon.core.error as err
import rxncon.core.reaction as rxn
import rxncon.core.state as sta


def component_from_string(component_string: str) -> com.Component:
    DOMAIN_DELIMITER = '_'

    DOMAIN_SUBDOMAIN_RESIDUE_REGEX = '^[\w:-]+\/[\w:-]+\([\w:-]+\)$'
    DOMAIN_RESIDUE_REGEX = '^[\w:-]+\([\w:-]+\)$'
    DOMAIN_SUBDOMAIN_REGEX = '^[\w:-]+\/[\w:-]+$'
    RESIDUE_REGEX = '^\([\w:-]+\)$'
    DOMAIN_REGEX = '^[\w:-]+$'

    items = component_string.split(DOMAIN_DELIMITER, maxsplit=1)

    if len(items) == 1:
        return com.Component(items[0], None, None, None)

    elif len(items) == 2:
        name = items[0]
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
            raise AssertionError

        return com.Component(name, domain, subdomain, residue)

    else:
        raise err.RxnConParseError('Could not parse component string '.format(component_string))


def reaction_from_string(reaction_string: str) -> rxn.Reaction:
    INPUT_REGEX = '^\[.+?\]$'

    if re.match(INPUT_REGEX, reaction_string):
        return rxn.OutputReaction(reaction_string)

    subject, verb, object = _reaction_string_to_subject_verb_object_strings(reaction_string)

    subject = component_from_string(subject)
    verb = rxn.Verb(verb)
    object = component_from_string(object)

    try:
        reaction_definition = rxn.VERB_REACTION_TABLE[verb]

    except KeyError:
        raise err.RxnConParseError('The verb {} could not be found in the reaction lookup table'.format(verb))

    category, directionality, influence, isomerism, modifier = reaction_definition

    return rxn.Reaction(subject, verb, object, category, directionality, influence, isomerism, modifier)


def _reaction_string_to_subject_verb_object_strings(reaction_string: str) -> Tuple[str, str, str]:
    known_verbs = [v.value for v in rxn.Verb]
    delimiters = ['_' + verb + '_' for verb in known_verbs]

    matches = 0
    subject_string, verb_string, object_string = None, None, None

    for delimiter in delimiters:
        splitted = reaction_string.split(delimiter)

        if len(splitted) > 1:
            matches += 1

            if len(splitted) != 2:
                raise err.RxnConParseError('Incorrect S-V-O in reaction string {}'.format(reaction_string))

            subject_string, verb_string, object_string = splitted[0], delimiter.strip('_'), splitted[1]

    if matches == 0:
        raise err.RxnConParseError('Unknown verb in reaction string {}'.format(reaction_string))

    elif matches > 1:
        raise err.RxnConParseError('Ambiguous verb in reaction string {}'.format(reaction_string))

    return subject_string, verb_string, object_string


def state_from_string(state_string: str) -> sta.State:
    OUTPUT_REGEX = '^\[.+?\]$'
    INTERACTION_REGEX = '^.+?--.+?$'
    MODIFIER_REGEX = '.+?-{.+?}'
    MODIFIER_VALUE_REGEX = '{.+?}'

    assert isinstance(state_string, str)

    if re.match(INTERACTION_REGEX, state_string):
        return _interaction_state_from_string(state_string)

    elif re.match(MODIFIER_REGEX, state_string):
        modifier_string = re.search(MODIFIER_VALUE_REGEX, state_string).group(0).strip('{}')

        # @todo Maybe this is a bit dirty.
        try:
            modifier = sta.StateModifier(modifier_string)
            return _covalent_modification_state_from_string(state_string, modifier)

        except ValueError:
            return _translocation_state_from_string(state_string)

    elif re.match(OUTPUT_REGEX, state_string):
        return sta.InputState(state_string)

    else:
        raise err.RxnConParseError('Could not parse state string {} into State'.format(state_string))


def _interaction_state_from_string(state_string: str) -> sta.InteractionState:
    component_strings = state_string.split('--')
    first_component = component_from_string(component_strings[0])
    second_component = component_from_string(component_strings[1])

    return sta.InteractionState(first_component, second_component)


def _covalent_modification_state_from_string(state_string: str, modifier: sta.StateModifier) -> sta.CovalentModificationState:
    substrate_string = state_string.split('-{')[0]
    substrate = component_from_string(substrate_string)

    return sta.CovalentModificationState(substrate, modifier)


def _translocation_state_from_string(state_string) -> sta.TranslocationState:
    # @todo Implement this.
    pass