from abc import ABCMeta, abstractmethod, abstractproperty
from enum import unique
import typecheck as tc
import re
import typing as tp

from rxncon.util.utils import OrderedEnum
import rxncon.core.specification as spec

import rxncon.syntax.specification_from_string as sfs


@unique
class StateModifier(OrderedEnum):
    unmodified = '0'
    phosphor   = 'p'
    ubiquitin  = 'ub'
    guanosintriphosphat = 'gtp'
    truncated  = 'truncated'


class StateDefinition:
    #SPEC_REGEX_GROUPED = '([\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?@[0-9]+?|[\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?\.[a-zA-Z]+?|[\\w]+?@[0-9]+?\.[a-zA-Z]+?)'  # we allow something like A.gene
    #SPEC_REGEX_UNGROUPED = '(?:[\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?@[0-9]+?|[\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?\.[a-zA-Z]+?|[\\w]+?@[0-9]+?\.[a-zA-Z]+?)'  # substring matched by the group cannot be retrieved after performing a match or referenced later in the pattern.
    SPEC_REGEX_GROUPED = '([\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?@[0-9]+?|[\w]+?)'
    SPEC_REGEX_UNGROUPED = '(?:[\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?@[0-9]+?|[\w]+?)'  # substring matched by the group cannot be retrieved after performing a match or referenced later in the pattern.
    def __init__(self, name, representation_def, variables_def, subspecification_def):

        self.name, self.representation_def, \
        self.variables_def, self.subspecification_of_def = name, representation_def, variables_def, subspecification_def

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'state-definition: name={0}, representation_def={1}'.format(self.name, self.representation_def)

    def matches_representation(self, representation):
        return re.match(self._to_matching_regex(), representation)

    def _to_base_regex(self):
        regex = '^{}$'.format(self.representation_def.replace('+', '\+'))
        regex = regex.replace('[', '\[')
        regex = regex.replace(']', '\]')
        return regex

    def _to_matching_regex(self):
        regex = self._to_base_regex()
        for var in self.variables_def.keys():
            regex = regex.replace(var, self.SPEC_REGEX_GROUPED)
        return regex

    def variables_from_representation(self, representation):
        assert self.matches_representation(representation)
        variables = { }
        for var, var_def in self.variables_def.items():
            var_regex = self._to_base_regex().replace(var, self.SPEC_REGEX_GROUPED)
            for other_var in self.variables_def.keys():
                if other_var != var:
                    var_regex = var_regex.replace(other_var, self.SPEC_REGEX_UNGROUPED)

            val_str = re.match(var_regex, representation).group(1)
            if self.variables_def[var] is StateModifier:
                value = state_modifier_from_string(val_str)
            elif self.variables_def[var] is spec.DomainResolution:
                domain, subdomain, residue = sfs.domain_resolution_from_string(val_str)
                value = spec.DomainResolution(domain, subdomain, residue)
            else:
                value = sfs.specification_from_string(val_str)

            variables[var] = value

        return variables

    def representation_from_variables(self, variables):
        representation = self.representation_def
        for var, val in variables.items():
            if val is StateModifier:
                representation = representation.replace(var, str(val.value))
            else:
                representation = representation.replace(var, str(val))

        return representation

    def superspec_from_definition(self):
        pass


def state_modifier_from_string(modifier: str):
    return StateModifier(modifier.lower())


STATE_DEFINITION = [
    StateDefinition('interaction-state',
                    '$x--$y',
                    {'$x': spec.Specification,
                     '$y': spec.Specification,},
                    ['component-state']
                    ),

    StateDefinition('self-interaction-state',
                    '$x--[$y]',
                    {'$x': spec.Specification,
                     '$y': spec.DomainResolution},
                    ['component-state']),


    StateDefinition('input-state',
                    '[$x]',
                    {'$x': spec.Specification},
                    []),

    StateDefinition('covalent-modification-state',
                    '$x-{$y}',
                    {'$x': spec.Specification,
                     '$y': StateModifier},
                    ['component-state']),

    StateDefinition('component-state',
                    '$x',
                    {'$x': spec.Specification},
                    [],
                    ),
]


class State():
    def __init__(self, definition: StateDefinition, variables):
        self.definition, self.variables = definition, variables

    def __repr__(self):
        return str(self)

    def __str__(self):
        return self.definition.representation_from_variables(self.variables)

    def __eq__(self, other):
        return self.definition == other.definition and self.variables == other.variables

    def components(self):
        return [value for value in self.variables.values() if isinstance(value, spec.Specification)]

    def modifier(self):
        return [value for value in self.variables.values() if isinstance(value, StateModifier)]

    def is_superset_of(self, other) -> bool:
        # A -> B -> C is C subset of A? True subsets of subsets
        # todo: multi-step inheritance
        # todo: neutral modifier
        return other.is_subset_of(self)

    def is_subset_of(self, other) -> bool:

        subspec = False
        if self.definition.name == other.definition.name or other.definition.name in self.definition.subspecification_of_def:
                for self_component in self.components():
                    for other_component in other.components():
                        if self_component.is_subspecification_of(other_component) and other_component.is_superspecification_of(self_component):
                            subspec = True
                        elif self.modifier() == other.modifier():
                            subspec = True
                        else:
                            return False
        return subspec

def representation_to_general_template(representation):
    return re.sub('@[0-9]+?', '', representation)

def definition_of_state(representation: str):
    the_definition = None
    for definition in STATE_DEFINITION:

        if definition.matches_representation(representation):
            assert not the_definition
            the_definition = definition

    assert the_definition
    return the_definition


def state_from_string(representation: str):
    the_definition = definition_of_state(representation)
    variables = the_definition.variables_from_representation(representation)

    return State(the_definition, variables)