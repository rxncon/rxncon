from abc import ABCMeta, abstractmethod, abstractproperty
from enum import unique
import typecheck as tc
import re
import typing as tp

from rxncon.util.utils import OrderedEnum
import rxncon.core.specification as spec

from rxncon.syntax.rxncon_from_string import specification_from_string


@unique
class StateModifier(OrderedEnum):
    unmodified = '0'
    phosphor   = 'p'
    ubiquitin  = 'ub'
    guanosintriphosphat = 'gtp'
    truncated  = 'truncated'



class StateDefinition():
    SPEC_REGEX_GROUPED = '([a-zA-Z0-9\/\[\]\(\)_]+?)'
    SPEC_REGEX_UNGROUPED = '[a-zA-Z0-9\/\[\]\(\)_]+?'

    def __init__(self, name, representation_def, variables_def, superspecification_def):

        self.name, self.representation_def, \
        self.variables_def, self.superspecification_of_def  = name, representation_def, variables_def, superspecification_def

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'state-definition: name={0}, representation_def={1}'.format(self.name, self.representation_def)

    def matches_representation(self, representation):
        return re.match(self._to_matching_regex(), representation)

    def _to_matching_regex(self):
        regex = '{}'.format(self.representation_def)
        for var in self.variables_def.keys():
            regex = regex.replace(var, self.SPEC_REGEX_GROUPED)
        return '^{}$'.format(regex)

    def variables_from_representation(self, representation):
        assert self.matches_representation(representation)
        variables = { }
        for var, var_def in self.variables_def.items():
            var_regex = '^{}$'.format(self.representation_def.replace(var, self.SPEC_REGEX_GROUPED))
            for other_var in self.variables_def.keys():
                if other_var != var:
                    var_regex = var_regex.replace(other_var, self.SPEC_REGEX_UNGROUPED)

            val_str = re.match(var_regex, representation).group(1)
            if self.variables_def[var] is StateModifier:
                value = state_modifier_from_string(val_str)
            else:
                value = specification_from_string(val_str)

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

# OUTPUT_REGEX = '^\[.+?\]$'
#     INTERACTION_REGEX = '^.+?--.+?$'
#     MODIFIER_REGEX = '.+?-{.+?}'
#     MODIFIER_VALUE_REGEX = '{.+?}'
#     COMPONENT_REGEX = '\w'

STATE_DEFINITION = [
    StateDefinition('interaction-state',
                    '$x--$y',
                    {'$x': spec.Specification,
                     '$y': spec.Specification},
                    ['interaction-state']
                    ),

    StateDefinition('covalent-modification-state',
                    '$x-{$y}',
                    {'$x': spec.Specification,
                     '$y': StateModifier},
                    ['covalent-modification-state']),

    StateDefinition('component-state',
                    '$x',
                    {'$x': spec.Specification},
                    ['interaction-state', 'covalent-modification-state']
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

    def is_superspecification_of(self, other) -> bool:
        superspec = False
        if other.definition.name == self.definition.name or other.definition.name in self.definition.superspecification_of_def:
            for self_component in self.components():
                for other_component in other.components():
                    if self_component.to_component_specification() == other_component.to_component_specification():
                            if self_component.is_superspecification_of(other_component) \
                                    and other_component.is_subspecification_of(self_component):
                                if other.definition.name in self.definition.superspecification_of_def:
                                    superspec = True
                                elif self.modifier() == other.modifier():
                                    superspec = True
                                else:
                                    return False
                            else:
                                return False
        return superspec

    def is_subspecification_of(self, other) -> bool:
        subspec = False
        if other.definition.name == self.definition.name or self.definition.name in other.definition.superspecification_of_def:
            for self_component in self.components():
                for other_component in other.components():
                    if self_component.to_component_specification() == other_component.to_component_specification():
                        if self_component.is_subspecification_of(other_component) \
                                and other_component.is_superspecification_of(self_component):
                            if self.definition.name in other.definition.superspecification_of_def:
                                subspec = True
                            elif self.modifier() == other.modifier():
                                subspec = True
                            else:
                                return False
                        else:
                            return False
        return subspec


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
