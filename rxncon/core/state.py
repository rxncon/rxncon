from enum import unique
import re
import typing as tp

from rxncon.util.utils import OrderedEnum
import rxncon.core.specification as spec

import rxncon.syntax.specification_from_string as sfs


@unique
class StateModifier(OrderedEnum):
    neutral = '0'
    phosphor   = 'p'
    ubiquitin  = 'ub'
    guanosintriphosphat = 'gtp'
    truncated  = 'truncated'


class StateDefinition:
    SPEC_REGEX_GROUPED = '([\\w]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?|[\w]+?|[\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?@[0-9]+?|[\w]+?)'
    SPEC_REGEX_UNGROUPED = '(?:[\\w]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?|[\w]+?|[\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?@[0-9]+?|[\w]+?)'  # substring matched by the group cannot be retrieved after performing a match or referenced later in the pattern.
    def __init__(self, name, representation_def, variables_def):

        self.name, self.representation_def, self.variables_def = name, representation_def, variables_def

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
        variables = {}
        for var, var_def in self.variables_def.items():
            var_regex = self._to_base_regex().replace(var, self.SPEC_REGEX_GROUPED)
            for other_var in self.variables_def.keys():
                if other_var != var:
                    var_regex = var_regex.replace(other_var, self.SPEC_REGEX_UNGROUPED)

            val_str = re.match(var_regex, representation).group(1)
            if self.variables_def[var][0] is StateModifier:
                value = state_modifier_from_string(val_str)
            elif self.variables_def[var][0] is spec.Domain:
                domain, subdomain, residue = sfs.domain_resolution_from_string(val_str)
                value = spec.Domain(domain, subdomain, residue)
            else:
                value = sfs.specification_from_string(val_str)

            variables[var] = value

        return variables

    def representation_from_variables(self, variables):
        representation = self.representation_def
        for var, val in variables.items():
            if isinstance(val, StateModifier):
                representation = representation.replace(var, str(val.value))
            else:
                representation = representation.replace(var, str(val))

        return representation



STATE_DEFINITION = [
    StateDefinition('interaction-state',
                    '$x--$y',
                    {'$x': (spec.Specification, spec.SpecificationResolution.domain),
                     '$y': (spec.Specification, spec.SpecificationResolution.domain)}),

    StateDefinition('self-interaction-state',
                    '$x--[$y]',
                    {'$x': (spec.Specification, spec.SpecificationResolution.domain),
                     '$y': (spec.Domain, spec.SpecificationResolution.domain) }),

    StateDefinition('input-state',
                    '[$x]',
                    {'$x': (spec.Specification, spec.SpecificationResolution.component)}),

    StateDefinition('covalent-modification-state',
                    '$x-{$y}',
                    {'$x': (spec.Specification, spec.SpecificationResolution.residue),
                     '$y': (StateModifier, StateModifier.neutral)}),

    #StateDefinition('component-state',
    #                '$x',
    #                {'$x': (spec.Specification, spec.SpecificationResolution.component)},
    #                [],
#                    ),
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

    def components(self) -> tp.List[spec.Specification]:
        return [value for value in self.variables.values() if isinstance(value, spec.Specification)]

    def modifier(self) -> tp.List[StateModifier]:
        return [value for value in self.variables.values() if isinstance(value, StateModifier)]

    @property
    def neutral_modifier(self) -> StateModifier:
        result = [self.definition.variables_def[var][1] for var, value in self.variables.items()
                  if isinstance(value, StateModifier)]
        assert len(result) <= 1
        if result:
            return result[0]

    def elemental_resolution(self, specification_to_find: spec.Specification) -> spec.SpecificationResolution:
        result = [self.definition.variables_def[var][1] for var, value in self.variables.items()
                  if isinstance(value, spec.Specification) and specification_to_find == value]
        assert len(result) <= 1
        if result:
            return result[0]

    def is_elemental(self) -> bool:
        return all(component.resolution == self.elemental_resolution(component) for component in self.components())

    def is_superset_of(self, other: 'State') -> bool:
        return other.is_subset_of(self)

    def is_subset_of(self, other: 'State') -> bool:
        subspec = False
        if self.definition.name == other.definition.name:# or other.definition.name in self.definition.subset_of_def:
            return self._component_subspecification_validation(other, subspec)
        else:
            return False
        #else:
        #    assert NotImplementedError

    def _component_subspecification_validation(self, other: 'State', subspec: bool):
        for self_component in self.components():
            found_component = False
            for other_component in other.components():
                #check whether the two components are equal
                if self_component.to_component_specification() == other_component.to_component_specification():
                    #check if one is the subspecification and the othter the superspecification
                    found_component = True
                    if self_component.is_subspecification_of(other_component) and other_component.is_superspecification_of(self_component):
                        # otherwise we have to check if the modifier are equal like A_[(r)]-{p} != A_[(r)]-{ub}
                        if self.modifier() == other.modifier():
                            subspec = True
                        else:
                            return False
                    else:
                        return False
            if not found_component:
                return False


        return subspec

def state_modifier_from_string(modifier: str) -> StateModifier:
    return StateModifier(modifier.lower())


def _get_STATE_DEFINITION(superset_name: str) -> StateDefinition:
    for definition in STATE_DEFINITION:
        if superset_name == definition.name:
            return definition

def definition_of_state(representation: str) -> StateDefinition:
    the_definition = None
    for definition in STATE_DEFINITION:

        if definition.matches_representation(representation):
            assert not the_definition
            the_definition = definition

    assert the_definition
    return the_definition


def state_from_string(representation: str) -> State:
    the_definition = definition_of_state(representation)
    variables = the_definition.variables_from_representation(representation)

    return State(the_definition, variables)