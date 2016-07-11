import re
from enum import unique
from typing import List, Dict
from typecheck import typecheck

from rxncon.util.utils import OrderedEnum


@unique
class StateModifier(OrderedEnum):
    neutral   = '0'
    phosphor  = 'p'
    ubiquitin = 'ub'
    gtp       = 'gtp'
    truncated = 'truncated'


class StateDef:
    SPEC_REGEX_GROUPED = '([\\w]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?|[\w]+?|[\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?@[0-9]+?|[\w]+?)'
    SPEC_REGEX_UNGROUPED = '(?:[\\w]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?|[\w]+?|[\\w]+?@[0-9]+?_[\\w\\/\\[\\]\\(\\)]+?|[\w]+?@[0-9]+?|[\w]+?)'

    @typecheck
    def __init__(self, name: str, representation_def: str, variables_def: Dict):
        self.name, self.representation_def, self.variables_def = name, representation_def, variables_def

    def __str__(self) -> str:
        return 'StateDef: name={0}, representation_def={1}'.format(self.name, self.representation_def)

    def matches_representation(self, representation: str) -> bool:
        return True if re.match(self._to_matching_regex(), representation) else False

    def variables_from_representation(self, representation: str) -> Dict:
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
            elif self.variables_def[var][0] is Domain:
                domain, subdomain, residue = _locus_items_from_string(val_str)
                value = Domain(domain, subdomain, residue)
            else:
                value = spec_from_string(val_str)

            variables[var] = value

        return variables

    def representation_from_variables(self, variables: Dict) -> str:
        representation = self.representation_def
        for var, val in variables.items():
            if val is StateModifier:
                representation = representation.replace(var, str(val.value))
            else:
                representation = representation.replace(var, str(val))

        return representation

    def _to_base_regex(self) -> str:
        return '^{}$'.format(self.representation_def
                             .replace('+', '\+')
                             .replace('[', '\[')
                             .replace(']', '\]'))

    def _to_matching_regex(self) -> str:
        regex = self._to_base_regex()
        for var in self.variables_def.keys():
            regex = regex.replace(var, self.SPEC_REGEX_GROUPED)
        return regex


STATE_DEFS = [
    StateDef(
        'interaction-state',
        '$x--$y',
        {
            '$x': (Spec, SpecificationResolution.domain),
            '$y': (Spec, SpecificationResolution.domain)
        }
    ),
    StateDef(
        'self-interaction-state',
        '$x--[$y]',
        {
            '$x': (Spec, SpecificationResolution.domain),
            '$y': (Domain, SpecificationResolution.domain)
        }
    ),
    StateDef(
        'input-state',
        '[$x]',
        {
            '$x': (Spec, SpecificationResolution.component)
        }
    ),
    StateDef(
        'covalent-modification-state',
        '$x-{$y}',
        {
            '$x': (Spec, SpecificationResolution.residue),
            '$y': (StateModifier, StateModifier.neutral)
        }
    ),
    StateDef(
        'component-state',
        '$x',
        {
            '$x': (Spec, SpecificationResolution.component)
        }
    ),
]


class State:
    def __init__(self, state_defs: List[StateDef], definition: StateDef, variables: Dict):
        self.state_defs, self.definition, self.variables = state_defs, definition, variables

    def __repr__(self):
        return str(self)

    def __str__(self):
        return self.definition.representation_from_variables(self.variables)

    def __eq__(self, other):
        return self.definition == other.definition and self.variables == other.variables

    @property
    def components(self) -> List[Spec]:
        return [value for value in self.variables.values() if isinstance(value, Spec)]

    @property
    def modifier(self) -> StateModifier:
        return next((value for value in self.variables.values() if isinstance(value, StateModifier)), None)

    @property
    def elemental_resolutions(self) -> List[SpecificationResolution]:
        return [self.definition.variables_def[var][1] for var, value in self.variables.items()
                if isinstance(value, Spec)]

    def is_elemental(self) -> bool:
        return all(component.has_resolution(elemental_resolution) for component, elemental_resolution
                   in zip(self.components, self.elemental_resolutions))

    def is_superset_of(self, other: 'State') -> bool:
        return other.is_subset_of(self)

    def is_subset_of(self, other: 'State') -> bool:
        if self.definition == other.definition:
            return all(x.is_subspec_of(y) for x, y in zip(self.variables.values(), other.variables.values())
                       if isinstance(x, Spec))
        else:
            return False

def state_modifier_from_string(modifier: str) -> StateModifier:
    return StateModifier(modifier.lower())


def state_from_string(state_defs: List[StateDef], representation: str) -> State:
    the_definition = next((state_def for state_def in state_defs
                           if state_def.matches_representation(representation)), None)

    assert the_definition, 'Could not match reaction {} with definition'.format(representation)
    variables = the_definition.variables_from_representation(representation)

    return State(state_defs, the_definition, variables)
