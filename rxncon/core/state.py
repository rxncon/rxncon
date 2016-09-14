import re
from enum import unique
from collections import OrderedDict
from typing import List, Dict, Optional, Any
from typecheck import typecheck

from rxncon.util.utils import OrderedEnum, members
from rxncon.core.spec import Spec, MolSpec, EmptyMolSpec, Locus, LocusResolution, locus_from_string, mol_spec_from_string, spec_from_string

FULLY_NEUTRAL_STATE_REPR = '0'

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
    def __init__(self, name: str, repr_def: str, variables_def: Dict[str, Any], target_spec_def: str,
                 neutral_states_def: List[str]):
        self.name, self.repr_def, self.variables_def, self.target_spec_def, self.neutral_states_def = \
            name, repr_def, variables_def, target_spec_def, neutral_states_def

    def __str__(self) -> str:
        return 'StateDef: name={0}, representation_def={1}'.format(self.name, self.repr_def)

    @typecheck
    def matches_repr(self, repr: str) -> bool:
        return True if re.match(self._to_matching_regex(), repr) else False

    @typecheck
    def variables_from_repr(self, representation: str) -> Dict[str, Any]:
        assert self.matches_repr(representation)

        variables = {}
        for var, var_def in self.variables_def.items():
            var_regex = self._to_base_regex().replace(var, self.SPEC_REGEX_GROUPED)
            for other_var in self.variables_def.keys():
                if other_var != var:
                    var_regex = var_regex.replace(other_var, self.SPEC_REGEX_UNGROUPED)

            val_str = re.match(var_regex, representation).group(1)

            if self.variables_def[var][0] is StateModifier:
                value = state_modifier_from_string(val_str)
            elif self.variables_def[var][0] is Locus:
                value = locus_from_string(val_str)
            elif self.variables_def[var][0] is MolSpec:
                value = mol_spec_from_string(val_str)
            else:
                raise Exception('Unknown variable type {}'.format(self.variables_def[var][0]))

            variables[var] = value

        return variables

    @typecheck
    def repr_from_variables(self, variables: Dict[str, Any]) -> str:
        representation = self.repr_def
        for var, val in variables.items():
            if isinstance(val, StateModifier):
                representation = representation.replace(var, str(val.value))
            else:
                representation = representation.replace(var, str(val))

        return representation

    @typecheck
    def target_from_variables(self, variables: Dict[str, Any]) -> Spec:
        target_str = self.target_spec_def
        for var_symbol, value in variables.items():
            for method in members(value):
                symbol_with_method = '{0}.{1}'.format(var_symbol, method)
                if symbol_with_method in target_str:
                    target_str = target_str.replace(symbol_with_method, str(getattr(value, method)()))

            if var_symbol in target_str:
                target_str = target_str.replace(var_symbol, str(value))

        return spec_from_string(target_str)

    @typecheck
    def neutral_states_from_variables(self, variables: Dict[str, Any]) -> List['State']:
        states = []

        for state_def in self.neutral_states_def:
            state_str = state_def
            for var, val in variables.items():
                state_str = state_str.replace(var, str(val))

            try:
                the_state = state_from_string(state_str)
            except SyntaxError:
                the_state = None

            if the_state:
                states.append(the_state)

        return states

    @typecheck
    def validate_variables(self, variables: Dict[str, Any]):
        if all(isinstance(x, EmptyMolSpec) for x in variables.values() if isinstance(x, Spec)):
            raise SyntaxError('All MolSpec are EmptyMolSpec.')

        for var, val in variables.items():
            if isinstance(val, MolSpec):
                if val.resolution > self.variables_def[var][1]:
                    raise SyntaxError('Resolution too high.')

    def _to_base_regex(self) -> str:
        return '^{}$'.format(self.repr_def
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
        'interaction-state',                          # name
        '$x--$y',                                     # representation definition
        {
            '$x': [MolSpec, LocusResolution.domain],  # variables definition
            '$y': [MolSpec, LocusResolution.domain]
        },
        '$x~$y',                                      # target spec
        ['$x--0', '$y--0']                            # neutral states
    ),
    StateDef(
        'self-interaction-state',
        '$x--[$y]',
        {
            '$x': [MolSpec, LocusResolution.domain],
            '$y': [Locus, LocusResolution.domain]
        },
        '$x~$x.to_component_spec_[$y]',
        ['$x--0', '$y--0']
    ),
    StateDef(
        'input-state',
        '[$x]',
        {
            '$x': [MolSpec, LocusResolution.component]
        },
        'global',
        []
    ),
    StateDef(
        'covalent-modification-state',
        '$x-{$y}',
        {
            '$x': [MolSpec, LocusResolution.residue],
            '$y': [StateModifier]
        },
        '$x',
        ['$x-{0}']
    )
]


class State:
    @typecheck
    def __init__(self, definition: StateDef, variables: Dict[str, Any]):
        # @todo fix ordering of variables
        self.definition = definition
        self.state_defs = STATE_DEFS
        self.variables = OrderedDict((k, v) for k, v in sorted(variables.items()))

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return self.definition.repr_from_variables(self.variables)

    @typecheck
    def __eq__(self, other: 'State') -> bool:
        return self.definition == other.definition and self.variables == other.variables

    @typecheck
    def __lt__(self, other: 'State'):
        return str(self) < str(other)

    @property
    @typecheck
    def target(self) -> Spec:
        return self.definition.target_from_variables(self.variables)

    @property
    @typecheck
    def is_elemental(self) -> bool:
        return all(spec.has_resolution(elemental_resolution) for spec, elemental_resolution
                   in zip(self.mol_specs, self._elemental_resolutions))

    @typecheck
    def is_superset_of(self, other: 'State') -> bool:
        return other.is_subset_of(self)

    @typecheck
    def is_subset_of(self, other: 'State') -> bool:
        if self.definition == other.definition:
            return all(x.is_subspec_of(y) for x, y in zip(self.mol_specs, other.mol_specs)) and \
                   all(x == y for x, y in zip(self._non_mol_spec_props, other._non_mol_spec_props))
        else:
            return False

    @property
    @typecheck
    def is_neutral(self) -> bool:
        return len(self.neutral_states) == 1 and self == self.neutral_states[0]

    @property
    @typecheck
    def neutral_states(self) -> List['State']:
        return self.definition.neutral_states_from_variables(self.variables)

    @property
    @typecheck
    def mol_specs(self) -> List[MolSpec]:
        return [x for x in self.variables.values() if isinstance(x, MolSpec) and not isinstance(x, EmptyMolSpec)]

    @property
    @typecheck
    def components(self) -> List[MolSpec]:
        return [x.to_component_spec() for x in self.mol_specs]

    @property
    @typecheck
    def _elemental_resolutions(self) -> List[LocusResolution]:
        return [self.definition.variables_def[var][1] for var, value in self.variables.items()
                if isinstance(value, MolSpec)]

    @property
    def _non_mol_spec_props(self):
        return [x for x in self.variables.values() if not isinstance(x, MolSpec)]


class FullyNeutralState(State):
    def __init__(self):
        pass

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return 'fully-neutral-state'

    def __repr__(self) -> str:
        return str(self)

    def __lt__(self, other: 'State'):
        return True

    def __eq__(self, other: 'State') -> bool:
        return isinstance(other, FullyNeutralState)

    def is_subset_of(self, other: 'State') -> bool:
        raise AssertionError

    @property
    def neutral_states(self) -> List['State']:
        raise AssertionError

    @property
    def components(self) -> List[MolSpec]:
        return []

    @property
    def target(self) -> Spec:
        raise AssertionError

    @property
    def is_elemental(self) -> bool:
        raise AssertionError

    @property
    def is_neutral(self) -> bool:
        raise AssertionError

    def is_superset_of(self, other: 'State') -> bool:
        raise AssertionError


@typecheck
def state_modifier_from_string(modifier: str) -> StateModifier:
    return StateModifier(modifier.lower())

@typecheck
def state_from_string(repr: str) -> State:
    if repr == FULLY_NEUTRAL_STATE_REPR:
        return FullyNeutralState()

    state_def = next((x for x in STATE_DEFS if x.matches_repr(repr)), None)

    if not state_def:
        raise SyntaxError('Could not match State {} with definition'.format(repr))

    variables = state_def.variables_from_repr(repr)
    # variables = state_def.order_variables(variables)

    state_def.validate_variables(variables)

    return State(state_def, variables)
