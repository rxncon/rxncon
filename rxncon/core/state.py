import re
from enum import unique
from collections import OrderedDict
from typing import List, Dict, Optional, Any
from copy import deepcopy

from rxncon.util.utils import OrderedEnum, members
from rxncon.core.spec import Spec, Locus, LocusResolution, locus_from_str, spec_from_str


FULLY_NEUTRAL_STATE = '0'
GLOBAL_STATE_NAME   = 'GLOBALSTATE'
SPEC_REGEX          = '([A-Za-z][A-Za-z0-9]*(?:@[\d]+)*(?:_\[[\w\/\(\)]+\])*)'
STR_REGEX           = '[\w]+'
LOCUS_REGEX         = '[\w\/\(\)]+'
MATCHING_REGEX      = '([\w@\(\)\[\]]+)'
NON_MATCHING_REGEX  = '(?:[\w@\(\)\[\]]+)'

@unique
class StateModifier(OrderedEnum):
    neutral   = '0'
    phosphor  = 'p'
    ubiquitin = 'ub'
    gtp       = 'gtp'
    truncated = 'truncated'
    cytosol   = 'cytosol'  # @note: hack for localisation.
    nucleus   = 'nucleus'  # @note: hack for localisation.
    out       = 'out'      # @note: hack for SPS
    sat 	  = 'sat'
    dup 	  = 'dup'
    spb 	  = 'spb'
    sep 	  = 'sep'
    lic 	  = 'lic'
    repinit	  = 'repinit'
    rep 	  = 'rep'
    seg 	  = 'seg'
    pol 	  = 'pol'
    bud	  	  = 'bud'
    iso 	  = 'iso'
    cyt 	  = 'cyt'
    myr       = 'myr'
    ac        = 'ac'

def state_modifier_from_str(modifier_str: str) -> StateModifier:
    try:
        return StateModifier(modifier_str.lower())
    except ValueError:
        valid_modifiers = [modifier.value for modifier in StateModifier.__members__.values()]  # type: ignore
        raise ValueError('Invalid StateModifier {}, valid modifiers are {}'
                         .format(modifier_str, ', '.join(valid_modifiers)))


TYPE_TO_REGEX = {
    Spec:          SPEC_REGEX,
    str:           STR_REGEX,
    StateModifier: STR_REGEX,
    Locus:         LOCUS_REGEX
}

TYPE_TO_CONSTRUCTOR = {
    Spec:          spec_from_str,
    str:           lambda x: x.strip(),
    StateModifier: state_modifier_from_str,
    Locus:         locus_from_str
}


class StateDef:
    def __init__(self, name: str, repr_def: str, vars_def: Dict[str, Any], neutral_states_def: List[str]) -> None:
        self.name, self.repr_def, self.vars_def, self.neutral_states_def = name, repr_def, vars_def, neutral_states_def

    def __hash__(self):
        return hash(str(self))

    def __str__(self) -> str:
        return 'StateDef<{0}>'.format(self.name)

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, StateDef):
            return NotImplemented
        return self.name == other.name and self.repr_def == other.repr_def and self.vars_def == other.vars_def and \
            self.neutral_states_def == other.neutral_states_def

    def matches_repr(self, repr: str) -> bool:
        return True if re.match(self._to_matching_regex(), repr) else False

    def vars_from_repr(self, repr: str) -> Dict[str, Any]:
        assert self.matches_repr(repr)

        variables = {}
        for var, var_def in self.vars_def.items():
            var_regex = self._to_base_regex().replace(var, MATCHING_REGEX)
            for other_var in self.vars_def.keys():
                if other_var != var:
                    var_regex = var_regex.replace(other_var, NON_MATCHING_REGEX)

            val_str = re.match(var_regex, repr).group(1)
            variables[var] = TYPE_TO_CONSTRUCTOR[var_def[0]](val_str)

        return variables

    def repr_from_vars(self, vars: Dict[str, Any]) -> str:
        repr = self.repr_def
        for var, val in vars.items():
            if isinstance(val, StateModifier):
                repr = repr.replace(var, str(val.value))
            else:
                repr = repr.replace(var, str(val))

        return repr

    def neutral_states_from_vars(self, variables: Dict[str, Any]) -> List['State']:
        states = []

        for state_def in self.neutral_states_def:
            states.append(state_from_str(self._fill_vals_into_vars(state_def, variables)))

        return states

    def validate_vars(self, vars: Dict[str, Any]):
        for var, val in vars.items():
            assert isinstance(val, self.vars_def[var][0])

            if isinstance(val, Spec):
                if val.resolution > self.vars_def[var][1]:
                    raise SyntaxError('Resolution too high.')

    def _to_base_regex(self) -> str:
        return '^{}$'.format(self.repr_def
                             .replace('+', '\+')
                             .replace('[', '\[')
                             .replace(']', '\]'))

    def _fill_vals_into_vars(self, str_with_vars: str, variables: Dict[str, Any]) -> str:
        for var, val in variables.items():
            for method in members(val):
                symbol_with_method = '{0}.{1}'.format(var, method)
                if symbol_with_method in str_with_vars:
                    str_with_vars = str_with_vars.replace(symbol_with_method, str(getattr(val, method)()))

            if var in str_with_vars:
                str_with_vars = str_with_vars.replace(var, str(val))

        return str_with_vars

    def _to_matching_regex(self) -> str:
        regex = self._to_base_regex()
        for var, var_def in self.vars_def.items():
            regex = regex.replace(var, TYPE_TO_REGEX[var_def[0]])
        return regex


STATE_DEFS = [
    StateDef(
        'interaction-state',                          # name
        '$x--$y',                                     # representation definition
        {
            '$x': [Spec, LocusResolution.domain],     # variables definition
            '$y': [Spec, LocusResolution.domain]
        },
        ['$x--0', '$y--0']                            # neutral states
    ),
    StateDef(
        'empty-binding-domain',
        '$x--0',
        {
            '$x': [Spec, LocusResolution.domain]
        },
        ['$x--0']
    ),
    StateDef(
        'self-interaction-state',
        '$x--[$y]',
        {
            '$x': [Spec, LocusResolution.domain],
            '$y': [Locus, LocusResolution.domain]
        },
        ['$x--0', '$x.to_component_spec_[$y]--0']
    ),
    StateDef(
        GLOBAL_STATE_NAME,
        '[$x]',
        {
            '$x': [str]
        },
        []
    ),
    StateDef(
        'covalent-modification-state',
        '$x-{$y}',
        {
            '$x': [Spec, LocusResolution.residue],
            '$y': [StateModifier]
        },
        ['$x-{0}']
    )
]


class State:
    def __init__(self, definition: StateDef, vars: Dict[str, Any]) -> None:
        self.definition = definition
        self.state_defs = STATE_DEFS
        self.vars = OrderedDict((k, v) for k, v in sorted(vars.items()))

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return self.definition.repr_from_vars(self.vars)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented
        return self.definition == other.definition and self.vars == other.vars

    def __lt__(self, other: 'State'):
        return str(self) < str(other)

    def __getitem__(self, item):
        return self.vars[item]

    @property
    def name(self):
        return self.definition.name

    @property
    def repr_def(self):
        return self.definition.repr_def

    def is_superset_of(self, other: 'State') -> bool:
        return other.is_subset_of(self)

    def is_subset_of(self, other: 'State') -> bool:
        if self.definition == other.definition:
            return len(self.specs) == len(other.specs) and \
                   all(x.is_subspec_of(y) for x, y in zip(self.specs, other.specs)) and \
                   all(x == y for x, y in zip(self._non_spec_props, other._non_spec_props))
        else:
            return False

    @property
    def is_structured(self) -> bool:
        return all(mol_spec.struct_index is not None for mol_spec in self.specs)

    @property
    def is_global(self) -> bool:
        return self.name == GLOBAL_STATE_NAME

    @property
    def is_elemental(self) -> bool:
        return self.is_global or all(spec.has_resolution(elemental_resolution) for spec, elemental_resolution
                                     in zip(self.specs, self._elemental_resolutions))

    def is_mutually_exclusive_with(self, state: 'State') -> bool:
        assert self.is_elemental
        assert state.is_elemental

        if self == state:
            return False

        return any(spec in state.specs for spec in self.specs)

    @property
    def is_neutral(self) -> bool:
        return len(self.neutral_states) == 1 and self == self.neutral_states[0]

    @property
    def neutral_states(self) -> List['State']:
        return self.definition.neutral_states_from_vars(self.vars)

    @property
    def specs(self) -> List[Spec]:
        return [x for x in self.vars.values() if isinstance(x, Spec)]

    def update_specs(self, updates: Dict[Spec, Spec]):
        new_vars = deepcopy(self.vars)

        for k, v in self.vars.items():
            try:
                new_vars[k] = updates[v]
            except KeyError:
                pass

        self.vars = new_vars

    def to_non_structured(self) -> 'State':
        non_struct_vars = deepcopy(self.vars)
        for k, v in non_struct_vars.items():
            if isinstance(v, Spec):
                non_struct_vars[k] = v.to_non_struct_spec()

        return State(self.definition, non_struct_vars)

    def to_structured_from_spec(self, spec: Spec) -> 'State':
        assert spec.struct_index is not None
        struct_vars = deepcopy(self.vars)
        for k, v in struct_vars.items():
            if isinstance(v, Spec):
                if v.struct_index is not None:
                    continue
                elif spec.to_non_struct_spec().to_component_spec() == v.to_component_spec():
                    v.struct_index = spec.struct_index
                    break

        return State(self.definition, struct_vars)

    def to_structured_from_state(self, state: 'State') -> 'State':
        other_vars = deepcopy(state.vars)
        struct_vars = deepcopy(self.vars)

        for k, v in struct_vars.items():
            if isinstance(v, Spec):
                if v.struct_index is not None:
                    continue

                for kp, vp in other_vars.items():
                    if isinstance(vp, Spec) and vp.struct_index is not None and vp.to_non_struct_spec().is_superspec_of(v):
                        v.struct_index = vp.struct_index
                        other_vars[kp] = vp.to_non_struct_spec()
                        break

        return State(self.definition, struct_vars)

    @property
    def components(self) -> List[Spec]:
        return [x.to_component_spec() for x in self.specs]

    @property
    def _elemental_resolutions(self) -> List[LocusResolution]:
        return [self.definition.vars_def[var][1] for var, value in self.vars.items()
                if isinstance(value, Spec)]

    @property
    def _non_spec_props(self):
        return [x for x in self.vars.values() if not isinstance(x, Spec)]


class FullyNeutralState(State):
    def __init__(self) -> None:
        pass

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return 'fully-neutral-state'

    def __repr__(self) -> str:
        return str(self)

    def __lt__(self, other: 'State'):
        return True

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented
        return isinstance(other, FullyNeutralState)

    def to_non_structured(self):
        return self

    @property
    def is_structured(self) -> bool:
        return False

    def is_subset_of(self, other: 'State') -> bool:
        raise AssertionError

    @property
    def neutral_states(self) -> List['State']:
        raise AssertionError

    @property
    def components(self) -> List[Spec]:
        return []

    @property
    def is_elemental(self) -> bool:
        raise AssertionError

    @property
    def is_neutral(self) -> bool:
        raise AssertionError

    def is_superset_of(self, other: 'State') -> bool:
        raise AssertionError


def matching_state_def(repr: str) -> Optional[StateDef]:
    return next((x for x in STATE_DEFS if x.matches_repr(repr)), None)


def state_from_str(repr: str) -> State:
    repr = repr.strip()

    if repr == FULLY_NEUTRAL_STATE:
        return FullyNeutralState()

    state_def = matching_state_def(repr)

    if not state_def:
        raise SyntaxError('Could not match State {} with definition'.format(repr))

    variables = state_def.vars_from_repr(repr)
    state_def.validate_vars(variables)

    return State(state_def, variables)
