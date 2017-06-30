"""Module containing the abstract class State, classes InteractionState, EmptyBindingState, SelfInteractionState,
ModificationState, GlobalState, FullyNeutralState. Also contains constructor functions state_from_str,
state_modifier_from_str. The function initialize_state_modifiers allows runtime addition of state modifiers, and
is initially run when importing the module."""

import re
from enum import Enum
from typing import List, Dict
from abc import ABCMeta, abstractmethod
import logging
from copy import deepcopy

from rxncon.core.spec import Spec, LocusResolution, locus_from_str, spec_from_str

LOGGER = logging.getLogger(__name__)

SPEC_REGEX = r'([A-Za-z][A-Za-z0-9]*(?:@[\d]+)*(?:_\[[\w\/\(\)]+\])*)'
STR_REGEX = r'([\w]+)'
LOCUS_REGEX = r'([\w\/\(\)]+)'

FULLY_NEUTRAL_STATE_REGEX = r'^0$'
INTERACTION_STATE_REGEX = r'^{}--{}$'.format(SPEC_REGEX, SPEC_REGEX)
EMPTY_BINDING_STATE_REGEX = r'^{}--0$'.format(SPEC_REGEX)
SELF_INTERACTION_STATE_REGEX = r'^{}--\[{}\]$'.format(SPEC_REGEX, LOCUS_REGEX)
GLOBAL_STATE_REGEX = r'^\[{}\]$'.format(STR_REGEX)
MODIFICATION_STATE_REGEX = r'^' + SPEC_REGEX + r'-{' + STR_REGEX + r'}$'

StateModifier = None

DEFAULT_STATE_MODIFIERS = {
    'neutral': '0',
    'phosphor': 'p',
    'ubiquitin': 'ub',
    'GTP': 'GTP',
    'truncated': 'truncated'
}


def initialize_state_modifiers(additional_modifiers: Dict[str, str]=None) -> None:
    """Initializes the StateModifier name to be an Enum. Without passing in additional_modifiers,
    this function uses just the key-value pairs in DEFAULT_STATE_MODIFIERS. If additional pairs
    are passed in, for example through the dedicated Excel sheet, they are appended (overwriting
    possibly, the defaults.)"""
    if not additional_modifiers:
        additional_modifiers = {}

    global StateModifier
    modifiers = {
        **{k.lower(): v.lower() for k, v in DEFAULT_STATE_MODIFIERS.items()},
        **{k.lower(): v.lower() for k, v in additional_modifiers.items()}
    }

    StateModifier = Enum('StateModifier', modifiers)  # type: ignore


initialize_state_modifiers()


def state_modifier_from_str(modifier_str: str) -> StateModifier:
    try:
        return StateModifier(modifier_str.lower())
    except ValueError:
        valid_modifiers = [modifier.value for modifier in StateModifier.__members__.values()]  # type: ignore
        raise ValueError('Invalid StateModifier {}, valid modifiers are {}'
                         .format(modifier_str, ', '.join(valid_modifiers)))


class State(metaclass=ABCMeta):
    """Abstract Base Class for all States."""
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented
        else:
            raise NotImplementedError

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented
        else:
            raise NotImplementedError

    def clone(self) -> 'State':
        return deepcopy(self)

    @property
    @abstractmethod
    def is_elemental(self) -> bool:
        return False

    @abstractmethod
    def is_superset_of(self, other: 'State') -> bool:
        return False

    @abstractmethod
    def is_subset_of(self, other: 'State') -> bool:
        return False

    @property
    def is_structured(self) -> bool:
        return all(spec.is_structured for spec in self.specs)

    @property
    def is_global(self) -> bool:
        return False

    @abstractmethod
    def is_mutually_exclusive_with(self, state: 'State') -> bool:
        return False

    @property
    @abstractmethod
    def is_neutral(self) -> bool:
        return False

    @property
    @abstractmethod
    def is_homodimer(self) -> bool:
        return False

    @property
    @abstractmethod
    def neutral_states(self) -> List['State']:
        return []

    @property
    @abstractmethod
    def specs(self) -> List[Spec]:
        return []

    @abstractmethod
    def update_specs(self, updates: Dict[Spec, Spec]) -> None:
        pass

    @abstractmethod
    def to_non_structured(self) -> 'State':
        pass

    @abstractmethod
    def to_structured_from_spec(self, spec: Spec) -> 'State':
        pass

    def to_structured_from_state(self, state: 'State') -> 'State':
        structured = self.clone()
        for spec in state.specs:
            structured = structured.to_structured_from_spec(spec)

        return structured

    @property
    @abstractmethod
    def components(self) -> List[Spec]:
        pass


class InteractionState(State):
    """A State A_[x]--B_[y]."""
    def __init__(self, first: Spec, second: Spec) -> None:
        self.first, self.second = sorted([first, second])  # type: Spec, Spec
        self._validate()

    @property
    def name(self) -> str:
        return '{}--{}'.format(str(self.first), str(self.second))

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.name

    def __hash__(self) -> int:
        return hash(str(self))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented
        else:
            return isinstance(other, InteractionState) and self.first == other.first and self.second == other.second

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented

        return str(self) < str(other)

    def _validate(self) -> None:
        if self.first.resolution > LocusResolution.domain or self.second.resolution > LocusResolution.domain:
            raise SyntaxError(
                'Resolution for InteractionState too high {} {}'.format(str(self.first), str(self.second)))

        if self.first.is_structured and self.second.is_structured:
            assert not self.first.struct_index == self.second.struct_index, \
                'Interaction state {} has identical struct indices!'.format(self)

    @property
    def specs(self) -> List[Spec]:
        return [self.first.clone(), self.second.clone()]

    @property
    def is_global(self) -> bool:
        return False

    @property
    def is_elemental(self) -> bool:
        return self.first.has_resolution(LocusResolution.domain) and self.second.has_resolution(LocusResolution.domain)

    def is_mutually_exclusive_with(self, state: 'State') -> bool:
        if self == state:
            return False
        elif isinstance(state, InteractionState) or isinstance(state, SelfInteractionState):
            return self.first in (state.first, state.second) or self.second in (state.first, state.second)
        elif isinstance(state, EmptyBindingState):
            return self.first == state.spec or self.second == state.spec
        else:
            return False

    @property
    def is_neutral(self) -> bool:
        return False

    @property
    def is_homodimer(self) -> bool:
        return self.first.to_component_spec().to_non_struct_spec() == self.second.to_component_spec().to_non_struct_spec()

    @property
    def neutral_states(self) -> List['State']:
        if self.is_homodimer:
            return [EmptyBindingState(self.first)]
        else:
            return [EmptyBindingState(self.first), EmptyBindingState(self.second)]

    def update_specs(self, updates: Dict[Spec, Spec]) -> None:
        try:
            self.first = updates[self.first]
        except KeyError:
            pass
        try:
            self.second = updates[self.second]
        except KeyError:
            pass

        self.first, self.second = sorted([self.first, self.second])
        self._validate()

    def to_non_structured(self) -> 'State':
        return InteractionState(self.first.to_non_struct_spec(), self.second.to_non_struct_spec())

    @property
    def components(self) -> List[Spec]:
        return [self.first.to_component_spec(), self.second.to_component_spec()]

    def is_subset_of(self, other: 'State') -> bool:
        return isinstance(other, InteractionState) and ((self.first.is_subspec_of(other.first) and
                                                         self.second.is_subspec_of(other.second)) or (
                                                        self.first.is_subspec_of(other.second) and
                                                        self.second.is_subspec_of(other.first)))

    def is_superset_of(self, other: 'State') -> bool:
        return self == other or other.is_subset_of(self)

    def to_structured_from_spec(self, spec: Spec) -> 'State':
        if not spec.is_structured:
            return self

        if not self.first.is_structured and self.first.to_component_spec() == spec.to_non_struct_spec().to_component_spec():
            return InteractionState(self.first.with_struct_from_spec(spec), self.second)
        elif not self.second.is_structured and self.second.to_component_spec() == spec.to_non_struct_spec().to_component_spec():
            return InteractionState(self.first, self.second.with_struct_from_spec(spec))
        else:
            return self

    def to_structured_from_state(self, state: 'State') -> 'State':
        if self.is_homodimer and state.is_homodimer and state.is_structured:
            assert isinstance(state, InteractionState)
            return InteractionState(self.first.with_struct_from_spec(state.first),
                                    self.second.with_struct_from_spec(state.second))
        else:
            return super().to_structured_from_state(state)


class EmptyBindingState(State):
    """A State A_[x]--0."""
    def __init__(self, spec: Spec) -> None:
        if spec.resolution > LocusResolution.domain:
            raise SyntaxError('Resolution for EmptyBindingState too high {}'.format(str(spec)))
        self.spec = spec

    @property
    def name(self) -> str:
        return '{}--0'.format(str(self.spec))

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.name

    def __hash__(self) -> int:
        return hash(str(self))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented
        else:
            return isinstance(other, EmptyBindingState) and self.spec == other.spec

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented

        return str(self) < str(other)

    @property
    def specs(self) -> List[Spec]:
        return [self.spec.clone()]

    @property
    def is_global(self) -> bool:
        return False

    @property
    def is_elemental(self) -> bool:
        return self.spec.has_resolution(LocusResolution.domain)

    def is_mutually_exclusive_with(self, state: 'State') -> bool:
        if self == state:
            return False
        elif isinstance(state, InteractionState) or isinstance(state, SelfInteractionState):
            return self.spec == state.first or self.spec == state.second
        else:
            return False

    @property
    def is_neutral(self) -> bool:
        return True

    @property
    def is_homodimer(self) -> bool:
        return False

    @property
    def neutral_states(self) -> List['State']:
        return [self]

    def update_specs(self, updates: Dict[Spec, Spec]) -> None:
        try:
            self.spec = updates[self.spec]
        except KeyError:
            pass

    def to_non_structured(self) -> 'State':
        return EmptyBindingState(self.spec.to_non_struct_spec())

    @property
    def components(self) -> List[Spec]:
        return [self.spec.to_component_spec()]

    def is_subset_of(self, other: 'State') -> bool:
        return isinstance(other, EmptyBindingState) and self.spec.is_subspec_of(other.spec)

    def is_superset_of(self, other: 'State') -> bool:
        return self == other or other.is_subset_of(self)

    def to_structured_from_spec(self, spec: Spec) -> 'State':
        if spec.is_structured and self.spec.to_non_struct_spec().to_component_spec() == spec.to_non_struct_spec().to_component_spec():
            return EmptyBindingState(self.spec.with_struct_from_spec(spec))
        else:
            return self


class SelfInteractionState(State):
    """A State A_[x]--[y]."""
    def __init__(self, first: Spec, second: Spec) -> None:
        self.first, self.second = sorted([first, second])  # type: Spec, Spec
        self._validate()

    @property
    def name(self) -> str:
        return '{}--[{}]'.format(str(self.first), str(self.second.locus))

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.name

    def __hash__(self) -> int:
        return hash(str(self))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented
        else:
            return isinstance(other, SelfInteractionState) and self.first == other.first and self.second == other.second

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented

        return str(self) < str(other)

    def _validate(self) -> None:
        assert self.first.to_component_spec() == self.second.to_component_spec()
        if self.first.resolution > LocusResolution.domain or self.second.resolution > LocusResolution.domain:
            raise SyntaxError(
                'Resolution for SelfInteractionState too high {} {}'.format(str(self.first), str(self.second)))

    @property
    def specs(self) -> List[Spec]:
        return [self.first.clone(), self.second.clone()]

    @property
    def is_global(self) -> bool:
        return False

    @property
    def is_elemental(self) -> bool:
        return self.first.has_resolution(LocusResolution.domain) and \
               self.second.has_resolution(LocusResolution.domain)

    def is_mutually_exclusive_with(self, state: 'State') -> bool:
        if self == state:
            return False
        elif isinstance(state, InteractionState) or isinstance(state, SelfInteractionState):
            return self.first == state.first or self.second == state.second
        elif isinstance(state, EmptyBindingState):
            return self.first == state.spec or self.second == state.spec
        else:
            return False

    @property
    def is_neutral(self) -> bool:
        return False

    @property
    def is_homodimer(self) -> bool:
        return False

    @property
    def neutral_states(self) -> List['State']:
        return [EmptyBindingState(self.first), EmptyBindingState(self.second)]

    def update_specs(self, updates: Dict[Spec, Spec]) -> None:
        try:
            self.first = updates[self.first]
        except KeyError:
            pass
        try:
            self.second = updates[self.second]
        except KeyError:
            pass

        self.first, self.second = sorted([self.first, self.second])
        self._validate()

    def to_non_structured(self) -> 'State':
        return SelfInteractionState(self.first.to_non_struct_spec(), self.second.to_non_struct_spec())

    @property
    def components(self) -> List[Spec]:
        return [self.first.to_component_spec()]

    def is_subset_of(self, other: 'State') -> bool:
        return isinstance(other, SelfInteractionState) and self.first.is_subspec_of(other.first) and \
               self.second.is_subspec_of(other.second)

    def is_superset_of(self, other: 'State') -> bool:
        return self == other or other.is_subset_of(self)

    def to_structured_from_spec(self, spec: Spec) -> 'State':
        if self.first.to_non_struct_spec().to_component_spec() == spec.to_non_struct_spec().to_component_spec():
            assert self.second.to_non_struct_spec().to_component_spec() == spec.to_non_struct_spec().to_component_spec()
            return SelfInteractionState(self.first.with_struct_from_spec(spec), self.second.with_struct_from_spec(spec))
        else:
            return self


class ModificationState(State):
    """A State A_[(r)]-{p}."""
    def __init__(self, spec: Spec, modifier: StateModifier) -> None:
        self.spec, self.modifier = spec, modifier

    @property
    def name(self) -> str:
        return '{}-{{{}}}'.format(str(self.spec), self.modifier.value)

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.name

    def __hash__(self) -> int:
        return hash(str(self))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented
        else:
            return isinstance(other, ModificationState) and self.spec == other.spec and self.modifier == other.modifier

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented

        return str(self) < str(other)

    @property
    def specs(self) -> List[Spec]:
        return [self.spec.clone()]

    @property
    def is_global(self) -> bool:
        return False

    @property
    def is_elemental(self) -> bool:
        return self.spec.has_resolution(LocusResolution.residue)

    def is_mutually_exclusive_with(self, state: 'State') -> bool:
        if self == state:
            return False
        elif isinstance(state, ModificationState):
            return self.spec == state.spec
        else:
            return False

    @property
    def is_neutral(self) -> bool:
        return self.modifier == StateModifier.neutral  # type: ignore

    @property
    def is_homodimer(self) -> bool:
        return False

    @property
    def neutral_states(self) -> List['State']:
        return [ModificationState(self.spec, StateModifier.neutral)]  # type: ignore

    def update_specs(self, updates: Dict[Spec, Spec]) -> None:
        try:
            self.spec = updates[self.spec]
        except KeyError:
            pass

    def to_non_structured(self) -> 'State':
        return ModificationState(self.spec.to_non_struct_spec(), self.modifier)

    @property
    def components(self) -> List[Spec]:
        return [self.spec.to_component_spec()]

    def is_subset_of(self, other: 'State') -> bool:
        return isinstance(other, ModificationState) and self.spec.is_subspec_of(
            other.spec) and self.modifier == other.modifier

    def is_superset_of(self, other: 'State') -> bool:
        return self == other or other.is_subset_of(self)

    def to_structured_from_spec(self, spec: Spec) -> 'State':
        if self.spec.to_non_struct_spec().to_component_spec() == spec.to_non_struct_spec().to_component_spec() and spec.is_structured:
            return ModificationState(self.spec.with_struct_from_spec(spec), self.modifier)
        else:
            return self


class GlobalState(State):
    """A State [Turgor]."""
    def __init__(self, name: str) -> None:
        self.name = '[{}]'.format(name.strip('[]'))

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.name

    def __hash__(self) -> int:
        return hash(str(self))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented
        else:
            return isinstance(other, GlobalState) and self.name == other.name

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented

        return str(self) < str(other)

    @property
    def specs(self) -> List[Spec]:
        return []

    @property
    def is_global(self) -> bool:
        return True

    @property
    def is_elemental(self) -> bool:
        return True

    def is_mutually_exclusive_with(self, state: 'State') -> bool:
        return False

    @property
    def is_neutral(self) -> bool:
        return False

    @property
    def is_homodimer(self) -> bool:
        return False

    @property
    def neutral_states(self) -> List['State']:
        return []

    def update_specs(self, updates: Dict[Spec, Spec]) -> None:
        pass

    def to_non_structured(self) -> 'State':
        return self

    @property
    def is_structured(self) -> bool:
        return True

    @property
    def components(self) -> List[Spec]:
        return []

    def is_subset_of(self, other: 'State') -> bool:
        return self == other

    def is_superset_of(self, other: 'State') -> bool:
        return self == other or other.is_subset_of(self)

    def to_structured_from_spec(self, spec: Spec) -> 'State':
        return self


class FullyNeutralState(State):
    """A State 0. This happens in synthesis reactions, and serves as a placeholder for the actual
    combination of neutral states it is going to be replaced with. This actual set of states is unknown
    until the entire system is known."""
    def __init__(self) -> None:
        self.name = 'FullyNeutralState'

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return 'fully-neutral-state'

    def __repr__(self) -> str:
        return str(self)

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented
        return True

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented
        return isinstance(other, FullyNeutralState)

    def to_non_structured(self) -> State:
        return self

    @property
    def is_structured(self) -> bool:
        return False

    def is_subset_of(self, other: 'State') -> bool:
        raise NotImplementedError

    @property
    def neutral_states(self) -> List['State']:
        raise NotImplementedError

    @property
    def components(self) -> List[Spec]:
        return []

    @property
    def is_elemental(self) -> bool:
        raise NotImplementedError

    @property
    def is_neutral(self) -> bool:
        raise NotImplementedError

    def is_superset_of(self, other: 'State') -> bool:
        raise NotImplementedError

    def to_structured_from_spec(self, spec: Spec) -> 'State':
        raise NotImplementedError

    @property
    def is_global(self) -> bool:
        raise NotImplementedError

    def is_mutually_exclusive_with(self, state: 'State') -> bool:
        raise NotImplementedError

    def to_structured_from_state(self, state: 'State') -> 'State':
        raise NotImplementedError

    @property
    def specs(self) -> List[Spec]:
        raise NotImplementedError

    def update_specs(self, updates: Dict[Spec, Spec]) -> None:
        raise NotImplementedError

    @property
    def is_homodimer(self) -> bool:
        raise NotImplementedError


def state_from_str(state_str: str) -> State:  # pylint: disable=too-many-return-statements
    state_str = state_str.strip()

    if re.match(GLOBAL_STATE_REGEX, state_str):
        return GlobalState(state_str)
    elif re.match(INTERACTION_STATE_REGEX, state_str):
        first_str, second_str = re.findall(INTERACTION_STATE_REGEX, state_str)[0]
        return InteractionState(spec_from_str(first_str), spec_from_str(second_str))
    elif re.match(EMPTY_BINDING_STATE_REGEX, state_str):
        spec_str = re.findall(EMPTY_BINDING_STATE_REGEX, state_str)[0]
        return EmptyBindingState(spec_from_str(spec_str))
    elif re.match(SELF_INTERACTION_STATE_REGEX, state_str):
        first_str, second_str = re.findall(SELF_INTERACTION_STATE_REGEX, state_str)[0]
        first_spec = spec_from_str(first_str)
        second_spec = spec_from_str(first_str).with_locus(locus_from_str(second_str))
        return SelfInteractionState(first_spec, second_spec)
    elif re.match(GLOBAL_STATE_REGEX, state_str):
        return GlobalState(state_str.strip('[]'))
    elif re.match(MODIFICATION_STATE_REGEX, state_str):
        spec_str, mod_str = re.findall(MODIFICATION_STATE_REGEX, state_str)[0]
        return ModificationState(spec_from_str(spec_str), state_modifier_from_str(mod_str))
    elif re.match(FULLY_NEUTRAL_STATE_REGEX, state_str):
        return FullyNeutralState()
    else:
        raise SyntaxError('Unable to parse State string {}'.format(state_str))
