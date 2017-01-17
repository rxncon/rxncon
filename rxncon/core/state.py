import re
from enum import unique, Enum
from typing import List, Dict
from abc import ABC, abstractproperty, abstractmethod
import logging

from rxncon.core.spec import Spec, Locus, LocusResolution, locus_from_str, spec_from_str
from rxncon.util.utils import current_function_name


logger = logging.getLogger(__name__)


SPEC_REGEX          = '([A-Za-z][A-Za-z0-9]*(?:@[\d]+)*(?:_\[[\w\/\(\)]+\])*)'
STR_REGEX           = '([\w]+)'
LOCUS_REGEX         = '([\w\/\(\)]+)'


FULLY_NEUTRAL_STATE_REGEX    = '^0$'
INTERACTION_STATE_REGEX      = '^{}--{}$'.format(SPEC_REGEX, SPEC_REGEX)
EMPTY_BINDING_STATE_REGEX    = '^{}--0$'.format(SPEC_REGEX)
SELF_INTERACTION_STATE_REGEX = '^{}--\[{}\]$'.format(SPEC_REGEX, LOCUS_REGEX)
GLOBAL_STATE_REGEX           = '^\[{}\]$'.format(STR_REGEX)
MODIFICATION_STATE_REGEX     = '^' + SPEC_REGEX + '-{' + STR_REGEX + '}$'


@unique
class StateModifier(Enum):
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


class State(ABC):
    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return self.name

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented
        return str(self) < str(other)

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

    @property
    def is_elemental(self) -> bool:
        return self.is_global or all(spec.has_resolution(elemental_resolution) for spec, elemental_resolution
                                     in zip(self.specs, self._elemental_resolutions))

    @abstractmethod
    def is_mutually_exclusive_with(self, state: 'State') -> bool:
        return False

    @abstractproperty
    def is_neutral(self) -> bool:
        return False

    @abstractproperty
    def is_homodimer(self):
        return False

    @abstractproperty
    def neutral_states(self) -> List['State']:
        return []

    @abstractproperty
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

    @abstractmethod
    def to_structured_from_state(self, state: 'State') -> 'State':
        pass

    @abstractmethod
    def components(self) -> List[Spec]:
        pass


class InteractionState(State):
    def __init__(self, first: Spec, second: Spec):
        if first.resolution > LocusResolution.domain or second.resolution > LocusResolution.domain:
            raise SyntaxError('Resolution for InteractionState too high {} {}'.format(str(first), str(second)))
        self.first, self.second = sorted([first, second])  # type: Spec, Spec
        self.repr_def = '$x--$y'

    @property
    def name(self):
        return '{}--{}'.format(str(self.first), str(self.second))

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.name

    def __hash__(self) -> int:
        return hash(str(self))

    def __eq__(self, other: object):
        if not isinstance(other, State):
            return NotImplemented
        else:
            return isinstance(other, InteractionState) and self.first == other.first and self.second == other.second

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented

        return str(self) < str(other)

    def __getitem__(self, item):
        if item == '$x':
            return self.first
        elif item == '$y':
            return self.second
        else:
            raise KeyError

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
            return self.first in (state.first, state.second) or self.second in (state.first, state.second)
        elif isinstance(state, EmptyBindingState):
            return self.first == state.spec or self.second == state.spec
        else:
            return False

    @property
    def is_neutral(self) -> bool:
        return False

    @property
    def is_homodimer(self):
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

    def to_non_structured(self) -> 'State':
        return InteractionState(self.first.to_non_struct_spec(), self.second.to_non_struct_spec())

    @property
    def is_structured(self) -> bool:
        return self.first.is_structured and self.second.is_structured

    @property
    def components(self) -> List[Spec]:
        return [self.first.to_component_spec(), self.second.to_component_spec()]

    def is_subset_of(self, other: 'State') -> bool:
        return isinstance(other, InteractionState) and self.first.is_subspec_of(other.first) and \
            self.second.is_subspec_of(other.second)

    def is_superset_of(self, other: 'State') -> bool:
        return self == other or other.is_subset_of(self)

    def to_structured_from_spec(self, spec: Spec) -> 'State':
        if not spec.is_structured:
            return self

        if self.is_homodimer:
            if not self.first.is_structured:
                return InteractionState(self.first.with_struct_from_spec(spec), self.second)
            elif not self.second.is_structured:
                return InteractionState(self.first, self.second.with_struct_from_spec(spec))
            else:
                logger.info('{} : Not structuring homodimer {} with spec {}, since pre-structured.'
                            .format(current_function_name(), str(self), str(spec)))
                return self
        else:
            if self.first.to_non_struct_spec().to_component_spec() == spec.to_non_struct_spec().to_component_spec():
                return InteractionState(self.first.with_struct_from_spec(spec), self.second)
            elif self.second.to_non_struct_spec().to_component_spec() == spec.to_non_struct_spec().to_component_spec():
                return InteractionState(self.first, self.second.with_struct_from_spec(spec))
            else:
                return self

    def to_structured_from_state(self, state: 'State') -> 'State':
        if self.is_homodimer and state.is_homodimer and state.is_structured:
            return InteractionState(self.first.with_struct_from_spec(state.first), self.second.with_struct_from_spec(state.second))
        else:
            for spec in state.specs:
                self = self.to_structured_from_spec(spec)

            return self


class EmptyBindingState(State):
    def __init__(self, spec: Spec):
        if spec.resolution > LocusResolution.domain:
            raise SyntaxError('Resolution for EmptyBindingState too high {}'.format(str(spec)))
        self.spec = spec

    @property
    def name(self):
        return '{}--0'.format(str(self.spec))

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.name

    def __hash__(self) -> int:
        return super().__hash__()

    def __eq__(self, other: object):
        if not isinstance(other, State):
            return NotImplemented
        else:
            return isinstance(other, EmptyBindingState) and self.spec == other.spec

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented

        return str(self) < str(other)

    def __getitem__(self, item):
        if item == '$x':
            return self.spec
        else:
            raise KeyError

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
    def is_homodimer(self):
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
    def is_structured(self) -> bool:
        return self.spec.is_structured

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

    def to_structured_from_state(self, state: 'State') -> 'State':
        for spec in state.specs:
            self = self.to_structured_from_spec(spec)

        return self


class SelfInteractionState(State):
    def __init__(self, first: Spec, second: Spec):
        assert first.to_component_spec() == second.to_component_spec()
        if first.resolution > LocusResolution.domain or second.resolution > LocusResolution.domain:
            raise SyntaxError('Resolution for SelfInteractionState too high {} {}'.format(str(first), str(second)))
        self.first, self.second = sorted([first, second])  # type: Spec, Spec

    @property
    def name(self):
        return '{}--[{}]'.format(str(self.first), str(self.second.locus))

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.name

    def __hash__(self) -> int:
        return super().__hash__()

    def __eq__(self, other: object):
        if not isinstance(other, State):
            return NotImplemented
        else:
            return isinstance(other, SelfInteractionState) and self.first == other.first and self.second == other.second

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented

        return str(self) < str(other)

    def __getitem__(self, item):
        if item == '$x':
            return self.first
        elif item == '$y':
            return self.second
        else:
            raise KeyError

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
    def is_homodimer(self):
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

    def to_non_structured(self) -> 'State':
        return SelfInteractionState(self.first.to_non_struct_spec(), self.second.to_non_struct_spec())

    @property
    def is_structured(self) -> bool:
        return self.first.is_structured and self.second.is_structured

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

    def to_structured_from_state(self, state: 'State') -> 'State':
        for spec in state.specs:
            self = self.to_structured_from_spec(spec)

        return self


class ModificationState(State):
    def __init__(self, spec: Spec, modifier: StateModifier):
        self.spec, self.modifier = spec, modifier

    @property
    def name(self):
        return '{}-{{{}}}'.format(str(self.spec), self.modifier.value)

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.name

    def __hash__(self) -> int:
        return super().__hash__()

    def __eq__(self, other: object):
        if not isinstance(other, State):
            return NotImplemented
        else:
            return isinstance(other, ModificationState) and self.spec == other.spec and self.modifier == other.modifier

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented

        return str(self) < str(other)

    def __getitem__(self, item):
        if item == '$x':
            return self.spec
        elif item == '$y':
            return self.modifier
        else:
            raise KeyError

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
        return self.modifier == StateModifier.neutral

    @property
    def is_homodimer(self):
        return False

    @property
    def neutral_states(self) -> List['State']:
        return [ModificationState(self.spec, StateModifier.neutral)]

    def update_specs(self, updates: Dict[Spec, Spec]) -> None:
        try:
            self.spec = updates[self.spec]
        except KeyError:
            pass

    def to_non_structured(self) -> 'State':
        return ModificationState(self.spec.to_non_struct_spec(), self.modifier)

    @property
    def is_structured(self) -> bool:
        return self.spec.is_structured

    @property
    def components(self) -> List[Spec]:
        return [self.spec.to_component_spec()]

    def is_subset_of(self, other: 'State') -> bool:
        return isinstance(other, ModificationState) and self.spec.is_subspec_of(other.spec) and self.modifier == other.modifier

    def is_superset_of(self, other: 'State') -> bool:
        return self == other or other.is_subset_of(self)

    def to_structured_from_spec(self, spec: Spec) -> 'State':
        if self.spec.to_non_struct_spec().to_component_spec() == spec.to_non_struct_spec().to_component_spec() and spec.is_structured:
            return ModificationState(self.spec.with_struct_from_spec(spec), self.modifier)
        else:
            return self

    def to_structured_from_state(self, state: 'State') -> 'State':
        for spec in state.specs:
            self = self.to_structured_from_spec(spec)

        return self


class GlobalState(State):
    def __init__(self, name: str):
        self.name = '[{}]'.format(name.strip('[]'))

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.name

    def __hash__(self) -> int:
        return super().__hash__()

    def __eq__(self, other: object):
        if not isinstance(other, State):
            return NotImplemented
        else:
            return isinstance(other, GlobalState) and self.name == other.name

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, State):
            return NotImplemented

        return str(self) < str(other)

    def __getitem__(self, item):
        if item == '$x':
            return self.name.strip('[]')
        else:
            raise KeyError

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
    def is_homodimer(self):
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

    def to_structured_from_state(self, state: 'State') -> 'State':
        return self


class FullyNeutralState(State):
    def __init__(self) -> None:
        pass

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

    def to_structured_from_spec(self, spec: Spec) -> 'State':
        raise AssertionError

    @property
    def is_global(self) -> bool:
        raise AssertionError

    def is_mutually_exclusive_with(self, state: 'State') -> bool:
        raise AssertionError

    def to_structured_from_state(self, state: 'State') -> 'State':
        raise AssertionError

    @property
    def specs(self) -> List[Spec]:
        raise AssertionError

    def update_specs(self, updates: Dict[Spec, Spec]) -> None:
        raise AssertionError

    @property
    def is_homodimer(self):
        raise AssertionError


def state_from_str(state_str: str) -> State:
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
