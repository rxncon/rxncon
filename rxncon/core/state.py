from abc import ABCMeta, abstractmethod, abstractproperty
from enum import Enum, unique
import typecheck as tc
from typing import List

import rxncon.core.specification as com
import rxncon.syntax.string_from_rxncon as sfr


@unique
class StateModifier(Enum):
    unmodified = 'u'
    phosphor   = 'p'
    ubiquitin  = 'ub'
    truncated  = 'truncated'


class State(metaclass=ABCMeta):

    def __repr__(self) -> str:
        return str(self)

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def is_superspecification_of(self, other) -> bool:
        pass

    @abstractmethod
    def is_subspecification_of(self, other) -> bool:
        pass

    @abstractproperty
    def components(self) -> List[com.Specification]:
        pass

class ComponentState(State):

    @tc.typecheck
    def __init__(self, component: com.Specification):
        self.component = component
        self._validation()

    def _validation(self):
        assert self.component.name is not None
        assert self.component.domain is None
        assert self.component.subdomain is None
        assert self.component.residue is None

    @tc.typecheck
    def __eq__(self, other):
        return isinstance(other, ComponentState) and self.component == other.component

    def __hash__(self) -> int:
        return hash("*comp-state-{}".format(self.component))

    def __repr__(self):
        return str(self)

    def __str__(self):
        return sfr.string_from_component_state(self)

    @tc.typecheck
    def is_superspecification_of(self, other: State) -> bool:
        if isinstance(other, ComponentState) and self.component.is_superspecification_of(other.component):
            return True
        elif isinstance(other, InterProteinInteractionState) or isinstance(other, IntraProteinInteractionState):
            return self.component.is_superspecification_of(other.first_component) or self.component.is_superspecification_of(other.second_component)
        elif isinstance(other, CovalentModificationState):
            return self.component.is_superspecification_of(other.substrate)
        elif isinstance(other, TranslocationState):
            return self.component.is_superspecification_of(other.substrate)
        else:
            raise NotImplementedError

    def is_subspecification_of(self, other: State) -> bool:
        return isinstance(other, ComponentState) and self.component.is_subspecification_of(other.component)

    @property
    @tc.typecheck
    def components(self) -> List[com.Specification]:
        return [self.component]

class CovalentModificationState(State):
    @tc.typecheck
    def __init__(self, substrate: com.Specification, modifier: StateModifier):
        self.substrate = substrate
        self.modifier = modifier

    @tc.typecheck
    def __eq__(self, other: State) -> bool:
        return isinstance(other, CovalentModificationState) and self.substrate == other.substrate and self.modifier == other.modifier

    def __hash__(self) -> int:
        return hash('*cov-mod-state-{}-{}*'.format(self.substrate, self.modifier))

    def __str__(self) -> str:
        return sfr.string_from_covalent_modification_state(self)

    @tc.typecheck
    def is_superspecification_of(self, other: State) -> bool:
        return isinstance(other, CovalentModificationState) and self.modifier == other.modifier and \
            self.substrate.is_superspecification_of(other.substrate)

    @tc.typecheck
    def is_subspecification_of(self, other: State) -> bool:
        return isinstance(other, CovalentModificationState) and self.modifier == other.modifier and \
            self.substrate.is_subspecification_of(other.substrate)

    @property
    @tc.typecheck
    def components(self) -> List[com.Specification]:
        return [self.substrate]


class InterProteinInteractionState(State):
    @tc.typecheck
    def __init__(self, first_component: com.Specification, second_component: com.Specification):
        self.first_component = first_component
        self.second_component = second_component

    @tc.typecheck
    def __eq__(self, other: State) -> bool:
        return isinstance(other, InterProteinInteractionState) and self.first_component == other.first_component and \
            self.second_component == other.second_component

    def __hash__(self) -> int:
        return hash('*ppi-state-{}-{}*'.format(self.first_component, self.second_component))

    def __str__(self) -> str:
        return sfr.string_from_inter_protein_interaction_state(self)

    @tc.typecheck
    def is_superspecification_of(self, other: State) -> bool:
        return isinstance(other, InterProteinInteractionState) and self.first_component.is_superspecification_of(other.first_component) \
            and self.second_component.is_superspecification_of(other.second_component)

    @tc.typecheck
    def is_subspecification_of(self, other: State) -> bool:
        return isinstance(other, InterProteinInteractionState) and self.first_component.is_subspecification_of(other.first_component) \
            and self.second_component.is_subspecification_of(other.second_component)

    @property
    @tc.typecheck
    def components(self) -> List[com.Specification]:
        return [self.first_component, self.second_component]


class IntraProteinInteractionState(State):
    @tc.typecheck
    def __init__(self, first_component: com.Specification, second_component: com.Specification):
        assert first_component.name == second_component.name
        self.first_component = first_component
        self.second_component = second_component

    @tc.typecheck
    def __eq__(self, other: State) -> bool:
        return isinstance(other, IntraProteinInteractionState) and self.first_component == other.first_component and \
            self.second_component == other.second_component

    def __hash__(self) -> int:
        return hash('*ipi-state-{}-{}*'.format(self.first_component, self.second_component))

    def __str__(self) -> str:
        return sfr.string_from_intra_protein_interaction_state(self)

    @tc.typecheck
    def is_superspecification_of(self, other: State) -> bool:
        return isinstance(other, IntraProteinInteractionState) and self.first_component.is_superspecification_of(other.first_component) \
            and self.second_component.is_superspecification_of(other.second_component)

    @tc.typecheck
    def is_subspecification_of(self, other: State) -> bool:
        return isinstance(other, IntraProteinInteractionState) and self.first_component.is_subspecification_of(other.first_component) \
            and self.second_component.is_subspecification_of(other.second_component)

    @property
    @tc.typecheck
    def components(self) -> List[com.Specification]:
        return [self.first_component, self.second_component]


class TranslocationState(State):
    @tc.typecheck
    def __init__(self, substrate: com.Specification, compartment: str):
        self.substrate = substrate
        self.compartment = compartment

    @tc.typecheck
    def __eq__(self, other: State) -> bool:
        return isinstance(other, TranslocationState) and self.substrate == other.substrate and \
            self.compartment == other.compartment

    def __hash__(self) -> int:
        return hash('*transloc-state-{}-{}*'.format(self.substrate, self.compartment))

    def __str__(self) -> str:
        return sfr.string_from_translocation_state(self)

    @tc.typecheck
    def is_superspecification_of(self, other: State) -> bool:
        return isinstance(other, TranslocationState) and self.compartment == other.compartment and \
            self.substrate.is_superspecification_of(other.substrate)

    @tc.typecheck
    def is_subspecification_of(self, other: State) -> bool:
        return isinstance(other, TranslocationState) and self.compartment == other.compartment and \
            self.substrate.is_subspecification_of(other.substrate)

    @property
    @tc.typecheck
    def components(self) -> List[com.Specification]:
        return [self.substrate]


class InputState(State):
    @tc.typecheck
    def __init__(self, name: str):
        self.name = name

    def __hash__(self) -> int:
        return hash('*input-state-{}*'.format(self.name))

    def __str__(self) -> str:
        return sfr.string_from_input_state(self)

    @tc.typecheck
    def is_superspecification_of(self, other: State) -> bool:
        raise NotImplementedError

    @tc.typecheck
    def is_subspecification_of(self, other: State) -> bool:
        raise NotImplementedError

    @property
    @tc.typecheck
    def components(self) -> List[com.Specification]:
        return []
