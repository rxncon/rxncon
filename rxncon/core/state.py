from abc import ABCMeta, abstractmethod
from enum import Enum, unique

import rxncon.core.component as com
import rxncon.syntax.string_from_rxncon as sfr


@unique
class StateModifier(Enum):
    undefined = None
    phosphor  = 'p'
    ubiquitin = 'ub'
    truncated = 'truncated'


class State:
    __metaclass__ = ABCMeta

    def __repr__(self):
        return str(self)

    @abstractmethod
    def is_subset_of(self, other: 'State') -> bool:
        pass

    @abstractmethod
    def is_superset_of(self, other: 'State') -> bool:
        pass


class CovalentModificationState(State):
    def __init__(self, substrate: com.Component, modifier: StateModifier):
        self.substrate = substrate
        self.modifier = modifier

    def __eq__(self, other: State) -> bool:
        assert isinstance(other, State)

        if isinstance(other, CovalentModificationState):
            return self.substrate == other.substrate and self.modifier == other.modifier

        else:
            return False

    def __hash__(self):
        return hash('*cov-mod-state-{}-{}*'.format(self.substrate, self.modifier))

    def __str__(self) -> str:
        return sfr.string_from_covalent_modification_state(self)


# @todo rename to InterProteinInteractionState, create new IntraProteinInteractionState
class InterProteinInteractionState(State):
    def __init__(self, first_component: com.Component, second_component: com.Component):
        self.first_component = first_component
        self.second_component = second_component

    def __eq__(self, other: State) -> bool:
        assert isinstance(other, State)

        if isinstance(other, InterProteinInteractionState):
            return self.first_component == other.first_component and self.second_component == other.second_component

        else:
            return False

    def __hash__(self):
        return hash('*ppi-state-{}-{}*'.format(self.first_component, self.second_component))

    def __str__(self) -> str:
        return sfr.string_from_inter_protein_interaction_state(self)


class IntraProteinInteractionState(State):
    def __init__(self, first_component: com.Component, second_component: com.Component):
        assert first_component.name == second_component.name
        self.first_component = first_component
        self.second_component = second_component

    def __eq__(self, other: State) -> bool:
        assert isinstance(other, State)

        if isinstance(other, IntraProteinInteractionState):
            return self.first_component == other.first_component and self.second_component == other.second_component

        else:
            return False

    def __hash__(self):
        return hash('*ipi-state-{}-{}*'.format(self.first_component, self.second_component))

    def __str__(self) -> str:
        return sfr.string_from_intra_protein_interaction_state(self)


class SynthesisDegradationState(State):
    def __init__(self, component: com.Component):
        self.component = component

    def __eq__(self, other: State) -> bool:
        assert isinstance(other, State)

        if isinstance(other, SynthesisDegradationState):
            return self.component == other.component

        else:
            return False

    def __hash__(self):
        return hash('*synth-deg-state-{}*'.format(self.component))

    def __str__(self) -> str:
        return sfr.string_from_synthesis_degradation_state(self)


class TranslocationState(State):
    def __init__(self, substrate: com.Component, compartment: StateModifier):
        self.substrate = substrate
        self.compartment = compartment

    def __eq__(self, other: State) -> bool:
        assert isinstance(other, State)

        if isinstance(other, TranslocationState):
            return self.substrate == other.substrate and self.compartment == other.compartment

        else:
            return False

    def __hash__(self):
        return hash('*transloc-state-{}-{}*'.format(self.substrate, self.compartment))

    def __str__(self) -> str:
        return sfr.string_from_translocation_state(self)


class InputState(State):
    def __init__(self, name: str):
        self.name = name

    def __hash__(self):
        return hash('*input-state-{}*'.format(self.name))
