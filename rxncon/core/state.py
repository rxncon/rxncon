from abc import ABCMeta, abstractmethod
from enum import Enum, unique
import typecheck as tc

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


class CovalentModificationState(State):
    @tc.typecheck
    def __init__(self, substrate: com.Component, modifier: StateModifier):
        self.substrate = substrate
        self.modifier = modifier

    @tc.typecheck
    def __eq__(self, other: State) -> bool:
        return isinstance(other, CovalentModificationState) and self.substrate == other.substrate and self.modifier == other.modifier

    def __hash__(self):
        return hash('*cov-mod-state-{}-{}*'.format(self.substrate, self.modifier))

    def __str__(self) -> str:
        return sfr.string_from_covalent_modification_state(self)


class InterProteinInteractionState(State):
    @tc.typecheck
    def __init__(self, first_component: com.Component, second_component: com.Component):
        self.first_component = first_component
        self.second_component = second_component

    @tc.typecheck
    def __eq__(self, other: State) -> bool:
        return isinstance(other, InterProteinInteractionState) and self.first_component == other.first_component and \
            self.second_component == other.second_component

    def __hash__(self):
        return hash('*ppi-state-{}-{}*'.format(self.first_component, self.second_component))

    def __str__(self) -> str:
        return sfr.string_from_inter_protein_interaction_state(self)


class IntraProteinInteractionState(State):
    @tc.typecheck
    def __init__(self, first_component: com.Component, second_component: com.Component):
        assert first_component.name == second_component.name
        self.first_component = first_component
        self.second_component = second_component

    @tc.typecheck
    def __eq__(self, other: State) -> bool:
        return isinstance(other, IntraProteinInteractionState) and self.first_component == other.first_component and \
            self.second_component == other.second_component

    def __hash__(self):
        return hash('*ipi-state-{}-{}*'.format(self.first_component, self.second_component))

    def __str__(self) -> str:
        return sfr.string_from_intra_protein_interaction_state(self)


class SynthesisDegradationState(State):
    @tc.typecheck
    def __init__(self, component: com.Component):
        self.component = component

    @tc.typecheck
    def __eq__(self, other: State) -> bool:
        return isinstance(other, SynthesisDegradationState) and self.component == other.component

    def __hash__(self):
        return hash('*synth-deg-state-{}*'.format(self.component))

    def __str__(self) -> str:
        return sfr.string_from_synthesis_degradation_state(self)


class TranslocationState(State):
    @tc.typecheck
    def __init__(self, substrate: com.Component, compartment: StateModifier):
        self.substrate = substrate
        self.compartment = compartment

    @tc.typecheck
    def __eq__(self, other: State) -> bool:
        return isinstance(other, TranslocationState) and self.substrate == other.substrate and \
            self.compartment == other.compartment

    def __hash__(self):
        return hash('*transloc-state-{}-{}*'.format(self.substrate, self.compartment))

    def __str__(self) -> str:
        return sfr.string_from_translocation_state(self)


class InputState(State):
    @tc.typecheck
    def __init__(self, name: str):
        self.name = name

    def __hash__(self):
        return hash('*input-state-{}*'.format(self.name))
