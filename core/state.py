from abc import ABCMeta
from enum import Enum, unique

import core.component as com


@unique
class StateModifier(Enum):
    undefined = None
    phosphor  = 'P'
    ubiquitin = 'Ub'
    truncated = 'Truncated'


class State:
    __metaclass__ = ABCMeta

    def __str__(self) -> str:
        return self.full_name


class CovalentModificationState(State):
    def __init__(self, substrate: com.Component, modifier: StateModifier):
        self.substrate = substrate
        self.modifier = modifier
        self.full_name = '{0}-{{{1}}}'.format(self.substrate.full_name, self.modifier.value)

    def __eq__(self, other: State) -> bool:
        assert isinstance(other, State)

        if isinstance(other, CovalentModificationState):
            return self.substrate == other.substrate and self.modifier == other.modifier
        else:
            return False


class InteractionState(State):
    def __init__(self, first_component: com.Component, second_component: com.Component):
        self.first_component = first_component
        self.second_component = second_component
        self.full_name = '{0}--{1}'.format(self.first_component.full_name, self.second_component.full_name)

    def __eq__(self, other: State) -> bool:
        assert isinstance(other, State)

        if isinstance(other, InteractionState):
            return self.first_component == other.first_component and self.second_component == other.second_component
        else:
            return False


class SynthesisDegradationState(State):
    def __init__(self, component: com.Component):
        self.component = component
        self.full_name = '{}'.format(self.component.full_name)

    def __eq__(self, other: State) -> bool:
        assert isinstance(other, State)

        if isinstance(other, SynthesisDegradationState):
            return self.component == other.component
        else:
            return False


class TranslocationState(State):
    def __init__(self, substrate: com.Component, compartment: StateModifier):
        self.substrate = substrate
        self.compartment = compartment
        self.full_name = '{0}-{{{1}}}'.format(self.substrate.full_name, self.compartment.value)

    def __eq__(self, other: State) -> bool:
        assert isinstance(other, State)

        if isinstance(other, TranslocationState):
            return self.substrate == other.substrate and self.compartment == other.compartment
        else:
            return False


class OutputState(State):
    def __init__(self, full_name: str):
        self.full_name = full_name
