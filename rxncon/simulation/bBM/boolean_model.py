from rxncon.venntastic.sets import Set as VennSet, ValueSet, Intersection
from rxncon.core.reaction import Reaction
from rxncon.core.state import State
from typecheck import typecheck
from typing import List, Union

class BooleanModel:
    @typecheck
    def __init__(self, update_rules: List['UpdateRule'], initial_conditions: List['InitialCondition']):
        self.update_rules = update_rules
        self.initial_conditions = initial_conditions
        self._validate()

    def _validate(self):
        pass


class InitialCondition:
    @typecheck
    def __init__(self, target: 'Target', value: bool):
        self.target = target
        self.value = value

    @typecheck
    def __eq__(self, other: 'InitialCondition'):
        return self.target == other.target and self.value == self.value

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "target: {0}, value: {1}".format(self.target, self.value)


class Target:
    def __hash__(self):
        return hash(str(self))

    def __repr__(self):
        return str(self)

    def __str__(self):
        return str(self.value)


class ReactionTarget(Target):
    @typecheck
    def __init__(self, reaction_parent: Reaction):
        self.reaction_parent = reaction_parent
        self.produced_targets    = [StateTarget(x) for x in reaction_parent.produced_states]
        self.consumed_targets    = [StateTarget(x) for x in reaction_parent.consumed_states]
        self.synthesised_targets = [StateTarget(x) for x in reaction_parent.synthesised_states]
        self.degraded_targets    = [StateTarget(x) for x in reaction_parent.degraded_states]

    @typecheck
    def produces(self, state_target: StateTarget):
        return state_target in self.produced_targets

    @typecheck
    def consumes(self, state_target: StateTarget):
        return state_target in self.consumed_targets

    @typecheck
    def synthesises(self, state_target: StateTarget):
        return state_target in self.synthesised_targets

    @typecheck
    def degrades(self, state_target: StateTarget):
        return state_target in self.degraded_targets


class StateTarget(Target):
    def __init__(self, state_parent: State):
        self.state_parent = state_parent

    def __eq__(self, other: 'StateTarget'):
        return self.state_parent == other.state_parent

    @typecheck
    def is_produced_by(self, reaction_target: ReactionTarget):
        return reaction_target.produces(self)

    def is_consumed_by(self, reaction_target: ReactionTarget):
        return reaction_target.consumes(self)

    def is_synthesised_by(self, reaction_target: ReactionTarget):
        return reaction_target.synthesises(self)

    def is_degraded_by(self, reaction_target: ReactionTarget):
        return reaction_target.degrades(self)

class UpdateRule:
    @typecheck
    def __init__(self, target: Target, factor: VennSet):
        self.target = target
        self.factor = factor
        self._validate()

    def _validate(self):
        assert self.factor.value_type == Target

    def __str__(self):
        return "target: {0}, factors: {1}".format(self.target, self.factor)
