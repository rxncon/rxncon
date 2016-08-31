from rxncon.venntastic.sets import Set as VennSet
from rxncon.core.reaction import Reaction
from rxncon.core.state import State
from rxncon.core.spec import MolSpec
from typecheck import typecheck
from typing import List

class BooleanModel:
    @typecheck
    def __init__(self, update_rules: List['UpdateRule'], initial_conditions: List['InitialCondition']):
        self.update_rules = update_rules
        self.initial_conditions = initial_conditions
        self._validate()

    def _validate(self):
        all_lhs_targets = []
        all_rhs_targets = []
        for rule in self.update_rules:
            all_lhs_targets.append(rule.target)
            all_rhs_targets += rule.factor_targets

        assert all(x in all_lhs_targets for x in all_rhs_targets)


class InitialCondition:
    @typecheck
    def __init__(self, target: 'Target', value: bool):
        self.target = target
        self.value = value

    @typecheck
    def __eq__(self, other: 'InitialCondition') -> bool:
        return self.target == other.target and self.value == self.value

    @typecheck
    def __repr__(self) -> str:
        return str(self)

    @typecheck
    def __str__(self) -> str:
        return "target: {0}, value: {1}".format(self.target, self.value)


class Target:
    @typecheck
    def __hash__(self) -> int:
        return hash(str(self))

    @typecheck
    def __repr__(self) -> str:
        return str(self)


class ReactionTarget(Target):
    @typecheck
    def __init__(self, reaction_parent: Reaction):
        self._reaction_parent = reaction_parent
        self.produced_targets    = [StateTarget(x) for x in reaction_parent.produced_states]
        self.consumed_targets    = [StateTarget(x) for x in reaction_parent.consumed_states]
        self.synthesised_targets = [StateTarget(x) for x in reaction_parent.synthesised_states]
        self.degraded_targets    = [StateTarget(x) for x in reaction_parent.degraded_states]

    @typecheck
    def __eq__(self, other: 'ReactionTarget'):
        #  Possibly more than one ReactionTarget from a single reaction_parent, so also check all its targets.
        return self._reaction_parent == other._reaction_parent and self.produced_targets == other.produced_targets and \
            self.consumed_targets == other.consumed_targets and self.synthesised_targets == other.synthesised_targets and \
            self.degraded_targets == other.degraded_targets

    @typecheck
    def __str__(self) -> str:
        return str(self._reaction_parent)

    @typecheck
    def produces(self, state_target: StateTarget) -> bool:
        return state_target in self.produced_targets

    @typecheck
    def consumes(self, state_target: StateTarget) -> bool:
        return state_target in self.consumed_targets

    @typecheck
    def synthesises(self, state_target: StateTarget) -> bool:
        return state_target in self.synthesised_targets

    @typecheck
    def degrades(self, state_target: StateTarget) -> bool:
        return state_target in self.degraded_targets

    @typecheck
    @property
    def components(self) -> List[MolSpec]:
        return list(set(self._reaction_parent.components_lhs + self._reaction_parent.components_rhs))


class StateTarget(Target):
    @typecheck
    def __init__(self, state_parent: State):
        self._state_parent = state_parent

    @typecheck
    def __eq__(self, other: 'StateTarget') -> bool:
        return self._state_parent == other._state_parent

    @typecheck
    def is_produced_by(self, reaction_target: ReactionTarget) -> bool:
        return reaction_target.produces(self)

    @typecheck
    def is_consumed_by(self, reaction_target: ReactionTarget) -> bool:
        return reaction_target.consumes(self)

    @typecheck
    def is_synthesised_by(self, reaction_target: ReactionTarget) -> bool:
        return reaction_target.synthesises(self)

    @typecheck
    def is_degraded_by(self, reaction_target: ReactionTarget) -> bool:
        return reaction_target.degrades(self)

    @typecheck
    @property
    def components(self) -> List[MolSpec]:
        return self._state_parent.components

class UpdateRule:
    @typecheck
    def __init__(self, target: Target, factor: VennSet):
        self.target = target
        self.factor = factor
        self._validate()

    @typecheck
    def factor_targets(self) -> List[Target]:
        return self.factor.values

    def _validate(self):
        assert self.factor.value_type == Target

    def __str__(self):
        return "target: {0}, factors: {1}".format(self.target, self.factor)
