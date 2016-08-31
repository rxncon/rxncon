from rxncon.venntastic.sets import Set as VennSet
from rxncon.core.reaction import Reaction
from rxncon.core.state import State
from rxncon.core.spec import Spec
from typecheck import typecheck
from typing import List, Union

class BooleanModel:
    @typecheck
    def __init__(self, rules: List["Rule"], init_conditions: List['InitialCondition']):
        self.rules = rules
        self.init_conditions = init_conditions
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
    @typecheck
    def __init__(self, value: Union[Reaction, State]):
        self.value = value

    @typecheck
    def __eq__(self, other: 'Target'):
        if isinstance(self.value, Reaction) and isinstance(other.value, Reaction) and self.value == other.value:
            return True
        elif isinstance(self.value, State) and isinstance(other.value, State) and self.value == other.value:
            return True
        elif isinstance(self.value, Spec) and isinstance(other.value, Spec) and self.value == other.value:
            return True
        else:
            return False

    def __hash__(self):
        return hash(str(self))

    def __repr__(self):
        return str(self)

    def __str__(self):
        return str(self.value)


class Rule:
    @typecheck
    def __init__(self, target: Target, factor: VennSet):
        self.target = target
        self.factor = factor
        self._validate()

    def _validate(self):
        assert self.target.value is not None
        assert isinstance(self.target, Target)
        assert self.factor.value is not None
        assert isinstance(self.factor, Factor)

    def __str__(self):
        return "target: {0}, factors: {1}".format(self.target, self.factor)
