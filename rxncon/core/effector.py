from abc import ABCMeta, abstractproperty
from typing import List, Optional

from rxncon.core.state import State


class Effector(metaclass=ABCMeta):
    @property
    def name(self) -> Optional[str]:
        try:
            return self._name
        except AttributeError:
            return None

    @name.setter
    def name(self, value: str):
        self._name = value

    @abstractproperty
    def states(self) -> List[State]:
        pass


class StateEffector(Effector):
    def __init__(self, expr: State):
        self.expr = expr

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'StateEffector({})'.format(str(self.expr))

    def __eq__(self, other: Effector) -> bool:
        return isinstance(other, StateEffector) and self.expr == other.expr and self.name == other.name

    @property
    def states(self) -> List[State]:
        return [self.expr]


class NotEffector(Effector):
    def __init__(self, expr: Effector):
        self.expr = expr

    def __str__(self) -> str:
        return 'NotEffector({})'.format(self.expr)

    def __eq__(self, other: Effector) -> bool:
        return isinstance(other, NotEffector) and self.expr == other.expr and self.name == other.name

    @property
    def states(self) -> List[State]:
        return self.expr.states


class AndEffector(Effector):
    def __init__(self, *exprs):
        self.exprs = exprs

    def __str__(self) -> str:
        if self.name:
            return 'AndEffector{0}({1})'.format(self.name, ','.join(str(x) for x in self.exprs))
        else:
            return 'AndEffector({0})'.format(','.join(str(x) for x in self.exprs))

    def __eq__(self, other: Effector) -> bool:
        return isinstance(other, AndEffector) and self.name == other.name and \
               self.exprs == other.exprs

    @property
    def states(self):
        return [state for x in self.exprs for state in x.states]


class OrEffector(Effector):
    def __init__(self, *exprs):
        self.exprs = exprs

    def __str__(self) -> str:
        if self.name:
            return 'OrEffector{0}({1})'.format(self.name, ','.join(str(x) for x in self.exprs))
        else:
            return 'OrEffector({0})'.format(','.join(str(x) for x in self.exprs))

    def __eq__(self, other: Effector) -> bool:
        return isinstance(other, OrEffector) and self.name == other.name and \
               self.exprs == other.exprs

    @property
    def states(self):
        return [state for x in self.exprs for state in x.states]
