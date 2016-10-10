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


class BinaryEffector(Effector, metaclass=ABCMeta):
    def __init__(self, left_expr: Effector, right_expr: Effector):
        self.left_expr = left_expr
        self.right_expr = right_expr

    def __eq__(self, other: Effector) -> bool:
        return type(self) == type(other) and self.name == other.name and \
            self.left_expr == other.left_expr and self.right_expr == other.right_expr

    @property
    def states(self) -> List[State]:
        return self.left_expr.states + self.right_expr.states


class AndEffector(BinaryEffector):
    def __str__(self) -> str:
        if self.name:
            return 'AndEffector{0}({1}, {2})'.format(self.name, self.left_expr, self.right_expr)
        else:
            return 'AndEffector({0}, {1})'.format(self.left_expr, self.right_expr)


class OrEffector(BinaryEffector):
    def __str__(self) -> str:
        if self.name:
            return 'OrEffector{0}({1}, {2})'.format(self.name, self.left_expr, self.right_expr)
        else:
            return 'OrEffector({0}, {1})'.format(self.left_expr, self.right_expr)
