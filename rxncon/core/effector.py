import re
from abc import ABCMeta, abstractproperty
from typing import List, Optional, Tuple

from rxncon.core.spec import spec_from_str
from rxncon.core.state import State
from rxncon.util.utils import OrderedEnum


BOOLEAN_CONTINGENCY_REGEX = '^<.*>$'


class BooleanOperator(OrderedEnum):
    op_and = 'and'
    op_or  = 'or'
    op_not = 'not'
    op_eqv = 'eqv'


class BooleanContingencyName:
    def __init__(self, name: str):
        assert re.match(BOOLEAN_CONTINGENCY_REGEX, name)
        self.name = name

    def __eq__(self, other: 'BooleanContingencyName') -> bool:
        return self.name == other.name

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return self.name


class QualifiedSpec:
    def __init__(self, qualified_spec_str: str):
        self.namespace = [BooleanContingencyName(x) for x in qualified_spec_str.split('.')[:-1]]
        self.spec      = spec_from_str(qualified_spec_str.split('.')[-1])
        self._name     = qualified_spec_str


    def __str__(self) -> str:
        return self._name

    def __repr__(self) -> str:
        return 'QualifiedSpec<{}>'.format(self._name)

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

    def to_flattened(self) -> 'Effector':
        pass

    @property
    def is_leaf(self) -> bool:
        raise NotImplementedError


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

    def to_flattened(self):
        return self

    @property
    def is_leaf(self) -> bool:
        return True


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

    def to_flattened(self):
        return NotEffector(self.expr.to_flattened())

    @property
    def is_leaf(self):
        return self.expr.is_leaf


class NaryEffector(Effector):
    def __init__(self, *exprs):
        self.exprs        = exprs
        self.equivalences = []     # type: List[Tuple[QualifiedSpec, QualifiedSpec]]

    @property
    def states(self) -> List[State]:
        return [state for x in self.exprs for state in x.states]

    def to_flattened(self):
        pass

    @property
    def is_leaf(self) -> bool:
        return False

class AndEffector(NaryEffector):
    def __str__(self) -> str:
        if self.name:
            return 'AndEffector{0}({1})'.format(self.name, ','.join(str(x) for x in self.exprs))
        else:
            return 'AndEffector({0})'.format(','.join(str(x) for x in self.exprs))

    def __eq__(self, other: Effector) -> bool:
        return isinstance(other, AndEffector) and self.name == other.name and \
               self.exprs == other.exprs


class OrEffector(NaryEffector):
    def __str__(self) -> str:
        if self.name:
            return 'OrEffector{0}({1})'.format(self.name, ','.join(str(x) for x in self.exprs))
        else:
            return 'OrEffector({0})'.format(','.join(str(x) for x in self.exprs))

    def __eq__(self, other: Effector) -> bool:
        return isinstance(other, OrEffector) and self.name == other.name and \
               self.exprs == other.exprs


