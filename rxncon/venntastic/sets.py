import functools
import itertools as itt
from typing import Dict, List

from pyeda.inter import And, Or, Not, expr
from pyeda.boolalg.expr import AndOp, OrOp, NotOp, Variable, Implies

SYMS = 'abcdefghijklmnopqrstuvwxyz'

class Set:
    def to_full_simplified_form(self) -> 'Set':
        val_to_sym = self._make_val_to_sym_dict()
        sym_to_val = {sym: val for val, sym in val_to_sym.items()}
        return pyeda_to_venn(self._to_pyeda_expr(val_to_sym).to_cnf().to_dnf(), sym_to_val)

    def is_equivalent_to(self, other: 'Set') -> bool:
        val_to_sym = self._make_val_to_sym_dict()
        val_to_sym = other._make_val_to_sym_dict(val_to_sym)
        return self._to_pyeda_expr(val_to_sym).equivalent(other._to_pyeda_expr(val_to_sym))

    def is_subset_of(self, other: 'Set') -> bool:
        val_to_sym = self._make_val_to_sym_dict()
        val_to_sym = other._make_val_to_sym_dict(val_to_sym)
        return Implies(self._to_pyeda_expr(val_to_sym), other._to_pyeda_expr(val_to_sym)).equivalent(expr(1))

    def is_superset_of(self, other: 'Set') -> bool:
        return other.is_subset_of(self)

    @property
    def values(self):
        return []

    @property
    def value_type(self):
        types = [type(x) for x in self.values]
        # assert all(types[0] == x for x in types)
        return types[0]

    def _to_pyeda_expr(self, val_to_sym):
        return None

    def _make_val_to_sym_dict(self, existing_dict=None):
        vals = list(set(self.values))

        if existing_dict:
            d = existing_dict
            first_new_sym_idx = ord(max(d.values())) - 96
            for i, val in enumerate(x for x in vals if x not in d.keys()):
                d[val] = SYMS[i + first_new_sym_idx]
        else:
            d = {}
            for i, val in enumerate(vals):
                d[val] = SYMS[i]

        return d


class UnarySet(Set):
    pass


class ValueSet(UnarySet):
    def __init__(self, value):
        assert isinstance(hash(value), int)
        self.value = value

    def __eq__(self, other: Set):
        if isinstance(other, ValueSet):
            return self.value == other.value
        else:
            return False

    def __hash__(self) -> int:
        return hash('*property-set-{}*'.format(hash(self.value)))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        if self.value:
            return 'ValueSet({})'.format(self.value)
        else:
            return 'UniversalSet'

    @property
    def values(self):
        return [self.value]

    def _to_pyeda_expr(self, val_to_sym):
        if self.value:
            return expr(val_to_sym[self.value])
        else:
            return expr(1)


class EmptySet(UnarySet):
    def __eq__(self, other: Set):
        return isinstance(other, EmptySet)

    def __hash__(self) -> int:
        return hash('*empty-set*')

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'EmptySet'

    def _to_pyeda_expr(self, val_to_sym):
        return expr(0)


class Complement(UnarySet):
    def __init__(self, expr: Set):
        self.expr = expr

    def __eq__(self, other: Set):
        if isinstance(other, Complement):
            return self.expr == other.expr

    def __hash__(self) -> int:
        return hash('*complement-{}*'.format(self.expr))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'Complement({})'.format(self.expr)

    @property
    def values(self):
        return self.expr.values

    def _to_pyeda_expr(self, val_to_sym):
        return Not(self.expr._to_pyeda_expr(val_to_sym))


class BinarySet(Set):
    def __init__(self, left_expr: Set, right_expr: Set):
        self.left_expr = left_expr
        self.right_expr = right_expr

    @property
    def values(self):
        return self.left_expr.values + self.right_expr.values


class Intersection(BinarySet):
    def __eq__(self, other: Set):
        if isinstance(other, Intersection):
            return self.left_expr == other.left_expr and self.right_expr == other.right_expr
        else:
            return False

    def __hash__(self) -> int:
        return hash('*intersection-{0}{1}*'.format(hash(self.left_expr), hash(self.right_expr)))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'Intersection({0}, {1})'.format(self.left_expr, self.right_expr)

    def _to_pyeda_expr(self, val_to_sym):
        return And(self.left_expr._to_pyeda_expr(val_to_sym), self.right_expr._to_pyeda_expr(val_to_sym))


class Union(BinarySet):
    def __eq__(self, other: Set):
        if isinstance(other, Union):
            return self.left_expr == other.left_expr and self.right_expr == other.right_expr
        else:
            return False

    def __hash__(self) -> int:
        return hash('*union-{0}{1}*'.format(hash(self.left_expr), hash(self.right_expr)))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'Union({0}, {1})'.format(self.left_expr, self.right_expr)

    def _to_pyeda_expr(self, val_to_sym):
        return Or(self.left_expr._to_pyeda_expr(val_to_sym), self.right_expr._to_pyeda_expr(val_to_sym))


class Difference(Set):
    def __new__(cls, *args, **kwargs) -> Set:
        assert len(args) == 2
        return Intersection(args[0], Complement(args[1]))


def UniversalSet() -> Set:
    return ValueSet(None)


def MultiUnion(*args):
    return nested_expression_from_list_and_binary_op(list(args), Union)


def MultiIntersection(*args):
    return nested_expression_from_list_and_binary_op(list(args), Intersection)


### PUBLIC FUNCTIONS ###
def nested_expression_from_list_and_binary_op(xs: List[Set], binary_op) -> Set:
    if binary_op == Intersection:
        unit = UniversalSet()
    elif binary_op == Union:
        unit = EmptySet()
    else:
        raise TypeError

    if len(xs) == 0:
        return unit
    elif len(xs) == 1:
        return xs[0]
    else:
        return functools.reduce(binary_op, xs[1:], xs[0])


def pyeda_to_venn(pyeda_expr, sym_to_val):
    if isinstance(pyeda_expr, Variable):
        return ValueSet(sym_to_val[pyeda_expr.name])
    elif isinstance(pyeda_expr, AndOp):
        return MultiIntersection(*(pyeda_to_venn(x, sym_to_val) for x in pyeda_expr.xs))
    elif isinstance(pyeda_expr, OrOp):
        return MultiUnion(*(pyeda_to_venn(x, sym_to_val) for x in pyeda_expr.xs))
    elif isinstance(pyeda_expr, NotOp):
        return Complement(pyeda_to_venn(pyeda_expr.x, sym_to_val))
    else:
        raise Exception


