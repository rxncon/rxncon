import functools
import operator
from typing import Dict, List, Any, Optional
from collections import OrderedDict
from itertools import product
import re

from pyeda.inter import And, Or, Not, expr
from pyeda.boolalg.expr import AndOp, OrOp, NotOp, Variable, Implies, Expression, Literal, NaryOp, \
    Complement as pyedaComplement

SYMS = [''.join(tup) for tup in product('ABCDEFGHIJKLMNOPQRSTUVWXYZ', repeat=2)]

class Set:
    def calc_solutions(self) -> List[Dict['Set', bool]]:
        val_to_sym = self._make_val_to_sym_dict()
        sym_to_val = {sym: val for val, sym in val_to_sym.items()}

        venn_solns = []
        for s in self._to_pyeda_expr(val_to_sym).satisfy_all():
            venn_solns.append({sym_to_val[sym.name]: bool(truth) for sym, truth in s.items()})

        return venn_solns

    def to_simplified_set(self) -> 'Set':
        val_to_sym = self._make_val_to_sym_dict()
        sym_to_val = {sym: val for val, sym in val_to_sym.items()}
        return venn_from_pyeda(self._to_pyeda_expr(val_to_sym).simplify(), sym_to_val)

    def to_dnf_set(self) -> 'Set':
        val_to_sym = self._make_val_to_sym_dict()
        sym_to_val = {sym: val for val, sym in val_to_sym.items()}
        return venn_from_pyeda(self._to_pyeda_expr(val_to_sym).to_dnf(), sym_to_val)

    def to_dnf_list(self) -> List['Set']:
        val_to_sym = self._make_val_to_sym_dict()
        sym_to_val = {sym: val for val, sym in val_to_sym.items()}
        dnf_set = self._to_pyeda_expr(val_to_sym).to_dnf()

        if isinstance(dnf_set, Literal):
            return [venn_from_pyeda(dnf_set, sym_to_val)]
        elif isinstance(dnf_set, NaryOp):
            return [venn_from_pyeda(x, sym_to_val) for x in dnf_set.xs]
        else:
            raise Exception

    def to_dnf_nested_list(self) -> List[List['Set']]:
        val_to_sym = self._make_val_to_sym_dict()
        sym_to_val = {sym: val for val, sym in val_to_sym.items()}
        dnf_set = self._to_pyeda_expr(val_to_sym).to_dnf()

        if isinstance(dnf_set, Literal):
            return [venn_from_pyeda(dnf_set, sym_to_val)]

        res = []
        for term in dnf_set.xs:
            if isinstance(term, Literal):
                res.append([venn_from_pyeda(term, sym_to_val)])
            elif isinstance(term, NaryOp):
                res.append([venn_from_pyeda(x, sym_to_val) for x in term.xs])
            else:
                raise Exception
        return res

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
        # assert all(types[0] in x.mro() for x in types)
        return types[0]

    def _to_pyeda_expr(self, val_to_sym: Dict[Any, str]) -> Expression:
        return None

    def _make_val_to_sym_dict(self, existing_dict: Optional[OrderedDict]=None) -> OrderedDict:
        vals = list(set(self.values))

        if existing_dict:
            d = existing_dict
            first_new_sym_idx = SYMS.index(next(reversed(d.values()))) + 1
            for i, val in enumerate(x for x in vals if x not in d.keys()):
                d[val] = SYMS[i + first_new_sym_idx]
        else:
            d = OrderedDict()
            for i, val in enumerate(vals):
                d[val] = SYMS[i]

        return d


class UnarySet(Set):
    pass


class ValueSet(UnarySet):
    def __init__(self, value):
        assert isinstance(hash(value), int)
        self.value = value

    def __eq__(self, other: Set) -> bool:
        if isinstance(other, ValueSet):
            return self.value == other.value
        else:
            return False

    def __hash__(self) -> int:
        return hash('*value-set-{}*'.format(hash(self.value)))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        if self.value:
            return '({})'.format(str(self.value))
        else:
            return 'UniversalSet'

    @property
    def values(self):
        return [self.value]

    def _to_pyeda_expr(self, val_to_sym: OrderedDict) -> Expression:
        if self.value:
            return expr(val_to_sym[self.value])
        else:
            return expr(1)


class EmptySet(UnarySet):
    def __eq__(self, other: Set) -> bool:
        return isinstance(other, EmptySet)

    def __hash__(self) -> int:
        return hash('*empty-set*')

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'EmptySet'

    def _to_pyeda_expr(self, val_to_sym: OrderedDict) -> Expression:
        return expr(0)


class Complement(UnarySet):
    def __init__(self, expr: Set):
        self.expr = expr

    def __eq__(self, other: Set) -> bool:
        if isinstance(other, Complement):
            return self.expr == other.expr
        else:
            return False

    def __hash__(self) -> int:
        return hash('*complement-{}*'.format(self.expr))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return '!({})'.format(str(self.expr))

    @property
    def values(self):
        return self.expr.values

    def _to_pyeda_expr(self, val_to_sym: OrderedDict) -> Expression:
        return Not(self.expr._to_pyeda_expr(val_to_sym))


class NarySet(Set):
    def __init__(self, *exprs: Set):
        self.exprs = exprs

    @property
    def values(self):
        return functools.reduce(operator.add, [expr.values for expr in self.exprs], [])


class Intersection(NarySet):
    def __eq__(self, other: Set) -> bool:
        if isinstance(other, Intersection):
            return all(my_expr == other_expr for my_expr, other_expr in zip(self.exprs, other.exprs))
        else:
            return False

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return '({})'.format(' & '.join(str(expr) for expr in self.exprs))

    def _to_pyeda_expr(self, val_to_sym: OrderedDict) -> Expression:
        return And(*(expr._to_pyeda_expr(val_to_sym) for expr in self.exprs))


class Union(NarySet):
    def __eq__(self, other: Set) -> bool:
        if isinstance(other, Union):
            return all(my_expr == other_expr for my_expr, other_expr in zip(self.exprs, other.exprs))
        else:
            return False

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return '({})'.format(' | '.join(str(expr) for expr in self.exprs))

    def _to_pyeda_expr(self, val_to_sym: OrderedDict) -> Expression:
        return Or(*(expr._to_pyeda_expr(val_to_sym) for expr in self.exprs))


class Difference(Set):
    def __new__(cls, *args, **kwargs) -> Set:
        assert len(args) == 2
        return Intersection(args[0], Complement(args[1]))


def UniversalSet() -> Set:
    return ValueSet(None)


def venn_from_pyeda(pyeda_expr, sym_to_val):
    if isinstance(pyeda_expr, Variable):
        return ValueSet(sym_to_val[pyeda_expr.name])
    elif isinstance(pyeda_expr, AndOp):
        return Intersection(*(venn_from_pyeda(x, sym_to_val) for x in pyeda_expr.xs))
    elif isinstance(pyeda_expr, OrOp):
        return Union(*(venn_from_pyeda(x, sym_to_val) for x in pyeda_expr.xs))
    elif isinstance(pyeda_expr, NotOp):
        return Complement(venn_from_pyeda(pyeda_expr.x, sym_to_val))
    elif isinstance(pyeda_expr, pyedaComplement):
        return Complement(ValueSet(sym_to_val[pyeda_expr.inputs[0].name]))
    else:
        raise Exception


def venn_from_str(venn_str, value_parser=lambda x: x):
    # The values have to be surrounded by a single space.
    BOOL_REGEX            = '[\(\)\|\&\~]+'
    pyeda_str             = ''
    pyeda_sym_to_val      = {}
    venn_sym_to_pyeda_sym = {}
    current_sym           = 0

    parts = venn_str.split()

    for part in parts:
        match = re.match(BOOL_REGEX, part)
        if match and match.string == part:
            pyeda_str += part
        else:
            if part in venn_sym_to_pyeda_sym.keys():
                pyeda_str += venn_sym_to_pyeda_sym[part]
            else:
                pyeda_sym = SYMS[current_sym]
                current_sym += 1
                venn_sym_to_pyeda_sym[part] = pyeda_sym
                pyeda_sym_to_val[pyeda_sym] = value_parser(part)
                pyeda_str += pyeda_sym

    return venn_from_pyeda(expr(pyeda_str), pyeda_sym_to_val)
