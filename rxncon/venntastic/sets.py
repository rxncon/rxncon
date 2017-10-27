import functools
import operator
from typing import Dict, List, Generic, Optional, TypeVar, MutableMapping, Callable, Any
from collections import OrderedDict
from itertools import product
import re
from copy import deepcopy

from pyeda.inter import And, Or, Not, Xor, expr
from pyeda.boolalg.expr import XorOp, AndOp, OrOp, NotOp, Variable, Implies, Expression, Literal, \
    Complement as pyedaComplement, One, Zero

SYMS = [''.join(tup) for tup in product('ABCDEFGHIJKLMNOPQRSTUVWXYZ', repeat=2)]

# Since all Set expressions except ValueSet are covariant, we make the 'T' type var covariant,
# and make ValueSet[T_inv] for invariant.
T = TypeVar('T', covariant=True)
T_inv = TypeVar('T_inv')


class Set(Generic[T]):
    def calc_solutions(self) -> List[Dict[T, bool]]:
        val_to_sym = self._make_val_to_sym_dict()
        sym_to_val = {sym: val for val, sym in val_to_sym.items()}

        venn_solns = []
        for s in self._to_pyeda_expr(val_to_sym).satisfy_all():
            venn_solns.append({sym_to_val[sym.name]: bool(truth) for sym, truth in s.items()})

        return venn_solns

    def eval_boolean_func(self, vars: Dict[T, bool]) -> bool:
        val_to_sym = self._make_val_to_sym_dict()
        try:
            evaluated = self._to_pyeda_expr(val_to_sym).restrict({expr(sym): vars[val] for val, sym in val_to_sym.items()})
        except KeyError as e:
            raise AssertionError('eval_boolean_func missing variable {}'.format(e.args[0]))

        if evaluated.is_one():
            return True
        elif evaluated.is_zero():
            return False
        else:
            raise AssertionError

    def to_simplified_set(self) -> 'Set[T]':
        val_to_sym = self._make_val_to_sym_dict()
        sym_to_val = {sym: val for val, sym in val_to_sym.items()}
        return venn_from_pyeda(self._to_pyeda_expr(val_to_sym).simplify(), sym_to_val)

    def to_dnf_set(self) -> 'Set[T]':
        val_to_sym = self._make_val_to_sym_dict()
        sym_to_val = {sym: val for val, sym in val_to_sym.items()}
        return venn_from_pyeda(self._to_pyeda_expr(val_to_sym).to_dnf(), sym_to_val)

    def to_dnf_list(self) -> List['Set[T]']:
        val_to_sym = self._make_val_to_sym_dict()
        sym_to_val = {sym: val for val, sym in val_to_sym.items()}
        dnf_set = self._to_pyeda_expr(val_to_sym).to_dnf()

        if dnf_set is One:
            return [UniversalSet()]
        elif dnf_set is Zero:
            return [EmptySet()]
        elif isinstance(dnf_set, Literal):
            return [venn_from_pyeda(dnf_set, sym_to_val)]
        elif isinstance(dnf_set, AndOp):
            return [venn_from_pyeda(dnf_set, sym_to_val)]
        elif isinstance(dnf_set, OrOp):
            return [venn_from_pyeda(term, sym_to_val) for term in dnf_set.xs]
        else:
            raise Exception

    def to_dnf_nested_list(self) -> List[List['Set[T]']]:
        val_to_sym = self._make_val_to_sym_dict()
        sym_to_val = {sym: val for val, sym in val_to_sym.items()}
        dnf_set = self._to_pyeda_expr(val_to_sym).to_dnf()

        if dnf_set is One:
            return [[UniversalSet()]]
        elif dnf_set is Zero:
            return [[EmptySet()]]
        elif isinstance(dnf_set, Literal):
            return [[venn_from_pyeda(dnf_set, sym_to_val)]]
        elif isinstance(dnf_set, AndOp):
            return [[venn_from_pyeda(x, sym_to_val) for x in dnf_set.xs]]

        res = []  # type: List[List[Set[T]]]
        for term in dnf_set.xs:
            if isinstance(term, Literal):
                res.append([venn_from_pyeda(term, sym_to_val)])
            elif isinstance(term, AndOp):
                res.append([venn_from_pyeda(x, sym_to_val) for x in term.xs])
            else:
                raise Exception
        return res

    def is_equivalent_to(self, other: 'Set[T]') -> bool:
        val_to_sym = self._make_val_to_sym_dict()
        val_to_sym = other._make_val_to_sym_dict(val_to_sym)
        return self._to_pyeda_expr(val_to_sym).equivalent(other._to_pyeda_expr(val_to_sym))

    def is_subset_of(self, other: 'Set[T]') -> bool:
        val_to_sym = self._make_val_to_sym_dict()
        val_to_sym = other._make_val_to_sym_dict(val_to_sym)
        return Implies(self._to_pyeda_expr(val_to_sym), other._to_pyeda_expr(val_to_sym)).equivalent(expr(1))

    def is_superset_of(self, other: 'Set[T]') -> bool:
        return other.is_subset_of(self)

    @property
    def values(self) -> List[T]:
        return []

    def _to_pyeda_expr(self, val_to_sym: MutableMapping[T, str]) -> Expression:
        return None

    def _make_val_to_sym_dict(self, existing_dict: Optional[MutableMapping[T, str]]=None) -> MutableMapping[T, str]:
        vals = []  # type: List[T]

        for v in self.values:
            found = False
            for existing in vals:
                try:
                    if v == existing:
                        found = True
                        break
                except AttributeError:
                    pass

            if not found:
                vals.append(v)

        if existing_dict:
            d = existing_dict
            first_new_sym_idx = SYMS.index(next(reversed(d.values()))) + 1  # type: ignore
            for i, val in enumerate(x for x in vals if x not in d.keys()):
                d[val] = SYMS[i + first_new_sym_idx]
        else:
            d = OrderedDict()
            for i, val in enumerate(vals):
                d[val] = SYMS[i]

        return d


class EmptySet(Set[Any]):
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Set):
            return NotImplemented
        return isinstance(other, EmptySet)

    def __hash__(self) -> int:
        return hash('*empty-set*')

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'EmptySet'

    def _to_pyeda_expr(self, val_to_sym: MutableMapping[Any, str]) -> Expression:
        return expr(0)


class UniversalSet(Set[Any]):
    def __init__(self) -> None:
        pass

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Set):
            return NotImplemented
        return isinstance(other, UniversalSet)

    def __hash__(self) -> int:
        return hash('UniversalSet')

    def __repr__(self) -> str:
        return 'UniversalSet'

    def __str__(self) -> str:
        return 'UniversalSet'

    def _to_pyeda_expr(self, val_to_sym: MutableMapping[Any, str]) -> Expression:
        return expr(1)


class UnarySet(Set[T], Generic[T]):
    pass


class ValueSet(UnarySet[T_inv], Generic[T_inv]):
    def __init__(self, value: T_inv) -> None:
        assert value is not None
        assert isinstance(hash(value), int)
        self.value = value

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Set):
            return NotImplemented
        elif isinstance(other, ValueSet):
            return self.value == other.value
        else:
            return False

    def __hash__(self) -> int:
        return hash('*value-set-{}*'.format(hash(self.value)))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        if self.value:
            return '{}'.format(str(self.value))
        else:
            return 'UniversalSet'

    @property
    def values(self) -> List[T_inv]:
        return [self.value]

    def _to_pyeda_expr(self, val_to_sym: MutableMapping[Any, str]) -> Expression:
        return expr(val_to_sym[self.value])


class Complement(UnarySet[T], Generic[T]):
    def __init__(self, expr: Set[T]) -> None:
        self.expr = expr

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Set):
            return NotImplemented
        elif isinstance(other, Complement):
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
    def values(self) -> List[T]:
        return self.expr.values

    def _to_pyeda_expr(self, val_to_sym: MutableMapping[T, str]) -> Expression:
        return Not(self.expr._to_pyeda_expr(val_to_sym))


class NarySet(Set[T], Generic[T]):
    def __new__(cls, *exprs: Set[T]) -> Set[T]:
        assert len(exprs) > 0
        if len(exprs) == 1:
            return deepcopy(exprs[0])
        else:
            return super().__new__(cls, *exprs)  # type: ignore

    def __init__(self, *exprs: Set[T]) -> None:
        self.exprs = exprs

    @property
    def values(self) -> List[T]:
        return functools.reduce(operator.add, [expr.values for expr in self.exprs], [])


class Intersection(NarySet[T], Generic[T]):
    def __new__(cls, *exprs: Set[T], **kwargs: Any) -> Set[T]:
        if len(exprs) == 0:
            return UniversalSet()
        else:
            return NarySet.__new__(cls, *exprs)  # type: ignore

    def __deepcopy__(self, memodict: Dict) -> Set[T]:
        return Intersection(*(deepcopy(x) for x in self.exprs))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Set):
            return NotImplemented
        elif isinstance(other, Intersection):
            return all(my_expr == other_expr for my_expr, other_expr in zip(self.exprs, other.exprs))
        else:
            return False

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return '({})'.format(' & '.join(str(expr) for expr in self.exprs))

    def _to_pyeda_expr(self, val_to_sym: MutableMapping[T, str]) -> Expression:
        return And(*(expr._to_pyeda_expr(val_to_sym) for expr in self.exprs))


class Union(NarySet[T], Generic[T]):
    def __new__(cls, *exprs: Set[T], **kwargs: Any) -> Set[T]:
        if len(exprs) == 0:
            return EmptySet()
        else:
            return NarySet.__new__(cls, *exprs)  # type: ignore

    def __deepcopy__(self, memodict: Dict) -> Set[T]:
        return Union(*(deepcopy(x) for x in self.exprs))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Set):
            return NotImplemented
        elif isinstance(other, Union):
            return all(my_expr == other_expr for my_expr, other_expr in zip(self.exprs, other.exprs))
        else:
            return False

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return '({})'.format(' | '.join(str(expr) for expr in self.exprs))

    def _to_pyeda_expr(self, val_to_sym: MutableMapping[T, str]) -> Expression:
        return Or(*(expr._to_pyeda_expr(val_to_sym) for expr in self.exprs))


class DisjunctiveUnion(NarySet[T], Generic[T]):
    def __new__(cls, *exprs: Set[T], **kwargs: Any) -> Set[T]:
        if len(exprs) == 0:
            return EmptySet()
        else:
            return super().__new__(cls, *exprs)  # type: ignore

    def __deepcopy__(self, memodict: Dict) -> Set[T]:
        return DisjunctiveUnion(*(deepcopy(x) for x in self.exprs))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Set):
            return NotImplemented
        elif isinstance(other, DisjunctiveUnion):
            return all(my_expr == other_expr for my_expr, other_expr in zip(self.exprs, other.exprs))
        else:
            return False

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return '({})'.format(' XOR '.join(str(expr) for expr in self.exprs))

    def _to_pyeda_expr(self, val_to_sym: MutableMapping[T, str]) -> Expression:
        return Xor(*(expr._to_pyeda_expr(val_to_sym) for expr in self.exprs))


class Difference(Set[T], Generic[T]):
    def __new__(cls, *args: Set[T]) -> Set[T]:
        assert len(args) == 2
        return Intersection(args[0], Complement(args[1]))


def venn_from_pyeda(pyeda_expr: Expression, sym_to_val: MutableMapping[str, T]) -> Set[T]:
    if pyeda_expr is One:
        return UniversalSet()
    elif pyeda_expr is Zero:
        return EmptySet()
    elif isinstance(pyeda_expr, Variable):
        return ValueSet(sym_to_val[pyeda_expr.name])
    elif isinstance(pyeda_expr, AndOp):
        return Intersection(*(venn_from_pyeda(x, sym_to_val) for x in pyeda_expr.xs))
    elif isinstance(pyeda_expr, OrOp):
        return Union(*(venn_from_pyeda(x, sym_to_val) for x in pyeda_expr.xs))
    elif isinstance(pyeda_expr, NotOp):
        return Complement(venn_from_pyeda(pyeda_expr.x, sym_to_val))
    elif isinstance(pyeda_expr, pyedaComplement):
        return Complement(ValueSet(sym_to_val[pyeda_expr.inputs[0].name]))
    elif isinstance(pyeda_expr, XorOp):
        return DisjunctiveUnion(*(venn_from_pyeda(x, sym_to_val) for x in pyeda_expr.xs))
    else:
        raise Exception


def venn_from_str(venn_str: str, value_parser: Callable[[str], T]) -> Set[T]:
    # The values have to be surrounded by a single space.
    BOOL_REGEX            = '[\(\)\|\&\~]+'
    pyeda_str             = ''
    pyeda_sym_to_val      = OrderedDict()  # type: MutableMapping[str, T]
    venn_sym_to_pyeda_sym = OrderedDict()  # type: MutableMapping[str, str]
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
