import functools
import itertools as itt
from typing import List, Dict, Optional


METHOD_COMPLEMENTS_EXPANDED = '_complements_expanded'
METHOD_UNIONS_MOVED_TO_LEFT = '_unions_moved_to_left'


class Set:
    def cardinality_form(self) -> 'Set':
        simplification_methods = [
            METHOD_COMPLEMENTS_EXPANDED,
            METHOD_UNIONS_MOVED_TO_LEFT
        ]

        return _call_method_list_until_stable(self, simplification_methods)

    def nested_list_form(self):
        return self.cardinality_form()._to_nested_list()

    @property
    def cardinality(self) -> Dict['PropertySet', int]:
        union_sets = self.cardinality_form()._unions_flattened()
        terms = {}

        for i in range(len(union_sets)):
            tuples = itt.combinations(union_sets, i + 1)
            tuple_intersections = []

            for t in tuples:
                tuple_intersections.append(nested_expression_from_list_and_binary_op(t, Intersection))

            tuple_intersections = [x.cardinality_form() for x in tuple_intersections]

            for t in tuple_intersections:
                if i % 2 == 0:
                    terms = _add_dicts(terms, t._partial_cardinality())
                else:
                    terms = _add_dicts(terms, _negate_dict(t._partial_cardinality()))

        return {k: v for k, v in terms.items() if v != 0}

    def is_equivalent_to(self, other: 'Set'):
        assert isinstance(other, Set)

        return self.cardinality_form() == other.cardinality_form()

    def is_superset_of(self, other: 'Set') -> Optional[bool]:
        return None

    def is_subset_of(self, other: 'Set') -> Optional[bool]:
        return None

    def _complements_expanded(self) -> 'Set':
        return self

    def _unions_moved_to_left(self) -> 'Set':
        return self

    def _to_nested_list(self):
        pass

    def _partial_cardinality(self) -> Dict['Set', int]:
        raise NotImplementedError


class UnarySet(Set):
    def _to_nested_list(self):
        return [[self]]


class PropertySet(UnarySet):
    def __init__(self, value):
        assert hash(value)
        self.value = value

    def __eq__(self, other: Set) -> bool:
        assert isinstance(other, Set)

        if isinstance(other, PropertySet):
            return self.value == other.value

        else:
            return False

    def __hash__(self):
        return hash('*property-set-{}*'.format(hash(self.value)))

    def __lt__(self, other: Set):
        if isinstance(other, PropertySet):
            if self.value is None:
                return True

            return self.value < other.value

        elif isinstance(other, EmptySet):
            return False

        elif isinstance(other, Complement) or isinstance(other, BinarySet):
            return True

        else:
            raise AssertionError

    def __repr__(self):
        return str(self)

    def __str__(self):
        if self.value:
            return 'Property({})'.format(self.value)

        else:
            return 'UniversalSet'

    def is_superset_of(self, other: Set):
        if self.value is None:
            return True

        elif isinstance(other, EmptySet):
            return True

        elif isinstance(other, Intersection):
            return self.is_superset_of(other.left_expr) or self.is_superset_of(other.right_expr)

        else:
            return self == other

    def is_subset_of(self, other: Set):
        if self == other:
            return True

        else:
            return other.is_superset_of(self)

    def _partial_cardinality(self) -> Dict['Set', int]:
        return {self: 1}


class EmptySet(UnarySet):
    def __eq__(self, other: Set) -> bool:
        return isinstance(other, EmptySet)

    def __hash__(self) -> int:
        return hash('*empty-set*')

    def __lt__(self, other: Set) -> bool:
        assert isinstance(other, Set)
        return False if isinstance(other, EmptySet) else True

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'EmptySet'

    def is_superset_of(self, other: Set) -> Optional[bool]:
        return isinstance(other, EmptySet)

    def is_subset_of(self, other: Set) -> Optional[bool]:
        return True

    def _partial_cardinality(self) -> Dict['Set', int]:
        return {}


class Complement(UnarySet):
    def __init__(self, expr: Set):
        self.expr = expr

    def __eq__(self, other: Set) -> bool:
        assert isinstance(other, Set)

        if isinstance(other, Complement):
            return self.expr == other.expr

        else:
            return False

    def __hash__(self) -> int:
        return hash('*complement-{}*'.format(self.expr))

    def __lt__(self, other: Set) -> bool:
        if isinstance(other, Complement):
            return self.expr < other.expr

        elif isinstance(other, EmptySet) or isinstance(other, PropertySet):
            return False

        else:
            return True

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'Complement({})'.format(self.expr)

    def is_superset_of(self, other: Set):
        if self == other:
            return True

        elif isinstance(other, Complement):
            return self.expr.is_subset_of(other.expr)

        else:
            return None

    def is_subset_of(self, other: Set):
        if self == other:
            return True

        elif self.is_superset_of(other) is None:
            return None

        else:
            return not self.is_superset_of(other)

    def _complements_expanded(self) -> bool:
        if isinstance(self.expr, Complement):
            return self.expr.expr._complements_expanded()

        elif isinstance(self.expr, EmptySet):
            return UniversalSet()

        elif self.expr == UniversalSet():
            return EmptySet()

        elif isinstance(self.expr, Union):
            return Intersection(Complement(self.expr.left_expr), Complement(self.expr.right_expr))._complements_expanded()

        elif isinstance(self.expr, Intersection):
            return Union(Complement(self.expr.left_expr), Complement(self.expr.right_expr))._complements_expanded()

        else:
            return Complement(self.expr._complements_expanded())

    def _partial_cardinality(self) -> Dict['Set', int]:
        if isinstance(self.expr, PropertySet):
            return {UniversalSet(): 1, self.expr: -1}

        else:
            raise AssertionError


class BinarySet(Set):
    def __init__(self, left_expr: Set, right_expr: Set):
        self.left_expr = left_expr
        self.right_expr = right_expr

    def __eq__(self, other: Set) -> bool:
        assert isinstance(other, Set)

        if isinstance(other, type(self)) and (self.left_expr == other.left_expr) and (self.right_expr == other.right_expr):
            return True

        else:
            return False

    def _complements_expanded(self) -> Set:
        return type(self)(self.left_expr._complements_expanded(), self.right_expr._complements_expanded())


class Intersection(BinarySet):
    def __lt__(self, other: Set) -> bool:
        if isinstance(other, EmptySet) or isinstance(other, PropertySet) or isinstance(other, Complement):
            return False

        elif isinstance(other, Intersection):
            if self.left_expr < other.left_expr:
                return True

            elif self.left_expr == other.left_expr:
                return self.right_expr < other.right_expr

            return False

        elif isinstance(other, Union):
            return True

        else:
            raise AssertionError

    def __hash__(self):
        return hash('*intersection-{0}{1}*'.format(hash(self.left_expr), hash(self.right_expr)))

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'Intersection({0}, {1})'.format(self.left_expr, self.right_expr)

    def is_superset_of(self, other: Set) -> Optional[bool]:
        if self == other:
            return True

        else:
            return other.is_subset_of(self)

    def is_subset_of(self, other: Set) -> Optional[bool]:
        if self == other:
            return True

        elif isinstance(other, Intersection):
            return self.left_expr.is_subset_of(other.left_expr) or self.right_expr.is_subset_of(other.right_expr) or\
                   self.right_expr.is_subset_of(other.left_expr) or self.left_expr.is_subset_of(other.right_expr)

        elif isinstance(other, PropertySet):
            return other.is_superset_of(self)

        else:
            return None

    def _unions_moved_to_left(self) -> Set:
        if isinstance(self.left_expr, Union):
            # Distributivity of intersection with respect to union
            return Union(Intersection(self.left_expr.left_expr, self.right_expr), Intersection(self.left_expr.right_expr, self.right_expr))\
                ._unions_moved_to_left()

        elif isinstance(self.right_expr, Union):
            # Distributivity of intersection with respect to union
            return Union(Intersection(self.left_expr, self.right_expr.left_expr), Intersection(self.left_expr, self.right_expr.right_expr))\
                ._unions_moved_to_left()

        else:
            return Intersection(self.left_expr._unions_moved_to_left(), self.right_expr._unions_moved_to_left())

    def _to_nested_list(self):
        return [self.left_expr._to_nested_list()[0] + self.right_expr._to_nested_list()[0]]

    def _partial_cardinality(self) -> Dict[PropertySet, int]:
        if isinstance(self.left_expr, Complement):
            # Remove single complement using
            # |!A ^ B| = |B| - |A ^ B|
            return _add_dicts(self.right_expr._partial_cardinality(),
                              _negate_dict(Intersection(self.left_expr.expr, self.right_expr)
                                           .cardinality_form()
                                           ._partial_cardinality()))

        elif isinstance(self.right_expr, Complement):
            # Remove single complement using
            # |A ^ !B| = |A| - |B ^ A|
            return _add_dicts(self.left_expr._partial_cardinality(),
                              _negate_dict(Intersection(self.right_expr.expr, self.left_expr)
                                           .cardinality_form()
                                           ._partial_cardinality()))

        else:
            raise AssertionError


class Union(BinarySet):
    def __hash__(self):
        return hash('*union-{0}{1}*'.format(hash(self.left_expr), hash(self.right_expr)))

    def __lt__(self, other: Set) -> bool:
        if isinstance(other, EmptySet) or isinstance(other, PropertySet) or \
                isinstance(other, Complement) or isinstance(other, Intersection):
            return False

        elif isinstance(other, Union):
            if self.left_expr < other.left_expr:
                return True

            elif self.left_expr == other.left_expr:
                return self.right_expr < other.right_expr

            return False

        else:
            raise AssertionError

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'Union({0}, {1})'.format(self.left_expr, self.right_expr)

    def is_superset_of(self, other: Set) -> Optional[bool]:
        if self == other:
            return True

        elif isinstance(other, Union):
            return self.left_expr.is_superset_of(other.left_expr) or self.right_expr.is_superset_of(other.right_expr) or\
                   self.right_expr.is_superset_of(other.left_expr) or self.left_expr.is_superset_of(other.right_expr)

        elif isinstance(other, PropertySet):
            return self.left_expr.is_superset_of(other) or self.right_expr.is_superset_of(other)

        else:
            return None

    def is_subset_of(self, other: Set) -> Optional[bool]:
        if self == other:
            return True

        else:
            return other.is_superset_of(self)

    def _unions_moved_to_left(self) -> Set:
        if isinstance(self.right_expr, Union) and not isinstance(self.left_expr, Union):
            # Use commutativity to move all Union expressions to the left
            return Union(self.right_expr, self.left_expr)._unions_moved_to_left()

        elif isinstance(self.left_expr, Union) and isinstance(self.right_expr, Union):
            # Use associativity to move all Union expressions to the left
            return Union(Union(Union(self.left_expr.left_expr, self.left_expr.right_expr), self.right_expr.left_expr),
                         self.right_expr.right_expr)._unions_moved_to_left()

        else:
            return Union(self.left_expr._unions_moved_to_left(), self.right_expr._unions_moved_to_left())

    def _to_nested_list(self):
        return self.left_expr._to_nested_list() + self.right_expr._to_nested_list()

    def _partial_cardinality(self) -> Dict[PropertySet, int]:
        raise AssertionError


class Difference(Set):
    def __new__(cls, *args, **kwargs) -> Set:
        assert len(args) == 2
        return Intersection(args[0], Complement(args[1]))


def UniversalSet():
    return PropertySet(None)


def nested_expression_from_list_and_binary_op(xs: List[Set], binary_op) -> Set:
    if binary_op == Intersection:
        unit = UniversalSet()

    elif binary_op == Union:
        unit = EmptySet()

    else:
        raise AssertionError

    if len(xs) == 0:
        return unit

    elif len(xs) == 1:
        return xs[0]

    else:
        return functools.reduce(binary_op, xs[1:], xs[0])


def boolean_function_from_nested_list_form(nested_list):
    set_to_argument_num = {}

    argument_num = 0
    and_clauses = []

    for term in nested_list:
        required_true = []
        required_false = []

        for item in term:
            if isinstance(item, PropertySet) and not item in set_to_argument_num.keys():
                set_to_argument_num[item] = argument_num
                argument_num += 1

            elif isinstance(item, Complement) and not item.expr in set_to_argument_num.keys():
                set_to_argument_num[item.expr] = argument_num
                argument_num += 1


            if isinstance(item, PropertySet):
                required_true.append(set_to_argument_num[item])

            elif isinstance(item, Complement):
                required_false.append(set_to_argument_num[item.expr])

            else:
                raise AssertionError

        and_clauses.append(BooleanAndClause(required_true, required_false))

    return BooleanFunction(and_clauses)


class BooleanFunction:
    def __init__(self, and_clauses):
        assert all([isinstance(x, BooleanAndClause) for x in and_clauses])
        self.and_clauses = and_clauses
        self.valid_statements = list(set([statement for and_clause in self.and_clauses for statement in and_clause.valid_statements]))

    def __call__(self, *args, **kwargs):
        actual_true = kwargs.get('true', [])
        actual_false = kwargs.get('false', [])

        if any([statements not in self.valid_statements for statements in actual_true + actual_false]):
            raise NameError

        return any([and_clause(**{'true': actual_true, 'false': actual_false} ) for and_clause in self.and_clauses])


class BooleanAndClause:
    def __init__(self, required_true=None, required_false=None):
        if required_true is None:
            required_true = []

        if required_false is None:
            required_false = []

        assert not any([x in required_false for x in required_true])
        assert not any([x in required_true for x in required_false])

        self.required_true = required_true
        self.required_false = required_false
        self.valid_statements = list(set(required_true + required_false))

    def __call__(self, *args, **kwargs):
        if any([x in self.required_false for x in kwargs['true']]) or any([x in self.required_true for x in kwargs['false']]):
            return False

        else:
            return True


def _add_dicts(x: Dict[PropertySet, int], y: Dict[PropertySet, int]) -> Dict[PropertySet, int]:
    res = {}
    for k, v in x.items():
        res[k] = v

    for k, v in y.items():
        if k in res:
            res[k] += v

        else:
            res[k] = v

    return res


def _negate_dict(x: Dict[PropertySet, int]) -> Dict[PropertySet, int]:
    res = {}
    for k, v in x.items():
        res[k] = -1 * v

    return res


def _call_method_list_until_stable(expr: Set, methods: List[str]):
    max_simplifications = 100
    i = 0

    previous_simplification = expr
    current_simplification = expr
    simplification_done = False

    while not simplification_done:
        i += 1
        if i > max_simplifications:
            raise RecursionError

        for method in methods:
            current_simplification = getattr(current_simplification, method)()

        if current_simplification == previous_simplification:
            simplification_done = True

        else:
            previous_simplification = current_simplification

    return current_simplification











