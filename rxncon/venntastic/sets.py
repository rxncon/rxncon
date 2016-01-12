import functools
import itertools as itt
from typing import List, Dict, Optional


METHOD_COMPLEMENTS_EXPANDED = '_complements_expanded'
METHOD_UNIONS_MOVED_TO_LEFT = '_unions_moved_to_left'
METHOD_INTERSECTIONS_SIMPLIFIED = '_intersections_simplified'
METHOD_UNIONS_SIMPLIFIED = '_unions_simplified'

class Set:
    def canonical_form(self) -> 'Set':
        simplification_methods = [
            METHOD_COMPLEMENTS_EXPANDED,
            METHOD_UNIONS_MOVED_TO_LEFT
        ]

        return _call_method_list_until_stable(self, simplification_methods)

    @property
    def cardinality(self) -> Dict['PropertySet', int]:
        union_sets = self.canonical_form()._unions_flattened()
        terms = {}

        for i in range(len(union_sets)):
            tuples = itt.combinations(union_sets, i + 1)
            tuple_intersections = []

            for t in tuples:
                tuple_intersections.append(flat_list_to_nested_expression(t, Intersection))

            tuple_intersections = [x.canonical_form() for x in tuple_intersections]

            for t in tuple_intersections:
                if i % 2 == 0:
                    terms = _add_dicts(terms, t._partial_cardinality())
                else:
                    terms = _add_dicts(terms, _negate_dict(t._partial_cardinality()))

        return {k: v for k, v in terms.items() if v != 0}

    def equivalent_forms(self):
        return [self]

    def is_equivalent_to(self, other: 'Set'):
        assert isinstance(other, Set)

        return self in other.equivalent_forms() or other in self.equivalent_forms()

    def is_superset_of(self, other: 'Set') -> Optional[bool]:
        return None

    def is_subset_of(self, other: 'Set') -> Optional[bool]:
        return None

    def _complements_expanded(self) -> 'Set':
        return self

    def _unions_moved_to_left(self) -> 'Set':
        return self

    def _partial_cardinality(self) -> Dict['Set', int]:
        raise NotImplementedError


class UnarySet(Set):
    def _unions_flattened(self) -> List[Set]:
        return [self]


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

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'Complement({})'.format(self.expr)

    def equivalent_forms(self):
        forms = [self]

        if isinstance(self.expr, Complement):
            forms += self.expr.expr.equivalent_forms()

        elif isinstance(self.expr, Union):
            forms += Intersection(Complement(self.expr.left_expr), Complement(self.expr.right_expr)).equivalent_forms()

        elif isinstance(self.expr, Intersection):
            forms += Union(Complement(self.expr.left_expr), Complement(self.expr.right_expr)).equivalent_forms()

        elif isinstance(self.expr, EmptySet):
            forms += [UniversalSet()]

        elif self.expr == UniversalSet():
            forms += [EmptySet()]

        return forms

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

    def flip(self):
        return type(self)(self.right_expr, self.left_expr)

    def equivalent_forms(self):
        X = type(self)

        if X == Union:
            Y = Intersection
            U = UniversalSet()
            S = EmptySet()

        elif X == Intersection:
            Y = Union
            U = EmptySet()
            S = UniversalSet()

        forms = []

        for left, right in itt.product(self.left_expr.equivalent_forms(), self.right_expr.equivalent_forms()):
            forms += [X(left, right),
                      X(left, right).flip()]

            if left in Complement(right).equivalent_forms() or right in Complement(left).equivalent_forms():
                forms += [U]

            if left == U or right == U:
                forms += [U]

            if left == S:
                forms += [right]

            if right == S:
                forms += [left]

            if left == right:
                forms += [left]

            if isinstance(left, X):
                forms += [X(left.left_expr, X(left.right_expr, right)),
                          X(left.left_expr, X(left.right_expr, right)).flip(),
                          X(left.left_expr, X(left.right_expr, right).flip()),
                          X(left.left_expr, X(left.right_expr, right).flip()).flip()]

            if isinstance(right, X):
                forms += [X(X(left, right.left_expr), right.right_expr),
                          X(X(left, right.left_expr), right.right_expr).flip(),
                          X(X(left, right.left_expr).flip(), right.right_expr),
                          X(X(left, right.left_expr).flip(), right.right_expr).flip()]

            if isinstance(left, Y):
                forms += [Y(X(left.left_expr, right), X(left.right_expr, right)),
                          Y(X(left.left_expr, right).flip(), X(left.right_expr, right)),
                          Y(X(left.left_expr, right), X(left.right_expr, right).flip()),
                          Y(X(left.left_expr, right).flip(), X(left.right_expr, right).flip()),
                          Y(X(left.left_expr, right), X(left.right_expr, right)).flip(),
                          Y(X(left.left_expr, right).flip(), X(left.right_expr, right)).flip(),
                          Y(X(left.left_expr, right), X(left.right_expr, right).flip()).flip(),
                          Y(X(left.left_expr, right).flip(), X(left.right_expr, right).flip()).flip()]

            if isinstance(right, Y):
                forms += [Y(X(left, right.left_expr), X(left, right.right_expr)),
                          Y(X(left, right.left_expr).flip(), X(left, right.right_expr)),
                          Y(X(left, right.left_expr), X(left, right.right_expr).flip()),
                          Y(X(left, right.left_expr).flip(), X(left, right.right_expr).flip()),
                          Y(X(left, right.left_expr), X(left, right.right_expr)).flip(),
                          Y(X(left, right.left_expr).flip(), X(left, right.right_expr)).flip(),
                          Y(X(left, right.left_expr), X(left, right.right_expr).flip()).flip(),
                          Y(X(left, right.left_expr).flip(), X(left, right.right_expr).flip()).flip()]

        return list(set(forms))


    def _complements_expanded(self) -> Set:
        return type(self)(self.left_expr._complements_expanded(), self.right_expr._complements_expanded())


class Intersection(BinarySet):
    def __hash__(self):
        return hash('*intersection-{0}{1}*'.format(hash(self.left_expr), hash(self.right_expr)))

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'Intersection({0}, {1})'.format(self.left_expr, self.right_expr)

    def is_superset_of(self, other: Set) -> Optional[bool]:
        if self == other:
            return True

        elif self.is_subset_of(other) is None:
            return None

        else:
            return not self.is_subset_of(other)

    def is_subset_of(self, other: Set) -> Optional[bool]:
        if self == other:
            return True

        elif self.left_expr.is_subset_of(other) or self.right_expr.is_subset_of(other):
            return True

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

    def _partial_cardinality(self) -> Dict[PropertySet, int]:
        if isinstance(self.left_expr, Complement):
            # Remove single complement using
            # |!A ^ B| = |B| - |A ^ B|
            return _add_dicts(self.right_expr._partial_cardinality(),
                              _negate_dict(Intersection(self.left_expr.expr, self.right_expr)
                                           .canonical_form()
                                           ._partial_cardinality()))

        elif isinstance(self.right_expr, Complement):
            # Remove single complement using
            # |A ^ !B| = |A| - |B ^ A|
            return _add_dicts(self.left_expr._partial_cardinality(),
                              _negate_dict(Intersection(self.right_expr.expr, self.left_expr)
                                           .canonical_form()
                                           ._partial_cardinality()))

        else:
            raise AssertionError


class Union(BinarySet):
    def __hash__(self):
        return hash('*union-{0}{1}*'.format(hash(self.left_expr), hash(self.right_expr)))

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'Union({0}, {1})'.format(self.left_expr, self.right_expr)

    def is_superset_of(self, other: Set):
        if self == other:
            return True

        elif self.left_expr.is_superset_of(other) or self.right_expr.is_superset_of(other):
            return True

        else:
            return None

    def is_subset_of(self, other: Set):
        if self == other:
            return True

        elif self.is_superset_of(other) is None:
            return None

        else:
            return False

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

    def _partial_cardinality(self) -> Dict[PropertySet, int]:
        raise AssertionError


class Difference(Set):
    def __new__(cls, *args, **kwargs) -> Set:
        assert len(args) == 2
        return Intersection(args[0], Complement(args[1]))


def UniversalSet():
    return PropertySet(None)


def flat_list_to_nested_expression(xs: List[Set], set_class) -> Set:
    if set_class == Intersection:
        unit = UniversalSet()

    elif set_class == Union:
        unit = EmptySet()

    else:
        raise AssertionError

    if len(xs) == 0:
        return unit

    elif len(xs) == 1:
        return xs[0]

    else:
        return functools.reduce(set_class, xs[1:], xs[0])


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


def _call_method_until_stable(expr: Set, method):
    max_simplifications = 100
    i = 0

    previous_simplification = expr
    current_simplification = expr
    simplification_done = False

    while not simplification_done:
        i += 1
        if i > max_simplifications:
            raise RecursionError

        current_simplification = getattr(current_simplification, method)()

        if current_simplification == previous_simplification:
            simplification_done = True

        else:
            previous_simplification = current_simplification

    return current_simplification


def _call_method_list_until_stable(expr: Set, methods: List[str]):
    previous_simplification = expr
    current_simplification = expr
    simplification_done = False

    while not simplification_done:
        for method in methods:
            current_simplification = _call_method_until_stable(current_simplification, method)

        if current_simplification == previous_simplification:
            simplification_done = True

        else:
            previous_simplification = current_simplification

    return current_simplification
