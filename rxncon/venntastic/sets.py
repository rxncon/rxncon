from functools import reduce
from itertools import combinations
from typing import List, Dict, Callable, Optional


METHOD_COMPLEMENTS_EXPANDED = '_complements_expanded'
METHOD_UNIONS_MOVED_TO_LEFT = '_unions_moved_to_left'
METHOD_UNIONS_ORDERED_SIMPLIFIED = '_unions_ordered_simplified'


class Set:
    def __hash__(self) -> int:
        return hash(str(self))

    def canonical_form(self) -> 'Set':
        simplification_methods = [
            METHOD_COMPLEMENTS_EXPANDED,
            METHOD_UNIONS_MOVED_TO_LEFT,
            METHOD_UNIONS_ORDERED_SIMPLIFIED
        ]

        return _call_method_until_stable(self, simplification_methods)

    @property
    def cardinality(self) -> Dict['PropertySet', int]:
        union_sets = self.canonical_form()._unions_flattened()
        terms = {}

        for i in range(len(union_sets)):
            tuples = combinations(union_sets, i + 1)
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

    def is_superset_of(self, other: 'Set') -> bool:
        return False

    def _complements_expanded(self) -> 'Set':
        return self

    def _unions_flattened(self) -> List['Set']:
        """
        Returns a flattened list given a left-simplified nested Union expression, for example:
            Union(Union(Union(a,b),c),d) -> [a, b, c, d]
        """
        return []

    def _unions_moved_to_left(self) -> 'Set':
        return self

    def _nodes_ordered(self, set_type, visited_sets: List['Set'],
                       postprocessing_func: Callable[[List['Set']], 'Set']) -> 'Set':
        """
        Recursively orders the nodes in a nested expression. If set_type is Union, the expression is of the form
        Union(Union(...)), if set_type is Intersection it is of the form Intersection(Intersection(...)).
        This works by first unpacking the nested expression into a list: elements to visited_sets are recursively added.

        The postprocessing_func maps the list of Sets to resultant Set. This allows us to use this method both in the
        case when we want to collapse all sub-expressions and in the case when we want to expand all sub-expressions.
        """
        return self

    def _unions_ordered_simplified(self, visited_sets: Optional[List['Set']]=None) -> 'Set':
        return self

    def _partial_cardinality(self) -> Dict['PropertySet', int]:
        raise AssertionError


class UnarySet(Set):
    def _unions_flattened(self) -> List[Set]:
        return [self]


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

    def _nodes_ordered(self, set_type, visited_sets: List[Set], postprocessing_func: Callable[[List[Set]], Set]) -> Set:
        if isinstance(self.left_expr, set_type):
            visited_sets.append(self.right_expr._nodes_ordered(set_type, [], postprocessing_func))
            return self.left_expr._nodes_ordered(set_type, visited_sets, postprocessing_func)

        else:
            visited_sets.append(self.right_expr)
            visited_sets.append(self.left_expr)
            visited_sets.reverse()

            processed = postprocessing_func(visited_sets)
            visited_sets.clear()

            return processed

    def _complements_expanded(self) -> Set:
        return type(self)(self.left_expr._complements_expanded(), self.right_expr._complements_expanded())


class PropertySet(Set):
    def __init__(self, value):
        assert hasattr(value, '__hash__')
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
            return self.value < other.value

        elif isinstance(other, Complement):
            return True

        elif isinstance(other, EmptySet):
            return False

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
        return self == other

    def is_subset_of(self, other: Set):
        return self == other

    def _partial_cardinality(self) -> Dict['PropertySet', int]:
        return {self: 1}


class EmptySet(UnarySet):
    """
    The empty set, which contains no elements.
    """
    def __eq__(self, other: Set) -> bool:
        return isinstance(other, EmptySet)

    def __lt__(self, other: Set) -> bool:
        if isinstance(other, PropertySet) or isinstance(other, Complement):
            return True

        else:
            return False

    def __repr__(self) -> str:
        return 'EmptySet'

    def __hash__(self) -> int:
        return hash(str(self))

    def is_superset_of(self, other: Set) -> bool:
        return isinstance(other, EmptySet)

    def _partial_cardinality(self) -> Dict[PropertySet, int]:
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

    def __lt__(self, other: Set) -> bool:
        if isinstance(other, Complement):
            return self.expr < other.expr

        else:
            return False

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return 'Complement({})'.format(self.expr)

    def is_superset_of(self, other: Set) -> bool:
        assert isinstance(other, Set)

        if isinstance(other, Complement) and other.expr.is_superset_of(self.expr):
            return True
        else:
            return False

    def is_subset_of(self, other: Set):
        return self == other

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

    def _partial_cardinality(self) -> Dict[PropertySet, int]:
        if isinstance(self.expr, PropertySet):
            return {UniversalSet(): 1, self.expr: -1}

        else:
            raise AssertionError


class Intersection(BinarySet):
    def __lt__(self, other: Set) -> bool:
        if isinstance(other, EmptySet) or isinstance(other, PropertySet) or isinstance(other, Complement):
            return True

        elif isinstance(other, Intersection):
            if self.left_expr < other.left_expr:
                return True

            elif not self.left_expr < other.left_expr and not other.left_expr < self.left_expr:
                return self.right_expr < other.right_expr

            else:
                return False

        else:
            return False

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return 'Intersection({0}, {1})'.format(self.left_expr, self.right_expr)

    def _unions_flattened(self) -> List[Set]:
        return [self]

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
    def __lt__(self, other: Set) -> bool:
        if isinstance(other, EmptySet) or isinstance(other, PropertySet) or \
                isinstance(other, Complement) or isinstance(other, Intersection):
            return True

        elif isinstance(other, Union):
            if self.left_expr < other.left_expr:
                return True

            elif not self.left_expr < other.left_expr and not other.left_expr < self.left_expr:
                return self.right_expr < other.right_expr

            else:
                return False

        else:
            raise AssertionError

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return 'Union({0}, {1})'.format(self.left_expr, self.right_expr)

    def _unions_flattened(self) -> List[Set]:
        if isinstance(self.right_expr, Union):
            raise AssertionError

        if isinstance(self.left_expr, Union):
            left_contribution = self.left_expr._unions_flattened()

        else:
            left_contribution = [self.left_expr]

        return left_contribution + [self.right_expr]

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

    def _unions_ordered_simplified(self, visited_sets: Optional[List[Set]]=None) -> Set:
        if visited_sets is None:
            visited_sets = []

        return self._nodes_ordered(Union, visited_sets, lambda x: _find_total_union(x))


class Difference(Set):
    """
    The difference of two Set expressions, by the identity A \ B = A ^ !B
    """
    def __new__(cls, *args, **kwargs) -> Set:
        assert len(args) == 2
        return Intersection(args[0], Complement(args[1]))

def UniversalSet():
    return PropertySet(None)


def flat_list_to_nested_expression(xs: List[Set], set_class) -> Set:
    """
    Reduces a list of Sets to a nested expression. E.g. if set_class is Intersection:
        [a, b, c] -> Intersection(Intersection(a, b), c)
    """
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
        return reduce(set_class, xs[1:], xs[0])


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
    """
    Flips the sign of the values of a dictionary.
    """
    res = {}
    for k, v in x.items():
        res[k] = -1 * v

    return res


def _has_mutual_complements(sets: List[Set], set_operation):
    """
    Returns whether a list of Sets, combined by a certain set_operation (Intersection or Union), contains two
    mutually complementary Sets.
    Used to simplify A ^ !A = 0 and A u !A = U.
    """
    if not (set_operation.__name__ == Union.__name__ or set_operation.__name__ == Intersection.__name__):
        raise AssertionError('_find_mutual_complements must have Union or Intersection as set_operation')

    expansion_methods = [
        METHOD_COMPLEMENTS_EXPANDED,
        METHOD_UNIONS_MOVED_TO_LEFT,
        METHOD_UNIONS_ORDERED_SIMPLIFIED
    ]

    for n in range(len(sets)):
        left_combis = list(combinations(sets, n + 1))

        for left_combi in left_combis:
            rest = [x for x in sets if x not in left_combi]

            for m in range(len(rest)):
                right_combis = list(combinations(rest, m + 1))

                for right_combi in right_combis:
                    lhs = _call_method_until_stable(
                        Complement(flat_list_to_nested_expression(left_combi, set_operation)),
                        expansion_methods
                    )

                    rhs = _call_method_until_stable(
                        flat_list_to_nested_expression(right_combi, set_operation),
                        expansion_methods
                    )

                    if lhs == rhs:
                        return True

    return False


def _find_total_intersection(sets: List[Set]) -> Set:
    """
    Returns the total intersection of a list of Sets.
    @todo We can simplify this even further by looking at subsets of sub-expressions:
          Intersection(a, b) = b if a.contains(b)
    """
    if EmptySet() in sets:
        return EmptySet()

    # Filter out the universal set and duplicates
    sets = list(set(sets))
    if sets == [UniversalSet()]:
        return UniversalSet()

    sets = [x for x in sets if x != UniversalSet()]

    if _has_mutual_complements(sets, Intersection):
        return EmptySet()

    property_sets = [x for x in sets if isinstance(x, PropertySet)]
    complements = [x for x in sets if isinstance(x, Complement)]

    properties = []

    for s in property_sets:
        properties += s.properties

    properties = list(set(properties))

    if len(complements) > 0 and len(properties) == 0:
        return flat_list_to_nested_expression(sorted(complements), Intersection)

    else:
        return flat_list_to_nested_expression(sorted(complements + [PropertySet(*properties)]), Intersection)


def _find_total_union(sets: List[Set]):
    """
    Returns the total union of a list of Set expressions.
    @todo We can simplify this even further by looking at subsets of sub-expressions:
          Union(a, b) = a if a.contains(b)
    """
    if UniversalSet() in sets:
        return UniversalSet()

    # Filter out the Empty set and duplicates
    sets = list(set(sets))
    if sets == [EmptySet()]:
        return EmptySet()

    sets = [x for x in sets if x != EmptySet()]

    if _has_mutual_complements(sets, Union):
        return UniversalSet()

    return flat_list_to_nested_expression(sorted(sets), Union)


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
