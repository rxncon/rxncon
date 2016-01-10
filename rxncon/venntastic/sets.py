from functools import reduce
from itertools import combinations
from typing import List, Dict, Callable, Iterable, Optional


class Set:
    """
    Parent class for the different types of Set expressions.
    """
    def __hash__(self) -> int:
        return hash(str(self))

    def simplified(self) -> 'Set':
        """
        Returns a simplified version of the Set, in which:
           * Complement(x) only appears for x a PropertySet
           * all Unions are moved to the left by the commutativity and associativity axioms
           * within a Union node, all Intersections are moved to the left by the commutativity and associativity axioms
           * all Intersections are collapsed into simpler expressions (PropertySets are combined and A ^ !A = 0)
           * all Unions are collapsed into simpler expressions (A u !A = U)
        """
        simplification_methods = [
            '_expand_complements',
            '_move_unions_to_left',
            '_move_intersections_to_left',
            '_order_collapse_intersections',
            '_order_collapse_unions'
        ]

        return _apply_until_stable(self, simplification_methods)

    def cardinality(self) -> Dict['PropertySet', int]:
        """
        Returns the cardinality of a Set expression as the sum of intersections of PropertySets. This uses the
        inclusion-exclusion principle. The return value is a dictionary in which the keys are the PropertySets and
        the values are the coefficient in which they appear in the sum.
        """
        union_sets = self.simplified()._flatten_unions()
        terms = {}

        for i in range(len(union_sets)):
            tuples = combinations(union_sets, i + 1)
            tuple_intersections = []

            for t in tuples:
                tuple_intersections.append(flat_list_to_nested_expression(t, Intersection))

            tuple_intersections = [x.simplified() for x in tuple_intersections]

            for t in tuple_intersections:
                if i % 2 == 0:
                    terms = _add_dicts(terms, t._partial_cardinality())
                else:
                    terms = _add_dicts(terms, _negate_dict(t._partial_cardinality()))

        return {k: v for k, v in terms.items() if v != 0}

    def is_superset_of(self, other: 'Set') -> bool:
        """
        Returns whether this set contains the other set, i.e. whether the other set is a subset of self.
        """
        return False

    def _expand_complements(self) -> 'Set':
        """
        Returns a Set in which the Complement expressions are expanded, using among others De Morgan's identity.
        """
        return self

    def _flatten_unions(self) -> List['Set']:
        """
        Returns a flattened list given a left-simplified nested Union expression, for example:
            Union(Union(Union(a,b),c),d) -> [a, b, c, d]
        """
        return []

    def _move_unions_to_left(self) -> 'Set':
        """
        Moves all Unions to the left, by using the mutual distributivity of Union and Intersection, the
        commutativity and the associativity of Unions
        """
        return self

    def _move_intersections_to_left(self) -> 'Set':
        """
        Moves all Intersections to the left, by using the commutativity and the associativity of Intersections
        """
        return self

    def _order_nodes(self, set_type, visited_sets: List['Set'],
                     postprocessing_func: Callable[[List['Set']], 'Set']) -> 'Set':
        """
        Recursively orders the nodes in a nested expression. If set_type is Union, the expression is of the form
        Union(Union(...)), if set_type is Intersection it is of the form Intersection(Intersection(...)).
        This works by first unpacking the nested expression into a list: elements to visited_sets are recursively added.

        The postprocessing_func maps the list of Sets to resultant Set. This allows us to use this method both in the
        case when we want to collapse all sub-expressions and in the case when we want to expand all sub-expressions.
        """
        return self

    def _order_expand_intersections(self, visited_sets: Optional[List['Set']]=None) -> 'Set':
        """
        Orders a nested expression and expands all the contained Intersections. Calls _order_nodes.
        """
        return self

    def _order_collapse_intersections(self, visited_sets: Optional[List['Set']]=None) -> 'Set':
        """
        Orders a nested expression and collapses all the contained Intersections. Calls _order_nodes.
        """

        return self

    def _order_collapse_unions(self, visited_sets: Optional[List['Set']]=None) -> 'Set':
        """
        Orders a nested expression and collapses all the contained Unions. Calls _order_nodes.
        """
        return self

    def _partial_cardinality(self) -> Dict['PropertySet', int]:
        """
        Returns the cardinality of a sub-expression.
        """
        raise AssertionError('Set._partial_cardinality was called')


class UnarySet(Set):
    """
    Parent class of PropertySet and Complement.
    """
    def _flatten_unions(self) -> List[Set]:
        return [self]


class BinarySet(Set):
    """
    Parent class of Union and Intersection.
    """
    def __init__(self, left_expr: Set, right_expr: Set):
        self.left_expr = left_expr
        self.right_expr = right_expr

    def __eq__(self, other: Set) -> bool:
        assert isinstance(other, Set)

        if isinstance(other, type(self)) and (self.left_expr == other.left_expr) and \
                (self.right_expr == other.right_expr):
            return True

        else:
            return False

    def _order_nodes(self, set_type, visited_sets: List[Set], postprocessing_func: Callable[[List[Set]], Set]) -> Set:
        if isinstance(self.left_expr, set_type):
            visited_sets.append(self.right_expr._order_nodes(set_type, [], postprocessing_func))
            return self.left_expr._order_nodes(set_type, visited_sets, postprocessing_func)

        else:
            visited_sets.append(self.right_expr)
            visited_sets.append(self.left_expr)
            visited_sets.reverse()

            processed = postprocessing_func(visited_sets)
            visited_sets.clear()

            return processed

    def _expand_complements(self) -> Set:
        return self.__class__(
            self.left_expr._expand_complements(),
            self.right_expr._expand_complements()
        )

    def _order_expand_intersections(self, visited_sets: Optional[List[Set]]=None) -> Set:
        return self.__class__(
            self.left_expr._order_expand_intersections([]),
            self.right_expr._order_expand_intersections([])
        )

    def _order_collapse_intersections(self, visited_sets: Optional[List[Set]]=None) -> Set:
        return self.__class__(
            self.left_expr._order_collapse_intersections([]),
            self.right_expr._order_collapse_intersections([])
        )


class PropertySet(UnarySet):
    """
    Set parametrized by the properties it has. The properties are integers, and logically the following equality holds:
    PropertySet(1, 2, 3) = Intersection(Intersection(PropertySet(1), PropertySet(2)), PropertySet(3))

    PropertySet() is the universal set containing all the elements.
    """
    def __init__(self, *args):
        self.properties = _unique(args)
        self.properties.sort()

    def __eq__(self, other: Set) -> bool:
        if isinstance(other, PropertySet):
            return self.properties == other.properties

        else:
            return False

    def __lt__(self, other) -> bool:
        if isinstance(other, Complement):
            return True

        elif isinstance(other, PropertySet):
            return self.properties < other.properties

        else:
            return False

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return '(' + ','.join(map(str, self.properties)) + ')'

    def is_superset_of(self, other: Set) -> bool:
        if isinstance(other, EmptySet):
            return True

        if not self.properties:
            return True

        if isinstance(other, PropertySet):
            if len(self.properties) > len(other.properties):
                return False
            elif len(self.properties) == len(other.properties):
                return self.properties == other.properties
            elif len(self.properties) < len(other.properties):
                for state in self.properties:
                    if state not in other.properties:
                        return False
                return True
        else:
            return False

    def to_list(self) -> List['PropertySet']:
        return [PropertySet(x) for x in self.properties]

    def _order_expand_intersections(self, visited_sets: Optional[List[Set]]=None) -> Set:
        if len(self.properties) > 1:
            return flat_list_to_nested_expression([PropertySet(x) for x in self.properties], Intersection)

        else:
            return self

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
        return 'Empty'

    def __hash__(self) -> int:
        return hash(str(self))

    def is_superset_of(self, other: Set) -> bool:
        return isinstance(other, EmptySet)

    def _partial_cardinality(self) -> Dict[PropertySet, int]:
        return {}


class Complement(UnarySet):
    """
    The complement of a Set expression.
    """
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
        return 'Complement(' + str(self.expr) + ')'

    def is_superset_of(self, other: Set) -> bool:
        assert isinstance(other, Set)

        if isinstance(other, Complement) and other.expr.is_superset_of(self.expr):
            return True
        else:
            return False

    def _expand_complements(self) -> bool:
        if isinstance(self.expr, Complement):
            # Complement of complement is no-op
            return self.expr.expr._expand_complements()

        elif isinstance(self.expr, EmptySet):
            # Complement of empty set is universal set
            return PropertySet()

        elif isinstance(self.expr, PropertySet) and self.expr == PropertySet():
            # Complement of universal set is empty set
            return EmptySet()

        elif isinstance(self.expr, Union):
            # De Morgan's law
            return Intersection(Complement(self.expr.left_expr), Complement(self.expr.right_expr))._expand_complements()

        elif isinstance(self.expr, Intersection):
            # De Morgan's law
            return Union(Complement(self.expr.left_expr), Complement(self.expr.right_expr))._expand_complements()

        else:
            return Complement(self.expr._expand_complements())

    def _partial_cardinality(self) -> Dict[PropertySet, int]:
        if isinstance(self.expr, PropertySet):
            return {PropertySet(): 1, self.expr: -1}

        else:
            raise AssertionError('Complement._partial_cardinality was called with expr not Set')


class Intersection(BinarySet):
    """
    The Intersection of two Sets.
    """
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
        return 'Intersection(' + str(self.left_expr) + ', ' + str(self.right_expr) + ')'

    def _flatten_unions(self) -> List[Set]:
        return [self]

    def _move_unions_to_left(self) -> Set:
        if isinstance(self.left_expr, Union):
            # Distributivity of intersection with respect to union
            return Union(
                Intersection(self.left_expr.left_expr, self.right_expr),
                Intersection(self.left_expr.right_expr, self.right_expr)
            )._move_unions_to_left()

        elif isinstance(self.right_expr, Union):
            # Distributivity of intersection with respect to union
            return Union(
                Intersection(self.left_expr, self.right_expr.left_expr),
                Intersection(self.left_expr, self.right_expr.right_expr)
            )._move_unions_to_left()

        else:
            return Intersection(
                self.left_expr._move_unions_to_left(),
                self.right_expr._move_unions_to_left()
            )

    def _move_intersections_to_left(self) -> Set:
        if isinstance(self.right_expr, Intersection) and not isinstance(self.left_expr, Intersection):
            # Use commutativity to move the Intersection expressions to the left
            return Intersection(self.right_expr, self.left_expr)._move_intersections_to_left()

        elif isinstance(self.left_expr, Intersection) and isinstance(self.right_expr, Intersection):
            # Use associativity to move the Intersection expressions to the left
            return Intersection(Intersection(Intersection(
                self.left_expr.left_expr, self.left_expr.right_expr), self.right_expr.left_expr),
                self.right_expr.right_expr)._move_intersections_to_left()
        else:
            return Intersection(
                self.left_expr._move_intersections_to_left(),
                self.right_expr._move_intersections_to_left()
            )

    def _order_expand_intersections(self, visited_sets: Optional[List[Set]]=None) -> Set:
        if visited_sets is None:
            visited_sets = []

        return self._order_nodes(Intersection,
                                 visited_sets,
                                 lambda x: flat_list_to_nested_expression(sorted(x), Intersection))

    def _order_collapse_intersections(self, visited_sets: Optional[List[Set]]=None) -> Set:
        if visited_sets is None:
            visited_sets = []

        return self._order_nodes(Intersection,
                                 visited_sets,
                                 lambda x: _find_total_intersection(x))

    def _partial_cardinality(self) -> Dict[PropertySet, int]:
        if isinstance(self.left_expr, Complement):
            # Remove single complement using
            # |!A ^ B| = |B| - |A ^ B|
            return _add_dicts(self.right_expr._partial_cardinality(),
                              _negate_dict(Intersection(self.left_expr.expr, self.right_expr)
                                           .simplified()
                                           ._partial_cardinality()))

        elif isinstance(self.right_expr, Complement):
            # Remove single complement using
            # |A ^ !B| = |A| - |B ^ A|
            return _add_dicts(self.left_expr._partial_cardinality(),
                              _negate_dict(Intersection(self.right_expr.expr, self.left_expr)
                                           .simplified()
                                           ._partial_cardinality()))

        else:
            raise AssertionError('Intersection._partial_cardinality called with left_expr: ' + str(self.left_expr) +
                                 ' right_expr: ' + str(self.right_expr))


class Union(BinarySet):
    """
    The union of two Set expressions.
    """
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
            raise AssertionError('Union.__lt__ non-exhaustive match')

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return 'Union(' + str(self.left_expr) + ', ' + str(self.right_expr) + ')'

    def _flatten_unions(self) -> List[Set]:
        if isinstance(self.right_expr, Union):
            raise AssertionError('Union._flatten found union on right hand side of tree')

        if isinstance(self.left_expr, Union):
            left_contribution = self.left_expr._flatten_unions()

        else:
            left_contribution = [self.left_expr]

        return left_contribution + [self.right_expr]

    def _move_unions_to_left(self) -> Set:
        if isinstance(self.right_expr, Union) and not isinstance(self.left_expr, Union):
            # Use commutativity to move all Union expressions to the left
            return Union(self.right_expr, self.left_expr)._move_unions_to_left()

        elif isinstance(self.left_expr, Union) and isinstance(self.right_expr, Union):
            # Use associativity to move all Union expressions to the left
            return Union(Union(Union(self.left_expr.left_expr, self.left_expr.right_expr), self.right_expr.left_expr),
                         self.right_expr.right_expr)._move_unions_to_left()

        else:
            return Union(self.left_expr._move_unions_to_left(), self.right_expr._move_unions_to_left())

    def _move_intersections_to_left(self) -> Set:
        return Union(self.left_expr._move_intersections_to_left(), self.right_expr._move_intersections_to_left())

    def _partial_cardinality(self) -> Dict[PropertySet, int]:
        raise AssertionError('Union._partial_cardinality was called')

    def _order_collapse_unions(self, visited_sets: List[Set]=[]) -> Set:
        return self._order_nodes(Union,
                                 visited_sets,
                                 lambda x: _find_total_union(x))


class Difference(Set):
    """
    The difference of two Set expressions, by the identity A \ B = A ^ !B
    """
    def __new__(cls, *args, **kwargs) -> Set:
        assert len(args) == 2
        return Intersection(args[0], Complement(args[1]))


def flat_list_to_nested_expression(xs: List[Set], set_class) -> Set:
    """
    Reduces a list of Sets to a nested expression. E.g. if set_class is Intersection:
        [a, b, c] -> Intersection(Intersection(a, b), c)
    """
    if set_class == Intersection:
        unit = PropertySet()

    elif set_class == Union:
        unit = EmptySet()

    else:
        raise AssertionError('list_to_nested_expression can only be called with Union or Intersection')

    if len(xs) == 0:
        return unit

    elif len(xs) == 1:
        return xs[0]

    else:
        return reduce(set_class, xs[1:], xs[0])


def property_difference(x: PropertySet, y: PropertySet) -> PropertySet:
    assert y.is_superset_of(x)
    return PropertySet(*[p for p in x.properties if p not in y.properties])


def _add_dicts(x: Dict[PropertySet, int], y: Dict[PropertySet, int]) -> Dict[PropertySet, int]:
    """
    Adds the values of two dictionaries by key. Used to add terms in the calculation of cardinalities.
    """
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
        '_expand_complements',
        '_move_unions_to_left',
        '_move_intersections_to_left',
        '_order_expand_intersections'
    ]

    for n in range(len(sets)):
        left_combis = list(combinations(sets, n + 1))

        for left_combi in left_combis:
            rest = [x for x in sets if x not in left_combi]

            for m in range(len(rest)):
                right_combis = list(combinations(rest, m + 1))

                for right_combi in right_combis:
                    lhs = _apply_until_stable(
                        Complement(flat_list_to_nested_expression(left_combi, set_operation)),
                        expansion_methods
                    )

                    rhs = _apply_until_stable(
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
    sets = _unique(sets)
    if sets == [PropertySet()]:
        return PropertySet()

    sets = [x for x in sets if x != PropertySet()]

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
    if PropertySet() in sets:
        return PropertySet()

    # Filter out the Empty set and duplicates
    sets = _unique(sets)
    if sets == [EmptySet()]:
        return EmptySet()

    sets = [x for x in sets if x != EmptySet()]

    if _has_mutual_complements(sets, Union):
        return PropertySet()

    return flat_list_to_nested_expression(sorted(sets), Union)


def _unique(xs: Iterable) -> List:
    """
    For an iterable, returns the unique items in the same order in which they appear.
    """
    unique = []
    [unique.append(x) for x in xs if x not in unique]

    return unique


def _apply_until_stable(expr: Set, method):
    """
    Applies a method or a list of methods to a Set expression until the expression does not change anymore.
    """
    if isinstance(method, str):
        max_simplifications = 100
        i = 0

        previous_simplification = expr
        current_simplification = expr
        simplification_done = False

        while not simplification_done:
            i += 1
            if i > max_simplifications:
                raise RecursionError('_apply_until_stable does not converge with method ' + method + '\n' +
                                     ' current_expression: ' + str(current_simplification) + '\nprevious_expression: ' +
                                     str(previous_simplification))

            if not hasattr(current_simplification, method):
                raise AssertionError('_apply_until_stable ' + str(current_simplification) +
                                     ' missing method ' + str(method))

            current_simplification = getattr(current_simplification, method)()

            if current_simplification == previous_simplification:
                simplification_done = True
            else:
                previous_simplification = current_simplification

        return current_simplification
    elif isinstance(method, list):
        previous_simplification = expr
        current_simplification = expr
        simplification_done = False

        while not simplification_done:
            for m in method:
                current_simplification = _apply_until_stable(current_simplification, m)

            if current_simplification == previous_simplification:
                simplification_done = True
            else:
                previous_simplification = current_simplification

        return current_simplification
    else:
        raise AssertionError('_apply_until_stable called with expression ' + str(expr) + ' and method ' + str(method))
