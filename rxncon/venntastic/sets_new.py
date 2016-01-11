from typing import Optional, List
import itertools as itt


class Set:
    @property
    def canonical_form(self):
        result = self._complements_expanded()._unions_moved_to_left()._to_unions_of_intersections()
        result = _simplify_unions_of_intersections(result)

        return result

    @property
    def cardinality(self):
        return

    def is_equivalent_to(self, other: 'Set'):
        return self.canonical_form == other.canonical_form

    def is_superset_of(self, other: 'Set'):
        assert isinstance(other, Set)
        return False

    def is_subset_of(self, other: 'Set'):
        assert isinstance(other, Set)
        return False

    def _complements_expanded(self) -> 'Set':
        raise AssertionError

    def _unions_moved_to_left(self) -> 'Set':
        raise AssertionError

    def _to_unions_of_intersections(self) -> 'Set':
        raise AssertionError

    def _intersections_flattened(self):
        raise AssertionError


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

    def _complements_expanded(self):
        return self

    def _unions_moved_to_left(self):
        return self

    def _to_unions_of_intersections(self):
        return [self]

    def _intersections_flattened(self):
        return [self]


class EmptySet(Set):
    def __init__(self):
        pass

    def __eq__(self, other: Set) -> bool:
        assert isinstance(other, Set)

        if isinstance(other, EmptySet):
            return True

        else:
            return False

    def __hash__(self):
        return hash('*empty-set*')

    def __lt__(self, other):
        if isinstance(other, PropertySet) or isinstance(other, Complement):
            return True

        return False

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'EmptySet'

    def is_superset_of(self, other: Set):
        return self == other

    def is_subset_of(self, other: Set):
        return self == other

    def _complements_expanded(self):
        return self

    def _unions_moved_to_left(self):
        return self

    def _to_unions_of_intersections(self):
        return [self]

    def _intersections_flattened(self):
        return [self]


class Complement(Set):
    def __init__(self, expr: Set):
        self.expr = expr

    def __eq__(self, other: Set) -> bool:
        assert isinstance(other, Set)

        if isinstance(other, Complement):
            return self.expr == other.expr

        else:
            return False

    def __hash__(self):
        return hash('*complement-{}*'.format(hash(self.expr)))

    def __lt__(self, other):
        if isinstance(other, Complement):
            return self.expr < other.expr

        elif isinstance(other, PropertySet) or isinstance(other, EmptySet):
            return False

        else:
            raise AssertionError

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'Complement({})'.format(self.expr)

    def is_superset_of(self, other: Set):
        if self.is_equivalent_to(other):
            return True

        elif isinstance(other, Complement):
            return self.expr.is_subset_of(other.expr)

        else:
            return None

    def is_subset_of(self, other: Set):
        if self.is_equivalent_to(other):
            return True

        elif self.is_superset_of(other) is None:
            return None

        else:
            return not self.is_superset_of(other)

    def _complements_expanded(self):
        if isinstance(self.expr, Complement):
            return self.expr.expr._complements_expanded()

        elif isinstance(self.expr, PropertySet) and self.expr == PropertySet(None):
            return EmptySet()

        elif isinstance(self.expr, PropertySet) and self.expr != PropertySet(None):
            return self

        elif isinstance(self.expr, EmptySet):
            return PropertySet(None)

        elif isinstance(self.expr, Union):
            return Intersection(Complement(self.expr.left_expr), Complement(self.expr.right_expr))._complements_expanded()

        elif isinstance(self.expr, Intersection):
            return Union(Complement(self.expr.left_expr), Complement(self.expr.right_expr))._complements_expanded()

        else:
            raise AssertionError

    def _unions_moved_to_left(self):
        return Complement(self.expr._unions_moved_to_left())

    def _to_unions_of_intersections(self):
        return [self]

    def _intersections_flattened(self):
        return [self]

class BinarySet(Set):
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

    def _complements_expanded(self):
        return type(self)(self.left_expr._complements_expanded(), self.right_expr._complements_expanded())


class Intersection(BinarySet):
    def __hash__(self):
        return hash('*intersection-{0}{1}*'.format(hash(self.left_expr), hash(self.right_expr)))

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'Intersection({0}, {1})'.format(self.left_expr, self.right_expr)

    def is_superset_of(self, other: Set):
        if self.is_equivalent_to(other):
            return True

        elif self.is_subset_of(other) is None:
            return None

        else:
            return not self.is_subset_of(other)

    def is_subset_of(self, other: Set) -> Optional[bool]:
        if self.is_equivalent_to(other):
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

    def _to_unions_of_intersections(self):
        return [self]

    def _intersections_flattened(self):
        return self.left_expr._intersections_flattened() + self.right_expr._intersections_flattened()


class Union(BinarySet):
    def __hash__(self):
        return hash('*union-{0}{1}*'.format(hash(self.left_expr), hash(self.right_expr)))

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'Union({0}, {1})'.format(self.left_expr, self.right_expr)

    def is_superset_of(self, other: Set):
        if self.is_equivalent_to(other):
            return True

        elif self.left_expr.is_superset_of(other) or self.right_expr.is_superset_of(other):
            return True

        else:
            return None

    def is_subset_of(self, other: Set):
        if self.is_equivalent_to(other):
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
            return Union(Union(Union(self.left_expr.left_expr, self.left_expr.right_expr), self.right_expr.left_expr), self.right_expr.right_expr)\
                ._unions_moved_to_left()

        else:
            return Union(self.left_expr._unions_moved_to_left(), self.right_expr._unions_moved_to_left())

    def _to_unions_of_intersections(self) -> List[Set]:
        if isinstance(self.right_expr, Union):
            raise AssertionError

        if isinstance(self.left_expr, Union):
            left_contribution = self.left_expr._to_unions_of_intersections()

        else:
            left_contribution = [self.left_expr._intersections_flattened()]

        return left_contribution + [self.right_expr._intersections_flattened()]

    def _intersections_flattened(self):
        raise AssertionError


class Difference(Set):
    def __new__(cls, *args, **kwargs) -> Set:
        assert len(args) == 2
        return Intersection(args[0], Complement(args[1]))


def _simplify_unions_of_intersections(union_factors):
    # self is a list of lists that describes a Set expression of the form
    # Union(Union(A, B), C), where the inner lists are Intersections of
    # the Sets in the list.
    result = []

    for union_factor in union_factors:
        # union_factor is now a list of PropertySet / Complement(PropertySet), which we simplify by:
        # (1) if this contains the empty set, the union_factor is itself the empty set, and can therefore be discarded.
        if EmptySet() in union_factor:
            continue

        # (2) check whether there are a PropertySet and its complement in the list,
        if any([x == Complement(y) or Complement(x) == y for x in union_factor for y in union_factor]):
            continue

        # (3) removing duplicates from the list,
        cleaned_factor = list(set(union_factor))

        # (4) removing items from the list that are strict supersets of others.
        #cleaned_factor = [x for x in cleaned_factor for y in cleaned_factor
        #                  if x is not y and x.is_superset_of(y)]


        result.append(cleaned_factor)

    return result




