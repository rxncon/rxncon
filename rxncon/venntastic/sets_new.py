

class Set:
    @property
    def canonical_form(self) -> 'Set':
        pass

    @property
    def cardinality(self):
        return

    def is_superset_of(self, other: 'Set'):
        assert isinstance(other, Set)
        return False

    def is_subset_of(self, other: 'Set'):
        assert isinstance(other, Set)
        return False


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

    def is_superset_of(self, other: Set):
        return self == other

    def is_subset_of(self, other: Set):
        return self == other


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

    def is_superset_of(self, other: Set):
        return self == other

    def is_subset_of(self, other: Set):
        return self == other


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

    def is_superset_of(self, other: Set):
        return False

    def is_subset_of(self, other: Set):
        return False


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


class Intersection(BinarySet):
    def __hash__(self):
        return hash('*intersection-{0}{1}*'.format(hash(self.left_expr), hash(self.right_expr)))


class Union(BinarySet):
    def __hash__(self):
        return hash('*union-{0}{1}*'.format(hash(self.left_expr), hash(self.right_expr)))


class Difference(Set):
    def __new__(cls, *args, **kwargs) -> Set:
        assert len(args) == 2
        return Intersection(args[0], Complement(args[1]))