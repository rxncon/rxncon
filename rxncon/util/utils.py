import functools
from enum import Enum
from rxncon.venntastic.sets import PropertySet, EmptySet, Union, Intersection, Complement

def compose(*functions):
    return functools.reduce(lambda f, g: lambda x: f(g(x)), functions, lambda x: x)

class OrderedEnum(Enum):

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "{0}: {1}".format(self.name, self.value)

    def __lt__(self, other):
        if self.value is not None and other.value is not None:
            return self.value < other.value
        elif other.value is None:
            return False
        elif self.value is None:
            return True
        else:
            raise NotImplementedError


def transform_set_expression(set_expression, leaf_transformer):
    if isinstance(set_expression, PropertySet):
        return PropertySet(leaf_transformer(set_expression.value))
    elif isinstance(set_expression, Complement):
        return Complement(transform_set_expression(set_expression.expr, leaf_transformer))
    elif isinstance(set_expression, EmptySet):
        return EmptySet()
    elif isinstance(set_expression, Union):
        return Union(transform_set_expression(set_expression.left_expr, leaf_transformer),
                     transform_set_expression(set_expression.right_expr, leaf_transformer))
    elif isinstance(set_expression, Intersection):
        return Intersection(transform_set_expression(set_expression.left_expr, leaf_transformer),
                            transform_set_expression(set_expression.right_expr, leaf_transformer))
    else:
        raise NotImplementedError

