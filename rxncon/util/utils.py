import functools
from typing import List
from enum import Enum
from rxncon.venntastic.sets import ValueSet, EmptySet, Union, Intersection, Complement

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
    if isinstance(set_expression, ValueSet):
        return ValueSet(leaf_transformer(set_expression.value))
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

def members(obj) -> List[str]:
    return [x for x in dir(obj) if not x.startswith('__')]

def elems_eq(first_list, second_list):
    if all(isinstance(x, list) for x in first_list) and all(isinstance(x, list) for x in second_list):
        uniq_first = [set(x) for x in first_list]
        uniq_second = [set(x) for x in second_list]

        return all(x in uniq_second for x in uniq_first) and all(x in uniq_first for x in uniq_second)

    else:
        return set(first_list) == set(second_list)

def all_type_eq(elems):
    assert elems
    the_type = type(elems[0])

    for elem in elems:
        if not isinstance(elem, the_type):
            return False

    return True

def all_eq(elems):
    assert elems
    the_val = elems[0]

    for elem in elems:
        if not elem == the_val:
            return False

    return True