import pytest
import itertools as itt
from typing import List

from rxncon.venntastic.sets_new import *


@pytest.fixture
def property_sets():
    properties = ['red', 1, 1.23, ('a', 3)]

    return [PropertySet(x) for x in properties]


@pytest.fixture
def set_expressions(property_sets):
    unary_expressions = property_sets + [Complement(x) for x in property_sets]

    binary_expressions = []

    for x, y in itt.product(unary_expressions, unary_expressions):
        binary_expressions.append(Intersection(x, y))
        binary_expressions.append(Union(x, y))

    ternary_expressions = []

    for x, y, z in itt.product(unary_expressions, unary_expressions, unary_expressions):
        ternary_expressions.append(Intersection(x, Intersection(y, z)))
        ternary_expressions.append(Intersection(x, Union(y, z)))
        ternary_expressions.append(Union(x, Intersection(y, z)))
        ternary_expressions.append(Union(x, Union(y, z)))

    ternary_expressions = ternary_expressions + [Complement(x) for x in ternary_expressions]

    return unary_expressions + binary_expressions + ternary_expressions


def test_property_set_comparison(set_expressions):
    for x in set_expressions:
        assert x == x


def test_canonical_form():
    p1, p2, p3, p4 = [PropertySet(x) for x in [1, 2, 3, 4]]

    x = Complement(Intersection(p4, Union(p3, Intersection(p2, p1))))

    print()
    print(x._complements_expanded())
    print(x._complements_expanded()._unions_moved_to_left())
    print(x._complements_expanded()._unions_moved_to_left()._unions_moved_to_left())
