from rxncon.venntastic.sets import *

import itertools as itt

def test_property_set_construction():
    x = PropertySet()
    assert x.properties == []

    x = PropertySet(1, 2, 3, 4, 5)
    assert x.properties == [1, 2, 3, 4, 5]

    x = PropertySet(1, 2, 2, 2, 3)
    assert x.properties == [1, 2, 3]


def test_property_set_comparison():
    x = PropertySet(1, 2, 3)
    y = PropertySet(1, 2, 3)
    assert x == y

    x = PropertySet(1, 2)
    y = PropertySet(2, 3)
    assert not x == y

    x = PropertySet(1)
    y = PropertySet(1, 2)
    assert not x == y


def test_property_set_hashable():
    x = PropertySet(1, 2, 3)
    y = PropertySet(1)

    dictionary = {x: 'bla', y: 'diebla'}

    assert dictionary[PropertySet(1, 2, 3)] == 'bla'
    assert dictionary[x] == 'bla'
    assert dictionary[PropertySet(1)] == 'diebla'
    assert dictionary[y] == 'diebla'


def test_property_set_contains_other_property_set():
    # The universal set contains itself
    assert PropertySet().contains(PropertySet())

    # The universal set contains the empty set
    assert PropertySet().contains(EmptySet())

    # The universal set contains any given set
    assert PropertySet().contains(PropertySet(1, 5, 6))

    # Set with fewer specified properties contains one with more specified properties
    assert PropertySet(1).contains(PropertySet(1, 2, 3))

    # Set with more specified properties does not contain one with fewer specified properties
    assert not PropertySet(2, 3).contains(PropertySet(2))

    # A 'larger' set with the 'wrong' properties does not contain a more specific set
    assert not PropertySet(1).contains(PropertySet(2, 3))

    # A set contains itself
    assert PropertySet(1).contains(PropertySet(1))

    # Empty set definitions
    assert PropertySet().contains(EmptySet())
    assert PropertySet(4).contains(EmptySet())
    assert EmptySet().contains(EmptySet())
    assert not EmptySet().contains(PropertySet(3))
    assert not EmptySet().contains(PropertySet())


def test_simplify_complement_empty_and_universal_set():
    assert Complement(EmptySet()).simplify() == PropertySet()
    assert Complement(PropertySet()).simplify() == EmptySet()


def test_simplify_respects_de_morgan_s_identities():
    x = PropertySet(1)
    y = PropertySet(4)

    assert Intersection(Complement(x), Complement(y)) == Complement(Union(x, y)).simplify()
    assert Union(Complement(x), Complement(y)) == Complement(Intersection(x, y)).simplify()


def test_simplify_complement_squares_to_no_op():
    assert PropertySet(1, 2, 3) == Complement(Complement(PropertySet(1, 2, 3))).simplify()

    assert EmptySet() == Complement(Complement(EmptySet())).simplify()

    assert PropertySet() == Complement(Complement(PropertySet())).simplify()


def test_simplify_intersection_with_complement_yields_empty_set():
    a = PropertySet(1, 2)
    b = Complement(a)

    assert EmptySet() == Intersection(a, b).simplify()
    assert EmptySet() == Intersection(b, a).simplify()


def test_simplify_intersection_of_property_sets():
    # Intersection with the empty set yields the empty set
    assert EmptySet() == Intersection(EmptySet(), PropertySet(1)).simplify()
    assert EmptySet() == Intersection(PropertySet(1), EmptySet()).simplify()

    # Intersection of a set with the universal set yields the set itself
    assert PropertySet(2, 3, 4) == Intersection(PropertySet(), PropertySet(2, 3, 4)).simplify()
    assert PropertySet(2, 3, 4) == Intersection(PropertySet(2, 3, 4), PropertySet()).simplify()

    # Intersection of set with itself yields the set
    assert EmptySet() == Intersection(EmptySet(), EmptySet()).simplify()
    assert PropertySet() == Intersection(PropertySet(), PropertySet()).simplify()
    assert PropertySet(2) == Intersection(PropertySet(2), PropertySet(2)).simplify()

    # Intersection of two property sets have combination of properties
    assert PropertySet(1, 2) == Intersection(PropertySet(1), PropertySet(2)).simplify()
    assert PropertySet(1, 2, 3, 4, 5) == Intersection(PropertySet(1, 2, 3), PropertySet(3, 4, 5)).simplify()


def test_simplify_intersection_of_unions():
    # Simplify should 'pull out' all unions
    assert Union(PropertySet(1, 3), PropertySet(2, 3)) == \
           Intersection(Union(PropertySet(1), PropertySet(2)), PropertySet(3)).simplify()

    assert Union(PropertySet(1, 3), PropertySet(2, 3)) == \
           Intersection(PropertySet(3), Union(PropertySet(1), PropertySet(2))).simplify()

    assert Union(Union(Union(PropertySet(1, 3), PropertySet(1, 4)), PropertySet(2, 3)), PropertySet(2, 4)) == \
           Intersection(Union(PropertySet(1), PropertySet(2)), Union(PropertySet(3), PropertySet(4))).simplify()


def test_simplify_intersection_of_unions_with_complements():
    # Test
    # (!(1 ^ 2) u 3) ^ !(4 u 5) = (!1 ^ !4 ^ !5) u (!2 ^ !4 ^ !5) u (3 ^ !4 ^ !5)
    a = Complement(Intersection(PropertySet(1), PropertySet(2)))
    b = PropertySet(3)
    c = Complement(Union(PropertySet(4), PropertySet(5)))

    abc = Intersection(Union(a, b), c)

    # Test whether simplify can still guarantee basic properties of union and intersection
    assert Union(abc, Complement(abc)).simplify() == PropertySet()
    assert Intersection(abc, Complement(abc)).simplify() == EmptySet()

    # Test whether simplify is idempotent
    assert abc.simplify() == abc.simplify().simplify()


def test_simplify_difference_respects_absolute_relative_complement_identities():
    xs = [
        PropertySet(1, 2, 3),
        Complement(Union(PropertySet(3, 4), PropertySet(4, 5)))
    ]

    ys = [
        PropertySet(2), Union(PropertySet(1), Union(PropertySet(2), PropertySet(3))),
        Union(Complement(PropertySet(1)), Union(PropertySet(2), Complement(PropertySet(3))))
    ]

    for x, y in itt.product(xs, ys):
        assert Intersection(x, Complement(y)).simplify() == Difference(x, y).simplify()
        assert Union(Complement(x), y).simplify() == Complement(Difference(x, y)).simplify()


def test_simplify_union_with_complement_yields_universal_set():
    a = PropertySet(1, 2)
    b = Complement(a)

    assert PropertySet() == Union(a, b).simplify()
    assert PropertySet() == Union(b, a).simplify()


def test_simplify_union_of_property_sets():
    assert PropertySet(1) == Union(EmptySet(), PropertySet(1)).simplify()

    assert PropertySet(1) == Union(PropertySet(1), EmptySet()).simplify()

    assert PropertySet() == Union(PropertySet(), PropertySet(1, 2)).simplify()


def test_simplify_union_of_unions():
    a = Union(PropertySet(1), Union(PropertySet(2), Union(PropertySet(3), PropertySet(4))))

    assert Union(Union(Union(PropertySet(1), PropertySet(2)), PropertySet(3)), PropertySet(4)) == a.simplify()


def test_cardinality_empty_and_universal_set():
    assert EmptySet().cardinality() == {}

    assert PropertySet().cardinality() == {PropertySet(): 1}


def test_cardinality_respects_complement_properties():
    assert Intersection(PropertySet(1, 2, 3, 4), Complement(PropertySet(1, 2, 3, 4))).cardinality() == {}

    assert Union(PropertySet(1, 2, 3, 4), Complement(PropertySet(1, 2, 3, 4))).cardinality() == {PropertySet(): 1}


def test_cardinality_inclusion_exclusion_principle():
    # Test the standard examples of the inclusion-exclusion principle:
    # |A u (B u C)| = |A| + |B| + |C| - |A^B| - |A^C| - |B^C| + |A^B^C|
    assert Union(PropertySet(1), Union(PropertySet(2), PropertySet(3))).cardinality() == \
        {PropertySet(1): 1, PropertySet(2): 1, PropertySet(3): 1, PropertySet(1, 2): -1, PropertySet(1, 3): -1,
         PropertySet(2, 3): -1, PropertySet(1, 2, 3): 1}

    # |(A u B) u C| = same
    assert Union(Union(PropertySet(1), PropertySet(2)), PropertySet(3)).cardinality() == \
        {PropertySet(1): 1, PropertySet(2): 1, PropertySet(3): 1, PropertySet(1, 2): -1, PropertySet(1, 3): -1,
         PropertySet(2, 3): -1, PropertySet(1, 2, 3): 1}

    # |A u (B u (C u D))| = |A| + |B| + |C| + |D| - |A^B| - |A^C| - |A^D| - |B^C| - |B^D| -|C^D| + |A^B^C| +
    #                       |B^C^D| + |A^C^D| + |A^B^D| - |A^B^C^D|
    assert Union(PropertySet(1), Union(PropertySet(2), Union(PropertySet(3), PropertySet(4)))).cardinality() == \
        {PropertySet(1): 1, PropertySet(2): 1, PropertySet(3): 1, PropertySet(4): 1,
         PropertySet(1, 2): -1, PropertySet(1, 3): -1, PropertySet(1, 4): -1, PropertySet(2, 3): -1,
         PropertySet(2, 4): -1, PropertySet(3, 4): -1,
         PropertySet(1, 2, 3): 1, PropertySet(2, 3, 4): 1, PropertySet(1, 3, 4): 1, PropertySet(1, 2, 4): 1,
         PropertySet(1, 2, 3, 4): -1}

    # |(A u B) u (C u D)| = same
    assert Union(Union(PropertySet(1), PropertySet(2)), Union(PropertySet(3), PropertySet(4))).cardinality() == \
        {PropertySet(1): 1, PropertySet(2): 1, PropertySet(3): 1, PropertySet(4): 1,
         PropertySet(1, 2): -1, PropertySet(1, 3): -1, PropertySet(1, 4): -1, PropertySet(2, 3): -1,
         PropertySet(2, 4): -1, PropertySet(3, 4): -1,
         PropertySet(1, 2, 3): 1, PropertySet(2, 3, 4): 1, PropertySet(1, 3, 4): 1, PropertySet(1, 2, 4): 1,
         PropertySet(1, 2, 3, 4): -1}


def test_cardinality_of_expressions_containing_complements():
    assert Union(Complement(PropertySet(1)), PropertySet(2)).cardinality() == \
        {PropertySet(): 1, PropertySet(1): -1, PropertySet(1, 2): 1}

    assert Union(Complement(PropertySet(1)), Union(PropertySet(2), Complement(PropertySet(3)))).cardinality() == \
        {PropertySet(): 1, PropertySet(1, 3): -1, PropertySet(1, 2, 3): 1}

    assert Intersection(Complement(PropertySet(1)), Intersection(PropertySet(2), Complement(PropertySet(3)))).cardinality() == \
        {PropertySet(2): 1, PropertySet(1, 2): -1, PropertySet(2, 3): -1, PropertySet(1, 2, 3): 1}

