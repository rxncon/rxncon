import pytest
import itertools as itt

from rxncon.venntastic.sets import *


def test_property_set_construction():
    assert PropertySet('a')
    assert PropertySet(1)
    assert PropertySet(1.45)

    with pytest.raises(TypeError):
        PropertySet([1, 2, 3])


def test_property_set_comparison():
    assert PropertySet(1) == PropertySet(1)
    assert not PropertySet('a') == PropertySet('b')


def test_property_set_dictionary_keys():
    x = PropertySet(1)
    y = PropertySet(2)

    dictionary = {x: 'bla', y: 'diebla'}

    assert dictionary[PropertySet(1)] == 'bla'
    assert dictionary[x] == 'bla'
    assert dictionary[PropertySet(2)] == 'diebla'
    assert dictionary[y] == 'diebla'


def test_superset_subset_for_unary_sets():
    assert UniversalSet() == PropertySet(None)

    assert UniversalSet().is_superset_of(UniversalSet())
    assert UniversalSet().is_subset_of(UniversalSet())

    assert UniversalSet().is_superset_of(EmptySet())
    assert UniversalSet().is_superset_of(PropertySet(1))

    assert PropertySet(2).is_superset_of(PropertySet(2))
    assert PropertySet(2).is_subset_of(PropertySet(2))

    assert PropertySet(2).is_subset_of(UniversalSet())
    assert PropertySet(2).is_superset_of(EmptySet())
    assert EmptySet().is_subset_of(PropertySet(2))

    assert EmptySet().is_subset_of(EmptySet())
    assert EmptySet().is_superset_of(EmptySet())

    assert not PropertySet(2).is_subset_of(PropertySet(3))
    assert not PropertySet(2).is_superset_of(PropertySet(3))


def test_equivalence_complement_empty_and_universal_set():
    assert Complement(EmptySet()).is_equivalent_to(UniversalSet())
    assert Complement(UniversalSet()).is_equivalent_to(EmptySet())


def test_equivalence_respects_de_morgan_s_identities():
    x = PropertySet(1)
    y = PropertySet(2)

    assert Intersection(Complement(x), Complement(y)).is_equivalent_to(Complement(Union(x, y)))
    assert Union(Complement(x), Complement(y)).is_equivalent_to(Complement(Intersection(x, y)))


def test_equivalence_complement_squares_to_no_op():
    assert PropertySet(1).is_equivalent_to(Complement(Complement(PropertySet(1))))
    assert EmptySet().is_equivalent_to(Complement(Complement(EmptySet())))
    assert UniversalSet().is_equivalent_to(Complement(Complement(UniversalSet())))


def test_equivalence_intersection_with_complement_yields_empty_set():
    a = PropertySet(1)
    b = Complement(a)

    assert EmptySet().is_equivalent_to(Intersection(a, b))
    assert EmptySet().is_equivalent_to(Intersection(b, a))


def test_equivalence_intersection_of_property_sets():
    assert EmptySet().is_equivalent_to(Intersection(EmptySet(), PropertySet(1)))
    assert EmptySet().is_equivalent_to(Intersection(PropertySet(1), EmptySet()))

    # Intersection of a set with the universal set yields the set itself
    assert PropertySet(1).is_equivalent_to(Intersection(UniversalSet(), PropertySet(1)))
    assert PropertySet(1).is_equivalent_to(Intersection(PropertySet(1), UniversalSet()))

    # Intersection of set with itself yields the set
    assert EmptySet().is_equivalent_to(Intersection(EmptySet(), EmptySet()))
    assert UniversalSet().is_equivalent_to(Intersection(UniversalSet(), UniversalSet()))
    assert PropertySet(1).is_equivalent_to(Intersection(PropertySet(1), PropertySet(1)))


def test_canonical_form_intersection_of_unions():
    # Simplify should 'pull out' all unions
    assert Union(Intersection(PropertySet(1), PropertySet(3)), Intersection(PropertySet(2), PropertySet(3)))\
           .is_equivalent_to(Intersection(Union(PropertySet(1), PropertySet(2)), PropertySet(3)))

    assert Union(Intersection(PropertySet(1), PropertySet(3)), Intersection(PropertySet(2), PropertySet(3)))\
           .is_equivalent_to(Intersection(PropertySet(3), Union(PropertySet(1), PropertySet(2))))

    # assert Union(Union(Union(Intersection(PropertySet(1), PropertySet(3)), Intersection(PropertySet(1), PropertySet(4))),
    #                    Intersection(PropertySet(2), PropertySet(3))), Intersection(PropertySet(2), PropertySet(4)))\
    #        .is_equivalent_to(Intersection(Union(PropertySet(1), PropertySet(2)), Union(PropertySet(3), PropertySet(4))))


# def test_canonical_form_intersection_of_unions_with_complements():
    # # Test
    # # (!(1 ^ 2) u 3) ^ !(4 u 5) = (!1 ^ !4 ^ !5) u (!2 ^ !4 ^ !5) u (3 ^ !4 ^ !5)
    # a = Complement(Intersection(PropertySet(1), PropertySet(2)))
    # b = PropertySet(3)
    # c = Complement(Union(PropertySet(4), PropertySet(5)))
    # #
    # abc = Intersection(Union(a, b), c)
    # abc = Union(abc, Complement(abc))
    #
    # print(abc.equivalent_forms())
    #
    # print(Union(c, Complement(c)).equivalent_forms())
    # print(len(Union(c, Complement(c)).equivalent_forms()))

    # Test whether simplify can still guarantee basic properties of union and intersection
    # assert Union(abc, Complement(abc)).canonical_form() == UniversalSet()
    # assert Intersection(abc, Complement(abc)).canonical_form() == EmptySet()
    #
    # # Test whether simplify is idempotent
    # assert abc.canonical_form() == abc.canonical_form().canonical_form()


def test_simplify_difference_respects_absolute_relative_complement_identities():
    xs = [
        PropertySet(1),
        Complement(Union(PropertySet(2), PropertySet(3)))
    ]

    ys = [
        PropertySet(2), Union(PropertySet(1), Union(PropertySet(2), PropertySet(3))),
        Union(Complement(PropertySet(1)), Union(PropertySet(2), Complement(PropertySet(3))))
    ]

    for x, y in itt.product(xs, ys):
        assert Intersection(x, Complement(y)).is_equivalent_to(Difference(x, y))
        assert Union(Complement(x), y).is_equivalent_to(Complement(Difference(x, y)))


def test_simplify_union_with_complement_yields_universal_set():
    a = PropertySet(1)
    b = Complement(a)

    assert UniversalSet().is_equivalent_to(Union(a, b))
    assert Union(b, a).is_equivalent_to(UniversalSet())



def test_simplify_union_of_unions():
    a = Union(PropertySet(1), Union(PropertySet(2), Union(PropertySet(3), PropertySet(4))))
    b = Union(Union(Union(PropertySet(1), PropertySet(2)), PropertySet(3)), PropertySet(4))

    print(len(a.equivalent_forms()))

    for x in a.equivalent_forms():
        if str(x).startswith('Union(Union(Union('):
            print(x)

#
# def test_cardinality_empty_and_universal_set():
#     assert EmptySet().cardinality() == {}
#
#     assert PropertySet().cardinality() == {PropertySet(): 1}
#
#
# def test_cardinality_respects_complement_properties():
#     assert Intersection(PropertySet(1, 2, 3, 4), Complement(PropertySet(1, 2, 3, 4))).cardinality() == {}
#
#     assert Union(PropertySet(1, 2, 3, 4), Complement(PropertySet(1, 2, 3, 4))).cardinality() == {PropertySet(): 1}
#
#
# def test_cardinality_inclusion_exclusion_principle():
#     # Test the standard examples of the inclusion-exclusion principle:
#     # |A u (B u C)| = |A| + |B| + |C| - |A^B| - |A^C| - |B^C| + |A^B^C|
#     assert Union(PropertySet(1), Union(PropertySet(2), PropertySet(3))).cardinality() == \
#         {PropertySet(1): 1, PropertySet(2): 1, PropertySet(3): 1, PropertySet(1, 2): -1, PropertySet(1, 3): -1,
#          PropertySet(2, 3): -1, PropertySet(1, 2, 3): 1}
#
#     # |(A u B) u C| = same
#     assert Union(Union(PropertySet(1), PropertySet(2)), PropertySet(3)).cardinality() == \
#         {PropertySet(1): 1, PropertySet(2): 1, PropertySet(3): 1, PropertySet(1, 2): -1, PropertySet(1, 3): -1,
#          PropertySet(2, 3): -1, PropertySet(1, 2, 3): 1}
#
#     # |A u (B u (C u D))| = |A| + |B| + |C| + |D| - |A^B| - |A^C| - |A^D| - |B^C| - |B^D| -|C^D| + |A^B^C| +
#     #                       |B^C^D| + |A^C^D| + |A^B^D| - |A^B^C^D|
#     assert Union(PropertySet(1), Union(PropertySet(2), Union(PropertySet(3), PropertySet(4)))).cardinality() == \
#         {PropertySet(1): 1, PropertySet(2): 1, PropertySet(3): 1, PropertySet(4): 1,
#          PropertySet(1, 2): -1, PropertySet(1, 3): -1, PropertySet(1, 4): -1, PropertySet(2, 3): -1,
#          PropertySet(2, 4): -1, PropertySet(3, 4): -1,
#          PropertySet(1, 2, 3): 1, PropertySet(2, 3, 4): 1, PropertySet(1, 3, 4): 1, PropertySet(1, 2, 4): 1,
#          PropertySet(1, 2, 3, 4): -1}
#
#     # |(A u B) u (C u D)| = same
#     assert Union(Union(PropertySet(1), PropertySet(2)), Union(PropertySet(3), PropertySet(4))).cardinality() == \
#         {PropertySet(1): 1, PropertySet(2): 1, PropertySet(3): 1, PropertySet(4): 1,
#          PropertySet(1, 2): -1, PropertySet(1, 3): -1, PropertySet(1, 4): -1, PropertySet(2, 3): -1,
#          PropertySet(2, 4): -1, PropertySet(3, 4): -1,
#          PropertySet(1, 2, 3): 1, PropertySet(2, 3, 4): 1, PropertySet(1, 3, 4): 1, PropertySet(1, 2, 4): 1,
#          PropertySet(1, 2, 3, 4): -1}
#
#
# def test_cardinality_of_expressions_containing_complements():
#     assert Union(Complement(PropertySet(1)), PropertySet(2)).cardinality() == \
#         {PropertySet(): 1, PropertySet(1): -1, PropertySet(1, 2): 1}
#
#     assert Union(Complement(PropertySet(1)), Union(PropertySet(2), Complement(PropertySet(3)))).cardinality() == \
#         {PropertySet(): 1, PropertySet(1, 3): -1, PropertySet(1, 2, 3): 1}
#
#     assert Intersection(Complement(PropertySet(1)), Intersection(PropertySet(2), Complement(PropertySet(3)))).cardinality() == \
#         {PropertySet(2): 1, PropertySet(1, 2): -1, PropertySet(2, 3): -1, PropertySet(1, 2, 3): 1}
#
