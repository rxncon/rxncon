import pytest
import itertools as itt
from typing import List

from rxncon.venntastic.sets import *


# Test the basic properties of the Set data structure
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


def test_nested_list_simplifies_mutual_complements():
    x1 = Intersection(PropertySet(1), PropertySet(2))
    x2 = Intersection(PropertySet(1), Complement(PropertySet(2)))

    z = Union(x1, x2)

    assert z.to_nested_list_form() == [[PropertySet(1)]]


def test_nested_list_universal_empty():
    assert EmptySet().to_nested_list_form() == [[EmptySet()]]
    assert UniversalSet().to_nested_list_form() == [[UniversalSet()]]


def test_union_list_universal_empty():
    assert EmptySet().to_union_list_form() == [EmptySet()]
    assert UniversalSet().to_union_list_form() == [UniversalSet()]


def test_complementary_expansion():
    class Z3Integer:
        def __init__(self, value):
            assert isinstance(value, int)
            self.value = value % 3

        def __eq__(self, other):
            return self.value == other.value

        def __hash__(self):
            return hash(self.value)

        def __str__(self):
            return 'Z3Int({0})'.format(self.value)

        def __repr__(self):
            return str(self)

        def complements(self):
            return [Z3Integer(x) for x in range(3) if x != self.value]

    p = Complement(PropertySet(Z3Integer(1)))

    assert p.simplified_form().is_equivalent_to(Union(PropertySet(Z3Integer(2)), PropertySet(Z3Integer(0))))


# Test the superset / subset relationships
def test_superset_subset_for_unary_sets():
    assert UniversalSet() == PropertySet(None)
    assert PropertySet(None).is_superset_of(PropertySet(None))

    # UniversalSet == UniversalSet
    assert UniversalSet().is_superset_of(UniversalSet())
    assert UniversalSet().is_subset_of(UniversalSet())

    # UniversalSet contains all other sets
    assert UniversalSet().is_superset_of(EmptySet())
    assert UniversalSet().is_superset_of(PropertySet(1))

    # EmptySet is contained in all sets
    assert EmptySet().is_subset_of(EmptySet())
    assert EmptySet().is_superset_of(EmptySet())
    assert EmptySet().is_subset_of(PropertySet(2))

    # PropertySets <-> PropertySets
    assert PropertySet(2).is_superset_of(PropertySet(2))
    assert PropertySet(2).is_subset_of(PropertySet(2))
    assert not PropertySet(2).is_subset_of(PropertySet(3))
    assert not PropertySet(2).is_superset_of(PropertySet(3))

    # PropertySet <-> UniversalSet
    assert PropertySet(1).is_subset_of(UniversalSet())
    assert not UniversalSet().is_subset_of(PropertySet(1))
    assert UniversalSet().is_superset_of(PropertySet(1))
    assert not PropertySet(1).is_superset_of(UniversalSet())

    # PropertySet <-> EmptySet
    assert PropertySet(2).is_superset_of(EmptySet())
    assert not EmptySet().is_superset_of(PropertySet(2))
    assert EmptySet().is_subset_of(PropertySet(1))
    assert not PropertySet(1).is_subset_of(EmptySet())


def test_superset_subset_for_flat_intersections():
    assert Intersection(PropertySet(1), PropertySet(2)).is_subset_of(PropertySet(1))
    assert Intersection(PropertySet(1), PropertySet(2)).is_subset_of(PropertySet(2))

    assert not Intersection(PropertySet(1), PropertySet(2)).is_equivalent_to(PropertySet(1))
    assert not Intersection(PropertySet(1), PropertySet(2)).is_equivalent_to(PropertySet(2))

    assert PropertySet(1).is_superset_of(Intersection(PropertySet(1), PropertySet(2)))
    assert PropertySet(2).is_superset_of(Intersection(PropertySet(1), PropertySet(2)))


def test_superset_subset_for_nested_intersections():
    x = PropertySet(1)
    y = PropertySet(2)
    z = PropertySet(3)
    xy = Intersection(x, y)
    yz = Intersection(y, z)
    xz = Intersection(x, z)
    xyz = Intersection(x, Intersection(y, z))

    assert xyz.is_subset_of(x)
    assert xyz.is_subset_of(y)
    assert xyz.is_subset_of(z)
    assert xyz.is_subset_of(xy)
    assert xyz.is_subset_of(yz)
    assert xyz.is_subset_of(xz)

    assert x.is_superset_of(xyz)
    assert y.is_superset_of(xyz)
    assert z.is_superset_of(xyz)
    assert xy.is_superset_of(xyz)
    assert yz.is_superset_of(xyz)
    assert xz.is_superset_of(xyz)


def test_superset_subset_for_flat_unions():
    assert Union(PropertySet(1), PropertySet(2)).is_superset_of(PropertySet(1))
    assert Union(PropertySet(1), PropertySet(2)).is_superset_of(PropertySet(2))

    assert not Union(PropertySet(1), PropertySet(2)).is_equivalent_to(PropertySet(1))
    assert not Union(PropertySet(1), PropertySet(2)).is_equivalent_to(PropertySet(2))

    assert PropertySet(1).is_subset_of(Union(PropertySet(1), PropertySet(2)))
    assert PropertySet(2).is_subset_of(Union(PropertySet(1), PropertySet(2)))


def test_superset_subset_for_nested_unions():
    x = PropertySet(1)
    y = PropertySet(2)
    z = PropertySet(3)
    xy = Union(x, y)
    yz = Union(y, z)
    xz = Union(x, z)
    xyz = Union(x, Union(y, z))

    assert xyz.is_superset_of(x)
    assert xyz.is_superset_of(y)
    assert xyz.is_superset_of(z)
    assert xyz.is_superset_of(xy)
    assert xyz.is_superset_of(yz)
    assert xyz.is_superset_of(xz)

    assert x.is_subset_of(xyz)
    assert y.is_subset_of(xyz)
    assert z.is_subset_of(xyz)
    assert xy.is_subset_of(xyz)
    assert yz.is_subset_of(xyz)
    assert xz.is_subset_of(xyz)


# Test basic set-theoretic identities for a generated pool of set expressions
def test_de_morgan_s_identities(sets):
    for x, y in itt.product(sets, sets):
        assert Intersection(Complement(x), Complement(y)).is_equivalent_to(Complement(Union(x, y)))
        assert Union(Complement(x), Complement(y)).is_equivalent_to(Complement(Intersection(x, y)))


def test_complement_squares_to_no_op(sets):
    for x in sets:
        assert x.is_equivalent_to(Complement(Complement(x)))


def test_intersection_properties(sets):
    for x in sets:
        assert EmptySet().is_equivalent_to(Intersection(x, Complement(x)))
        assert EmptySet().is_equivalent_to(Intersection(Complement(x), x))

        assert EmptySet().is_equivalent_to(Intersection(EmptySet(), x))
        assert EmptySet().is_equivalent_to(Intersection(x, EmptySet()))

        assert x.is_equivalent_to(Intersection(UniversalSet(), x))
        assert x.is_equivalent_to(Intersection(x, UniversalSet()))

        assert x.is_equivalent_to(Intersection(x, x))


def test_union_properties(sets):
    for x in sets:
        assert UniversalSet().is_equivalent_to(Union(x, Complement(x)))
        assert UniversalSet().is_equivalent_to(Union(Complement(x), x))

        assert x.is_equivalent_to(Union(EmptySet(), x))
        assert x.is_equivalent_to(Union(x, EmptySet()))

        assert UniversalSet().is_equivalent_to(Union(UniversalSet(), x))
        assert UniversalSet().is_equivalent_to(Union(x, UniversalSet()))

        assert x.is_equivalent_to(Union(x, x))


def test_absolute_relative_complement_identities(sets):
    for x, y in itt.product(sets, sets):
        assert Intersection(x, Complement(y)).is_equivalent_to(Difference(x, y))
        assert Union(Complement(x), y).is_equivalent_to(Complement(Difference(x, y)))


def test_distributive_properties(sets):
    for x, y, z in itt.product(sets, sets, sets):
        assert Union(x, Intersection(y, z)).is_equivalent_to(Intersection(Union(x, y), Union(x, z)))
        assert Intersection(x, Union(y, z)).is_equivalent_to(Union(Intersection(x, y), Intersection(x, z)))


def test_absorption_properties(sets):
    for x, y in itt.product(sets, sets):
        assert x.is_equivalent_to(Union(x, Intersection(x, y)))
        assert x.is_equivalent_to(Intersection(x, Union(x, y)))


def test_is_equivalent_to():
    assert UniversalSet().is_equivalent_to(UniversalSet())
    assert EmptySet().is_equivalent_to(EmptySet())

    assert not UniversalSet().is_equivalent_to(PropertySet(1))
    assert not PropertySet(1).is_equivalent_to(UniversalSet())

    assert UniversalSet().is_equivalent_to(Union(PropertySet(1), Complement(PropertySet(1))))


# Test the cardinality calculator
def test_cardinality_empty_and_universal_set():
    assert EmptySet().cardinality == {}
    assert UniversalSet().cardinality == {(UniversalSet(),): 1}


def test_cardinality_respects_complement_properties(sets):
    for x in sets:
        assert Intersection(x, Complement(x)).cardinality == {}
        assert Union(x, Complement(x)).cardinality == {(UniversalSet(),): 1}


def test_cardinality_property_sets():
    # Test the standard examples of the inclusion-exclusion principle:
    # |A u (B u C)| = |A| + |B| + |C| - |A^B| - |A^C| - |B^C| + |A^B^C|
    assert Union(PropertySet(1), Union(PropertySet(2), PropertySet(3))).cardinality == \
        {(PropertySet(1),): 1, (PropertySet(2),): 1, (PropertySet(3),): 1, (PropertySet(1), PropertySet(2)): -1, (PropertySet(1), PropertySet(3)): -1,
         (PropertySet(2), PropertySet(3)): -1, (PropertySet(1), PropertySet(2), PropertySet(3)): 1}

    # |(A u B) u C| = same
    assert Union(Union(PropertySet(1), PropertySet(2)), PropertySet(3)).cardinality == \
        {(PropertySet(1),): 1, (PropertySet(2),): 1, (PropertySet(3),): 1, (PropertySet(1), PropertySet(2)): -1, (PropertySet(1), PropertySet(3)): -1,
         (PropertySet(2), PropertySet(3)): -1, (PropertySet(1), PropertySet(2), PropertySet(3)): 1}

    # |A u (B u (C u D))| = |A| + |B| + |C| + |D| - |A^B| - |A^C| - |A^D| - |B^C| - |B^D| -|C^D| + |A^B^C| +
    #                       |B^C^D| + |A^C^D| + |A^B^D| - |A^B^C^D|
    assert Union(PropertySet(1), Union(PropertySet(2), Union(PropertySet(3), PropertySet(4)))).cardinality == \
        {(PropertySet(1),): 1, (PropertySet(2),): 1, (PropertySet(3),): 1, (PropertySet(4),): 1,
         (PropertySet(1), PropertySet(2)): -1, (PropertySet(1), PropertySet(3)): -1, (PropertySet(1), PropertySet(4)): -1,
         (PropertySet(2), PropertySet(3)): -1, (PropertySet(2), PropertySet(4)): -1, (PropertySet(3), PropertySet(4)): -1,
         (PropertySet(1), PropertySet(2), PropertySet(3)): 1, (PropertySet(2), PropertySet(3), PropertySet(4)): 1,
         (PropertySet(1), PropertySet(3), PropertySet(4)): 1, (PropertySet(1), PropertySet(2), PropertySet(4)): 1,
         (PropertySet(1), PropertySet(2), PropertySet(3), PropertySet(4)): -1}

    # |(A u B) u (C u D)| = same
    assert Union(Union(PropertySet(1), PropertySet(2)), Union(PropertySet(3), PropertySet(4))).cardinality == \
        {(PropertySet(1),): 1, (PropertySet(2),): 1, (PropertySet(3),): 1, (PropertySet(4),): 1,
         (PropertySet(1), PropertySet(2)): -1, (PropertySet(1), PropertySet(3)): -1, (PropertySet(1), PropertySet(4)): -1,
         (PropertySet(2), PropertySet(3)): -1, (PropertySet(2), PropertySet(4)): -1, (PropertySet(3), PropertySet(4)): -1,
         (PropertySet(1), PropertySet(2), PropertySet(3)): 1, (PropertySet(2), PropertySet(3), PropertySet(4)): 1,
         (PropertySet(1), PropertySet(3), PropertySet(4)): 1, (PropertySet(1), PropertySet(2), PropertySet(4)): 1,
         (PropertySet(1), PropertySet(2), PropertySet(3), PropertySet(4)): -1}


def test_cardinality_property_sets_and_complements():
    assert Union(Complement(PropertySet(1)), PropertySet(2)).cardinality == \
        {(UniversalSet(),): 1, (PropertySet(1),): -1, (PropertySet(1), PropertySet(2)): 1}

    assert Union(Complement(PropertySet(1)), Union(PropertySet(2), Complement(PropertySet(3)))).cardinality == \
        {(UniversalSet(),): 1, (PropertySet(1), PropertySet(3)): -1, (PropertySet(1), PropertySet(2), PropertySet(3)): 1}

    assert Intersection(Complement(PropertySet(1)), Intersection(PropertySet(2), Complement(PropertySet(3)))).cardinality == \
        {(PropertySet(2),): 1, (PropertySet(1), PropertySet(2)): -1, (PropertySet(2), PropertySet(3)): -1,
         (PropertySet(1), PropertySet(2), PropertySet(3)): 1}


# Test the boolean functions that are used to determine super/subset relations
def test_boolean_function_single_and_clause():
    bool_func = BooleanFunction([BooleanAndClause(required_true=[0, 2], required_false=[1, 3])])

    assert bool_func(true=[0, 2], false=[1, 3])
    assert not bool_func(true=[0, 1, 2, 3])
    assert not bool_func(true=[1, 3], false=[0, 2])


def test_boolean_function_multiple_and_clauses():
    bool_func = BooleanFunction([BooleanAndClause(required_true=[0], required_false=[1]),
                                 BooleanAndClause(required_true=[1], required_false=[0])])

    assert not bool_func(true=[0, 1])
    assert bool_func(true=[0], false=[1])
    assert not bool_func(false=[0, 1])
    assert bool_func(true=[1], false=[0])

    bool_func = BooleanFunction([BooleanAndClause(required_true=['a']),
                                 BooleanAndClause(required_false=['a'])])

    assert bool_func(true=['a'])
    assert bool_func(false=['a'])


def test_boolean_function_always_true():
    # Non-trivial boolean function with one input.
    bool_func = BooleanFunction([BooleanAndClause(required_true=[0])])

    assert bool_func(true=[0])
    assert not bool_func(false=[0])
    assert not bool_func()

    assert not BooleanFunctionAlwaysTrue().implies(bool_func)
    assert bool_func.implies(BooleanFunctionAlwaysTrue())

    # Trivial 'always true' boolean function with one input.
    bool_func = BooleanFunction([BooleanAndClause(required_true=[0]),
                                 BooleanAndClause(required_false=[0])])

    assert bool_func.implies(BooleanFunctionAlwaysTrue())
    assert BooleanFunctionAlwaysTrue().implies(bool_func)


def test_boolean_function_always_false():
    # Non-trivial boolean function with one input.
    bool_func = BooleanFunction([BooleanAndClause(required_true=[0])])

    assert BooleanFunctionAlwaysFalse().implies(bool_func)
    assert not bool_func.implies(BooleanFunctionAlwaysFalse())


def test_boolean_function_called_with_wrong_arguments():
    bool_func = BooleanFunction([BooleanAndClause(required_true=['a']),
                                 BooleanAndClause(required_false=['a'])])

    with pytest.raises(NameError):
        bool_func(true=['b'])


def test_boolean_function_from_nested_list_form():
    bool_func = boolean_function_from_nested_list_form([[PropertySet(1), Complement(PropertySet(2))],
                                                        [PropertySet(2), Complement(PropertySet(1))]])

    assert bool_func(true=[PropertySet(1)], false=[PropertySet(2)])
    assert bool_func(true=[PropertySet(2)], false=[PropertySet(1)])

    assert not bool_func(true=[PropertySet(1), PropertySet(2)])
    assert not bool_func(false=[PropertySet(1), PropertySet(2)])

    assert not bool_func(true=[PropertySet(1)])
    assert not bool_func(true=[PropertySet(2)])
    assert not bool_func(false=[PropertySet(1)])
    assert not bool_func(false=[PropertySet(2)])


def test_boolean_function_from_nested_list_form_corner_cases():
    assert isinstance(boolean_function_from_nested_list_form([[PropertySet(1), Complement(PropertySet(1))]]),
                      BooleanFunctionAlwaysFalse)

    with pytest.raises(Exception):
        boolean_function_from_nested_list_form([])

    with pytest.raises(Exception):
        boolean_function_from_nested_list_form('a')


def test_generate_boolean_value_lists():
    value_generator = generate_boolean_value_lists([1, 2, 3])

    expected_boolean_lists = [
        ([1, 2, 3], []),
        ([1, 2], [3]),
        ([1, 3], [2]),
        ([1], [2, 3]),
        ([2, 3], [1]),
        ([2], [1, 3]),
        ([3], [1, 2]),
        ([], [1, 2, 3])
    ]

    actual_boolean_lists = list(value_generator)

    assert all(expected in actual_boolean_lists for expected in expected_boolean_lists)
    assert all(actual in expected_boolean_lists for actual in actual_boolean_lists)


def test_boolean_function_implies_equivalent_functions():
    first_bool_func = BooleanFunction([BooleanAndClause(required_true=[0], required_false=[1]),
                                       BooleanAndClause(required_true=[1], required_false=[0])])

    second_bool_func = BooleanFunction([BooleanAndClause(required_true=[1], required_false=[0]),
                                        BooleanAndClause(required_true=[0], required_false=[1])])

    assert first_bool_func.implies(second_bool_func)
    assert second_bool_func.implies(first_bool_func)


def test_boolean_function_implies_subset():
    and_func = BooleanFunction([BooleanAndClause(required_true=[0, 1])])
    zero_func = BooleanFunction([BooleanAndClause(required_true=[0])])
    one_func = BooleanFunction([BooleanAndClause(required_true=[1])])

    assert and_func.implies(zero_func)
    assert and_func.implies(one_func)
    assert not zero_func.implies(and_func)
    assert not one_func.implies(and_func)


# Test the "Gram-Schmidt" algorithm
def test_manual_gram_schmidt_overlaps_simplify_correctly():
    ab = PropertySet('ab')
    ac = PropertySet('ac')
    be = PropertySet('be')
    bf = PropertySet('bf')

    set1 = Intersection(ab, be)
    set2 = Intersection(Intersection(ab, bf), Complement(set1))
    set3 = Intersection(Intersection(Intersection(ac, be), Complement(set2)), Complement(set1))
    set4 = Intersection(Intersection(Intersection(Intersection(ac, bf), Complement(set3)), Complement(set2)), Complement(set1))

    assert set2.is_equivalent_to(Intersection(Intersection(ab, bf), Complement(be)))
    assert set3.is_equivalent_to(Intersection(Intersection(ac, be), Complement(ab)))
    assert set4.is_equivalent_to(Intersection(Intersection(Intersection(ac, bf), Complement(be)), Complement(ab)))

    assert len(set1.to_nested_list_form()) == 1
    assert len(set2.to_nested_list_form()) == 1
    assert len(set3.to_nested_list_form()) == 1
    assert len(set4.to_nested_list_form()) == 1


def test_gram_schmidt():
    ab = PropertySet('ab')
    ac = PropertySet('ac')
    be = PropertySet('be')
    bf = PropertySet('bf')

    ab_be = Intersection(ab, be)
    ab_bf = Intersection(ab, bf)
    ac_be = Intersection(ac, be)
    ac_bf = Intersection(ac, bf)

    set1, set2, set3, set4 = gram_schmidt_disjunctify([ab_be, ab_bf, ac_be, ac_bf])

    assert set1.is_equivalent_to(ab_be)
    assert set2.is_equivalent_to(Intersection(Intersection(ab, bf), Complement(be)))
    assert set3.is_equivalent_to(Intersection(Intersection(ac, be), Complement(ab)))
    assert set4.is_equivalent_to(Intersection(Intersection(Intersection(ac, bf), Complement(be)), Complement(ab)))


@pytest.fixture
def sets():
    return [
        EmptySet(),
        PropertySet(1),
        UniversalSet(),
        Union(PropertySet(1), PropertySet(2)),
        Intersection(PropertySet(1), PropertySet(2)),
        Intersection(PropertySet(1), Complement(PropertySet(2))),
        Union(Intersection(PropertySet(1), PropertySet(2)), PropertySet(3)),
        Union(Intersection(PropertySet(1), PropertySet(2)), Intersection(PropertySet(3), PropertySet(4))),
        Union(Complement(Union(PropertySet(1), Complement(PropertySet(2)))), Intersection(PropertySet(3), PropertySet(4)))
    ]
