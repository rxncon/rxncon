import pytest
import itertools as itt
from typing import List
from copy import deepcopy


from rxncon.venntastic.sets import ValueSet, Union, Intersection, Complement, EmptySet, UniversalSet, Difference, venn_from_str, Set, \
    DisjunctiveUnion


def test_property_set_construction() -> None:
    assert ValueSet('a')
    assert ValueSet(1)
    assert ValueSet(1.45)

    with pytest.raises(TypeError):
        ValueSet([1, 2, 3])


def test_property_set_comparison() -> None:
    assert ValueSet(1) == ValueSet(1)
    assert not ValueSet('a') == ValueSet('b')


def test_property_set_dictionary_keys() -> None:
    x = ValueSet(1)
    y = ValueSet(2)

    dictionary = {x: 'bla', y: 'diebla'}

    assert dictionary[ValueSet(1)] == 'bla'
    assert dictionary[x] == 'bla'
    assert dictionary[ValueSet(2)] == 'diebla'
    assert dictionary[y] == 'diebla'


def test_simplifies() -> None:
    x1 = Intersection(ValueSet(1), ValueSet(2))
    x2 = Intersection(ValueSet(1), Complement(ValueSet(2)))

    assert Union(x1, x2).is_equivalent_to(ValueSet(1))


def test_parser() -> None:
    assert venn_from_str('1 & 2', int).is_equivalent_to(Intersection(ValueSet(1), ValueSet(2)))
    assert venn_from_str('1 & 2', str).is_equivalent_to(Intersection(ValueSet('1'), ValueSet('2')))
    assert venn_from_str('( 1 | 2 ) & 3', int).is_equivalent_to(Intersection(ValueSet(3), Union(ValueSet(1), ValueSet(2))))
    assert venn_from_str('~ 1', int).is_equivalent_to(Complement(ValueSet(1)))
    assert venn_from_str('~( 1 | 2 )', int).is_equivalent_to(Intersection(Complement(ValueSet(1)), Complement(ValueSet(2))))


def test_values() -> None:
    assert set(venn_from_str('1 | 2 | 3 | 4 | 2', int).values) == {1, 2, 3, 4}


def test_dnf_form() -> None:
    assert venn_from_str('1 & 2', int).to_dnf_set().is_equivalent_to(venn_from_str('1 & 2', int))
    assert venn_from_str('1 | 2', int).to_dnf_set().is_equivalent_to(venn_from_str('1 | 2', int))
    assert UniversalSet().to_dnf_set().is_equivalent_to(UniversalSet())
    assert EmptySet().to_dnf_set().is_equivalent_to(EmptySet())


def test_nary_sets_constructor() -> None:
    assert Union() == EmptySet()
    assert Intersection() == UniversalSet()
    assert DisjunctiveUnion() == EmptySet()

    assert Union(venn_from_str('1', int)).is_equivalent_to(venn_from_str('1', int))
    assert Intersection(venn_from_str('1', int)).is_equivalent_to(venn_from_str('1', int))
    assert DisjunctiveUnion(venn_from_str('1', int)).is_equivalent_to(venn_from_str('1', int))


def test_deepcopy(sets: List[Set]) -> None:
    for x in sets:
        assert x.is_equivalent_to(deepcopy(x))


def test_eval_boolean_func_OR() -> None:
    f = venn_from_str('1 | 2', int)
    assert f.eval_boolean_func({1: True, 2: True}) is True
    assert f.eval_boolean_func({1: True, 2: False}) is True
    assert f.eval_boolean_func({1: False, 2: True}) is True
    assert f.eval_boolean_func({1: False, 2: False}) is False

    with pytest.raises(AssertionError):
        f.eval_boolean_func({1: True})

    with pytest.raises(AssertionError):
        f.eval_boolean_func({1: True, 3: False})


def test_eval_boolean_func_AND() -> None:
    f = venn_from_str('1 & 2', int)
    assert f.eval_boolean_func({1: True, 2: True}) is True
    assert f.eval_boolean_func({1: True, 2: False}) is False
    assert f.eval_boolean_func({1: False, 2: True}) is False
    assert f.eval_boolean_func({1: False, 2: False}) is False

    with pytest.raises(AssertionError):
        f.eval_boolean_func({1: True})

    with pytest.raises(AssertionError):
        f.eval_boolean_func({1: True, 3: False})


def test_list_form() -> None:
    assert venn_from_str('1', int).to_dnf_list() == [venn_from_str('1', int)]

    assert venn_from_str('1 & 2', int).to_dnf_list() == [venn_from_str('1 & 2', int)]
    assert set(venn_from_str('1 | 2', int).to_dnf_list()) == {ValueSet(1), ValueSet(2)}

    x = venn_from_str('1 & ( 2 | 3 )', int)
    assert any(elem.is_equivalent_to(venn_from_str('1 & 2', int)) for elem in x.to_dnf_list())
    assert any(elem.is_equivalent_to(venn_from_str('1 & 3', int)) for elem in x.to_dnf_list())

    assert UniversalSet().to_dnf_list() == [UniversalSet()]
    assert EmptySet().to_dnf_list() == [EmptySet()]


def test_nested_list_form() -> None:
    assert venn_from_str('1', int).to_dnf_nested_list() == [[venn_from_str('1', int)]]

    x = venn_from_str('1 & ( 2 | 3 )', int)
    assert [ValueSet(1), ValueSet(2)] in x.to_dnf_nested_list()
    assert [ValueSet(1), ValueSet(3)] in x.to_dnf_nested_list()

    assert venn_from_str('1 & 2', int).to_dnf_nested_list() == [[ValueSet(1), ValueSet(2)]]
    assert venn_from_str('1 | 2', int).to_dnf_nested_list() == [[ValueSet(1)], [ValueSet(2)]]

    assert UniversalSet().to_dnf_nested_list() == [[UniversalSet()]]
    assert EmptySet().to_dnf_nested_list() == [[EmptySet()]]


def test_calc_solutions() -> None:
    # Contradiction should give no solutions.
    assert venn_from_str('( a ) & ~( a )', str).calc_solutions() == []
    # Tautology should give empty dict as solution.
    assert venn_from_str('( a ) | ~( a )', str).calc_solutions() == [{}]

    assert {'a': True, 'b': False} in venn_from_str('a | ~( b )', str).calc_solutions()
    assert {'a': True, 'b': True} in venn_from_str('a | ~( b )', str).calc_solutions()
    assert {'a': False, 'b': False} in venn_from_str('a | ~( b )', str).calc_solutions()


# Test the superset / subset relationships
def test_superset_subset_for_unary_sets() -> None:
    # UniversalSet == UniversalSet
    assert UniversalSet().is_superset_of(UniversalSet())
    assert UniversalSet().is_subset_of(UniversalSet())

    # UniversalSet contains all other sets
    assert UniversalSet().is_superset_of(EmptySet())
    assert UniversalSet().is_superset_of(ValueSet(1))

    # EmptySet is contained in all sets
    assert EmptySet().is_subset_of(EmptySet())
    assert EmptySet().is_superset_of(EmptySet())
    assert EmptySet().is_subset_of(ValueSet(2))

    # PropertySets <-> PropertySets
    assert ValueSet(2).is_superset_of(ValueSet(2))
    assert ValueSet(2).is_subset_of(ValueSet(2))
    assert not ValueSet(2).is_subset_of(ValueSet(3))
    assert not ValueSet(2).is_superset_of(ValueSet(3))

    # PropertySet <-> UniversalSet
    assert ValueSet(1).is_subset_of(UniversalSet())
    assert not UniversalSet().is_subset_of(ValueSet(1))
    assert UniversalSet().is_superset_of(ValueSet(1))
    assert not ValueSet(1).is_superset_of(UniversalSet())

    # PropertySet <-> EmptySet
    assert ValueSet(2).is_superset_of(EmptySet())
    assert not EmptySet().is_superset_of(ValueSet(2))
    assert EmptySet().is_subset_of(ValueSet(1))
    assert not ValueSet(1).is_subset_of(EmptySet())


def test_superset_subset_for_flat_intersections() -> None:
    assert Intersection(ValueSet(1), ValueSet(2)).is_subset_of(ValueSet(1))
    assert Intersection(ValueSet(1), ValueSet(2)).is_subset_of(ValueSet(2))

    assert not Intersection(ValueSet(1), ValueSet(2)).is_equivalent_to(ValueSet(1))
    assert not Intersection(ValueSet(1), ValueSet(2)).is_equivalent_to(ValueSet(2))

    assert ValueSet(1).is_superset_of(Intersection(ValueSet(1), ValueSet(2)))
    assert ValueSet(2).is_superset_of(Intersection(ValueSet(1), ValueSet(2)))


def test_superset_subset_for_nested_intersections() -> None:
    x = ValueSet(1)
    y = ValueSet(2)
    z = ValueSet(3)
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

    assert not xyz.is_superset_of(x)
    assert not xyz.is_superset_of(y)
    assert not xyz.is_superset_of(z)
    assert not xyz.is_superset_of(xy)
    assert not xyz.is_superset_of(yz)
    assert not xyz.is_superset_of(xz)

    assert not x.is_subset_of(xyz)
    assert not y.is_subset_of(xyz)
    assert not z.is_subset_of(xyz)
    assert not xy.is_subset_of(xyz)
    assert not yz.is_subset_of(xyz)
    assert not xz.is_subset_of(xyz)


def test_superset_subset_for_flat_unions() -> None:
    assert Union(ValueSet(1), ValueSet(2)).is_superset_of(ValueSet(1))
    assert Union(ValueSet(1), ValueSet(2)).is_superset_of(ValueSet(2))

    assert not Union(ValueSet(1), ValueSet(2)).is_equivalent_to(ValueSet(1))
    assert not Union(ValueSet(1), ValueSet(2)).is_equivalent_to(ValueSet(2))

    assert ValueSet(1).is_subset_of(Union(ValueSet(1), ValueSet(2)))
    assert ValueSet(2).is_subset_of(Union(ValueSet(1), ValueSet(2)))


def test_superset_subset_for_nested_unions() -> None:
    x = ValueSet(1)
    y = ValueSet(2)
    z = ValueSet(3)
    xy = Union(x, y)
    yz = Union(y, z)
    xz = Union(x, z)
    xyz = Union(x, Union(y, z))

    for a in [x, y, z, xy, yz, xz, xyz]:
        for b in [x, y, z, xy, yz, xz, xyz]:
            if a != b:
                assert not a.is_equivalent_to(b)

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
def test_de_morgan_s_identities(sets: List[Set]) -> None:
    for x, y in itt.product(sets, sets):
        assert Intersection(Complement(x), Complement(y)).is_equivalent_to(Complement(Union(x, y)))
        assert Union(Complement(x), Complement(y)).is_equivalent_to(Complement(Intersection(x, y)))


def test_complement_squares_to_no_op(sets: List[Set]) -> None:
    for x in sets:
        assert x.is_equivalent_to(Complement(Complement(x)))


def test_intersection_properties(sets: List[Set]) -> None:
    for x in sets:
        assert EmptySet().is_equivalent_to(Intersection(x, Complement(x)))
        assert EmptySet().is_equivalent_to(Intersection(Complement(x), x))

        assert EmptySet().is_equivalent_to(Intersection(EmptySet(), x))
        assert EmptySet().is_equivalent_to(Intersection(x, EmptySet()))

        assert x.is_equivalent_to(Intersection(UniversalSet(), x))
        assert x.is_equivalent_to(Intersection(x, UniversalSet()))

        assert x.is_equivalent_to(Intersection(x, x))


def test_union_properties(sets: List[Set]) -> None:
    for x in sets:
        assert UniversalSet().is_equivalent_to(Union(x, Complement(x)))
        assert UniversalSet().is_equivalent_to(Union(Complement(x), x))

        assert Union(x, Complement(x)).is_equivalent_to(UniversalSet())
        assert Union(Complement(x), x).is_equivalent_to(UniversalSet())

        assert x.is_equivalent_to(Union(EmptySet(), x))
        assert x.is_equivalent_to(Union(x, EmptySet()))

        assert UniversalSet().is_equivalent_to(Union(UniversalSet(), x))
        assert UniversalSet().is_equivalent_to(Union(x, UniversalSet()))

        assert x.is_equivalent_to(Union(x, x))


def test_absolute_relative_complement_identities(sets: List[Set]) -> None:
    for x, y in itt.product(sets, sets):
        assert Intersection(x, Complement(y)).is_equivalent_to(Difference(x, y))
        assert Union(Complement(x), y).is_equivalent_to(Complement(Difference(x, y)))


def test_distributive_properties(sets: List[Set]) -> None:
    for x, y, z in itt.product(sets, sets, sets):
        assert Union(x, Intersection(y, z)).is_equivalent_to(Intersection(Union(x, y), Union(x, z)))
        assert Intersection(x, Union(y, z)).is_equivalent_to(Union(Intersection(x, y), Intersection(x, z)))


def test_absorption_properties(sets: List[Set]) -> None:
    for x, y in itt.product(sets, sets):
        assert x.is_equivalent_to(Union(x, Intersection(x, y)))
        assert x.is_equivalent_to(Intersection(x, Union(x, y)))


def test_disjunctive_union(sets: List[Set]) -> None:
    for x, y in itt.product(sets, sets):
        assert DisjunctiveUnion(x, y).is_equivalent_to(Union(Difference(x, y), Difference(y, x)))
        assert Union(Difference(x, y), Difference(y, x)).is_equivalent_to(DisjunctiveUnion(x, y))

        assert DisjunctiveUnion(x, y).is_equivalent_to(Difference(Union(x, y), Intersection(x, y)))
        assert Difference(Union(x, y), Intersection(x, y)).is_equivalent_to(DisjunctiveUnion(x, y))

def test_is_equivalent_to() -> None:
    assert UniversalSet().is_equivalent_to(UniversalSet())
    assert EmptySet().is_equivalent_to(EmptySet())

    assert not UniversalSet().is_equivalent_to(ValueSet(1))
    assert not ValueSet(1).is_equivalent_to(UniversalSet())

    assert UniversalSet().is_equivalent_to(Union(ValueSet(1), Complement(ValueSet(1))))


# # Test the cardinality calculator
# def test_cardinality_empty_and_universal_set():
#     assert EmptySet().cardinality == {}
#     assert UniversalSet().cardinality == {(UniversalSet(),): 1}
#
#
# def test_cardinality_respects_complement_properties(sets):
#     for x in sets:
#         assert Intersection(x, Complement(x)).cardinality == {}
#         assert Union(x, Complement(x)).cardinality == {(UniversalSet(),): 1}
#
#
# def test_cardinality_property_sets():
#     # Test the standard examples of the inclusion-exclusion principle:
#     # |A u (B u C)| = |A| + |B| + |C| - |A^B| - |A^C| - |B^C| + |A^B^C|
#     assert Union(ValueSet(1), Union(ValueSet(2), ValueSet(3))).cardinality == \
#            {(ValueSet(1),): 1, (ValueSet(2),): 1, (ValueSet(3),): 1, (ValueSet(1), ValueSet(2)): -1, (ValueSet(1), ValueSet(3)): -1,
#             (ValueSet(2), ValueSet(3)): -1, (ValueSet(1), ValueSet(2), ValueSet(3)): 1}
#
#     # |(A u B) u C| = same
#     assert Union(Union(ValueSet(1), ValueSet(2)), ValueSet(3)).cardinality == \
#            {(ValueSet(1),): 1, (ValueSet(2),): 1, (ValueSet(3),): 1, (ValueSet(1), ValueSet(2)): -1, (ValueSet(1), ValueSet(3)): -1,
#             (ValueSet(2), ValueSet(3)): -1, (ValueSet(1), ValueSet(2), ValueSet(3)): 1}
#
#     # |A u (B u (C u D))| = |A| + |B| + |C| + |D| - |A^B| - |A^C| - |A^D| - |B^C| - |B^D| -|C^D| + |A^B^C| +
#     #                       |B^C^D| + |A^C^D| + |A^B^D| - |A^B^C^D|
#     assert Union(ValueSet(1), Union(ValueSet(2), Union(ValueSet(3), ValueSet(4)))).cardinality == \
#            {(ValueSet(1),): 1, (ValueSet(2),): 1, (ValueSet(3),): 1, (ValueSet(4),): 1,
#             (ValueSet(1), ValueSet(2)): -1, (ValueSet(1), ValueSet(3)): -1, (ValueSet(1), ValueSet(4)): -1,
#             (ValueSet(2), ValueSet(3)): -1, (ValueSet(2), ValueSet(4)): -1, (ValueSet(3), ValueSet(4)): -1,
#             (ValueSet(1), ValueSet(2), ValueSet(3)): 1, (ValueSet(2), ValueSet(3), ValueSet(4)): 1,
#             (ValueSet(1), ValueSet(3), ValueSet(4)): 1, (ValueSet(1), ValueSet(2), ValueSet(4)): 1,
#             (ValueSet(1), ValueSet(2), ValueSet(3), ValueSet(4)): -1}
#
#     # |(A u B) u (C u D)| = same
#     assert Union(Union(ValueSet(1), ValueSet(2)), Union(ValueSet(3), ValueSet(4))).cardinality == \
#            {(ValueSet(1),): 1, (ValueSet(2),): 1, (ValueSet(3),): 1, (ValueSet(4),): 1,
#             (ValueSet(1), ValueSet(2)): -1, (ValueSet(1), ValueSet(3)): -1, (ValueSet(1), ValueSet(4)): -1,
#             (ValueSet(2), ValueSet(3)): -1, (ValueSet(2), ValueSet(4)): -1, (ValueSet(3), ValueSet(4)): -1,
#             (ValueSet(1), ValueSet(2), ValueSet(3)): 1, (ValueSet(2), ValueSet(3), ValueSet(4)): 1,
#             (ValueSet(1), ValueSet(3), ValueSet(4)): 1, (ValueSet(1), ValueSet(2), ValueSet(4)): 1,
#             (ValueSet(1), ValueSet(2), ValueSet(3), ValueSet(4)): -1}
#
#
# def test_cardinality_property_sets_and_complements():
#     assert Union(Complement(ValueSet(1)), ValueSet(2)).cardinality == \
#            {(UniversalSet(),): 1, (ValueSet(1),): -1, (ValueSet(1), ValueSet(2)): 1}
#
#     assert Union(Complement(ValueSet(1)), Union(ValueSet(2), Complement(ValueSet(3)))).cardinality == \
#            {(UniversalSet(),): 1, (ValueSet(1), ValueSet(3)): -1, (ValueSet(1), ValueSet(2), ValueSet(3)): 1}
#
#     assert Intersection(Complement(ValueSet(1)), Intersection(ValueSet(2), Complement(ValueSet(3)))).cardinality == \
#            {(ValueSet(2),): 1, (ValueSet(1), ValueSet(2)): -1, (ValueSet(2), ValueSet(3)): -1,
#             (ValueSet(1), ValueSet(2), ValueSet(3)): 1}
#



@pytest.fixture
def sets() -> List[Set]:
    return [
        EmptySet(),
        ValueSet(1),
        UniversalSet(),
        Union(ValueSet(1), ValueSet(2)),
        Intersection(ValueSet(1), ValueSet(2)),
        Intersection(ValueSet(1), Complement(ValueSet(2))),
        Union(Intersection(ValueSet(1), ValueSet(2)), ValueSet(3)),
        Union(Intersection(ValueSet(1), ValueSet(2)), Intersection(ValueSet(3), ValueSet(4))),
        Union(Complement(Union(ValueSet(1), Complement(ValueSet(2)))), Intersection(ValueSet(3), ValueSet(4)))
    ]
