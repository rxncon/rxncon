import pytest
import itertools as itt
from typing import List

from rxncon.venntastic.sets import *

# # Test the basic properties of the Set datastructure
# def test_property_set_construction():
#     assert PropertySet('a')
#     assert PropertySet(1)
#     assert PropertySet(1.45)
#
#     with pytest.raises(TypeError):
#         PropertySet([1, 2, 3])
#
#
# def test_property_set_comparison():
#     assert PropertySet(1) == PropertySet(1)
#     assert not PropertySet('a') == PropertySet('b')
#
#
# def test_property_set_dictionary_keys():
#     x = PropertySet(1)
#     y = PropertySet(2)
#
#     dictionary = {x: 'bla', y: 'diebla'}
#
#     assert dictionary[PropertySet(1)] == 'bla'
#     assert dictionary[x] == 'bla'
#     assert dictionary[PropertySet(2)] == 'diebla'
#     assert dictionary[y] == 'diebla'
#
#
# def test_property_set_ordering():
#     assert PropertySet(1) < PropertySet(2)
#     assert PropertySet('a') < PropertySet('b')
#
#
# def test_complement_ordering():
#     assert Complement(PropertySet(1)) < Complement(PropertySet(2))
#     assert Complement(PropertySet('a')) < Complement(PropertySet('b'))
#
#
# def test_intersection_ordering():
#     assert Intersection(PropertySet(1), PropertySet(2)) < Intersection(PropertySet(2), PropertySet(3))
#     assert Intersection(PropertySet(1), PropertySet(2)) < Intersection(Intersection(PropertySet(1), PropertySet(2)), PropertySet(3))
#
#
# def test_union_ordering():
#     assert Union(PropertySet(1), PropertySet(2)) < Union(PropertySet(2), PropertySet(3))
#     assert Union(PropertySet(1), PropertySet(2)) < Union(Union(PropertySet(1), PropertySet(2)), PropertySet(3))
#
#
# def test_mutual_ordering():
#     assert EmptySet() < UniversalSet()
#     assert UniversalSet() < PropertySet(1)
#     assert PropertySet(1) < Complement(PropertySet(1))
#     assert Complement(PropertySet(1)) < Intersection(PropertySet(1), PropertySet(2))
#     assert Intersection(PropertySet(1), PropertySet(2)) < Union(PropertySet(1), PropertySet(2))
#
#
# # Test the superset / subset relationships
# def test_superset_subset_for_unary_sets():
#     assert UniversalSet() == PropertySet(None)
#
#     assert UniversalSet().is_superset_of(UniversalSet())
#     assert UniversalSet().is_subset_of(UniversalSet())
#
#     assert UniversalSet().is_superset_of(EmptySet())
#     assert UniversalSet().is_superset_of(PropertySet(1))
#
#     assert PropertySet(2).is_superset_of(PropertySet(2))
#     assert PropertySet(2).is_subset_of(PropertySet(2))
#
#     assert PropertySet(2).is_subset_of(UniversalSet())
#     assert PropertySet(2).is_superset_of(EmptySet())
#     assert EmptySet().is_subset_of(PropertySet(2))
#
#     assert EmptySet().is_subset_of(EmptySet())
#     assert EmptySet().is_superset_of(EmptySet())
#
#     assert not PropertySet(2).is_subset_of(PropertySet(3))
#     assert not PropertySet(2).is_superset_of(PropertySet(3))
#
#
# def test_superset_subset_for_flat_intersections():
#     assert Intersection(PropertySet(1), PropertySet(2)).is_subset_of(PropertySet(1))
#     assert Intersection(PropertySet(1), PropertySet(2)).is_subset_of(PropertySet(2))
#
#     assert not Intersection(PropertySet(1), PropertySet(2)).is_equivalent_to(PropertySet(1))
#     assert not Intersection(PropertySet(1), PropertySet(2)).is_equivalent_to(PropertySet(2))
#
#     assert PropertySet(1).is_superset_of(Intersection(PropertySet(1), PropertySet(2)))
#     assert PropertySet(2).is_superset_of(Intersection(PropertySet(1), PropertySet(2)))
#
#
# def test_superset_subset_for_nested_intersections():
#     x = PropertySet(1)
#     y = PropertySet(2)
#     z = PropertySet(3)
#     xy = Intersection(x, y)
#     yz = Intersection(y, z)
#     xz = Intersection(x, z)
#     xyz = Intersection(x, Intersection(y, z))
#
#     assert xyz.is_subset_of(x)
#     assert xyz.is_subset_of(y)
#     assert xyz.is_subset_of(z)
#     assert xyz.is_subset_of(xy)
#     assert xyz.is_subset_of(yz)
#     assert xyz.is_subset_of(xz)
#
#     assert x.is_superset_of(xyz)
#     assert y.is_superset_of(xyz)
#     assert z.is_superset_of(xyz)
#     assert xy.is_superset_of(xyz)
#     assert yz.is_superset_of(xyz)
#     assert xz.is_superset_of(xyz)
#
#
# def test_superset_subset_for_flat_unions():
#     assert Union(PropertySet(1), PropertySet(2)).is_superset_of(PropertySet(1))
#     assert Union(PropertySet(1), PropertySet(2)).is_superset_of(PropertySet(2))
#
#     assert not Union(PropertySet(1), PropertySet(2)).is_equivalent_to(PropertySet(1))
#     assert not Union(PropertySet(1), PropertySet(2)).is_equivalent_to(PropertySet(2))
#
#     assert PropertySet(1).is_subset_of(Union(PropertySet(1), PropertySet(2)))
#     assert PropertySet(2).is_subset_of(Union(PropertySet(1), PropertySet(2)))
#
#
# def test_superset_subset_for_nested_unions():
#     x = PropertySet(1)
#     y = PropertySet(2)
#     z = PropertySet(3)
#     xy = Union(x, y)
#     yz = Union(y, z)
#     xz = Union(x, z)
#     xyz = Union(x, Union(y, z))
#
#     assert xyz.is_superset_of(x)
#     assert xyz.is_superset_of(y)
#     assert xyz.is_superset_of(z)
#     assert xyz.is_superset_of(xy)
#     assert xyz.is_superset_of(yz)
#     assert xyz.is_superset_of(xz)
#
#     assert x.is_subset_of(xyz)
#     assert y.is_subset_of(xyz)
#     assert z.is_subset_of(xyz)
#     assert xy.is_subset_of(xyz)
#     assert yz.is_subset_of(xyz)
#     assert xz.is_subset_of(xyz)
#
#
#
# # Test basic set-theoretic identities for a generated pool of set expressions
# def test_de_morgan_s_identities(sets):
#     for x, y in itt.product(sets, sets):
#         assert Intersection(Complement(x), Complement(y)).is_equivalent_to(Complement(Union(x, y)))
#         assert Union(Complement(x), Complement(y)).is_equivalent_to(Complement(Intersection(x, y)))
#
#
# def test_complement_squares_to_no_op(sets):
#     for x in sets:
#         assert x.is_equivalent_to(Complement(Complement(x)))
#
#
# def test_intersection_properties(sets):
#     for x in sets:
#         assert EmptySet().is_equivalent_to(Intersection(x, Complement(x)))
#         assert EmptySet().is_equivalent_to(Intersection(Complement(x), x))
#
#         assert EmptySet().is_equivalent_to(Intersection(EmptySet(), x))
#         assert EmptySet().is_equivalent_to(Intersection(x, EmptySet()))
#
#         assert x.is_equivalent_to(Intersection(UniversalSet(), x))
#         assert x.is_equivalent_to(Intersection(x, UniversalSet()))
#
#         assert x.is_equivalent_to(Intersection(x, x))
#

def test_union_properties(sets):
    for x in sets:
        assert UniversalSet().is_equivalent_to(Union(x, Complement(x)))
        # assert UniversalSet().is_equivalent_to(Union(Complement(x), x))
        #
        # assert x.is_equivalent_to(Union(EmptySet(), x))
        # assert x.is_equivalent_to(Union(x, EmptySet()))
        #
        # assert UniversalSet().is_equivalent_to(Union(UniversalSet(), x))
        # assert UniversalSet().is_equivalent_to(Union(x, UniversalSet()))
        #
        # assert x.is_equivalent_to(Union(x, x))


# def test_absolute_relative_complement_identities(sets):
#     for x, y in itt.product(sets, sets):
#         assert Intersection(x, Complement(y)).is_equivalent_to(Difference(x, y))
#         assert Union(Complement(x), y).is_equivalent_to(Complement(Difference(x, y)))
#
#
# def test_distributive_properties(sets):
#     for x, y, z in itt.product(sets, sets, sets):
#         assert Union(x, Intersection(y, z)).is_equivalent_to(Intersection(Union(x, y), Union(x, z)))
#         assert Intersection(x, Union(y, z)).is_equivalent_to(Union(Intersection(x, y), Intersection(x, z)))
#
#
# def test_absorption_properties(sets):
#     for x, y in itt.product(sets, sets):
#         assert x.is_equivalent_to(Union(x, Intersection(x, y)))
#         assert x.is_equivalent_to(Intersection(x, Union(x, y)))


@pytest.fixture
def sets():
    sets = [
        # EmptySet(),
        # PropertySet(1),
        # UniversalSet(),
        # Union(PropertySet(1), PropertySet(2)),
        # Intersection(PropertySet(1), PropertySet(2)),
        # Union(Intersection(PropertySet(1), PropertySet(2)), PropertySet(3)),
        Union(Intersection(PropertySet(1), PropertySet(2)), Intersection(PropertySet(3), PropertySet(4)))
    ]



    return sets









