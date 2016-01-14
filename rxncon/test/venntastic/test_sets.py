# import pytest
# import itertools as itt
#
# from rxncon.venntastic.sets import *
#
#

#
# def test_equivalence_complement_empty_and_universal_set():
#     # assert Complement(EmptySet()).is_equivalent_to(UniversalSet())
#     assert Complement(UniversalSet()).is_equivalent_to(EmptySet())
#
#

#


# # #
# # # def test_cardinality_empty_and_universal_set():
# # #     assert EmptySet().cardinality() == {}
# # #
# # #     assert PropertySet().cardinality() == {PropertySet(): 1}
# # #
# # #

# # #
# # # def test_cardinality_inclusion_exclusion_principle():
# # #     # Test the standard examples of the inclusion-exclusion principle:
# # #     # |A u (B u C)| = |A| + |B| + |C| - |A^B| - |A^C| - |B^C| + |A^B^C|
# # #     assert Union(PropertySet(1), Union(PropertySet(2), PropertySet(3))).cardinality() == \
# # #         {PropertySet(1): 1, PropertySet(2): 1, PropertySet(3): 1, PropertySet(1, 2): -1, PropertySet(1, 3): -1,
# # #          PropertySet(2, 3): -1, PropertySet(1, 2, 3): 1}
# # #
# # #     # |(A u B) u C| = same
# # #     assert Union(Union(PropertySet(1), PropertySet(2)), PropertySet(3)).cardinality() == \
# # #         {PropertySet(1): 1, PropertySet(2): 1, PropertySet(3): 1, PropertySet(1, 2): -1, PropertySet(1, 3): -1,
# # #          PropertySet(2, 3): -1, PropertySet(1, 2, 3): 1}
# # #
# # #     # |A u (B u (C u D))| = |A| + |B| + |C| + |D| - |A^B| - |A^C| - |A^D| - |B^C| - |B^D| -|C^D| + |A^B^C| +
# # #     #                       |B^C^D| + |A^C^D| + |A^B^D| - |A^B^C^D|
# # #     assert Union(PropertySet(1), Union(PropertySet(2), Union(PropertySet(3), PropertySet(4)))).cardinality() == \
# # #         {PropertySet(1): 1, PropertySet(2): 1, PropertySet(3): 1, PropertySet(4): 1,
# # #          PropertySet(1, 2): -1, PropertySet(1, 3): -1, PropertySet(1, 4): -1, PropertySet(2, 3): -1,
# # #          PropertySet(2, 4): -1, PropertySet(3, 4): -1,
# # #          PropertySet(1, 2, 3): 1, PropertySet(2, 3, 4): 1, PropertySet(1, 3, 4): 1, PropertySet(1, 2, 4): 1,
# # #          PropertySet(1, 2, 3, 4): -1}
# # #
# # #     # |(A u B) u (C u D)| = same
# # #     assert Union(Union(PropertySet(1), PropertySet(2)), Union(PropertySet(3), PropertySet(4))).cardinality() == \
# # #         {PropertySet(1): 1, PropertySet(2): 1, PropertySet(3): 1, PropertySet(4): 1,
# # #          PropertySet(1, 2): -1, PropertySet(1, 3): -1, PropertySet(1, 4): -1, PropertySet(2, 3): -1,
# # #          PropertySet(2, 4): -1, PropertySet(3, 4): -1,
# # #          PropertySet(1, 2, 3): 1, PropertySet(2, 3, 4): 1, PropertySet(1, 3, 4): 1, PropertySet(1, 2, 4): 1,
# # #          PropertySet(1, 2, 3, 4): -1}
# # #
# # #
# # # def test_cardinality_of_expressions_containing_complements():
# # #     assert Union(Complement(PropertySet(1)), PropertySet(2)).cardinality() == \
# # #         {PropertySet(): 1, PropertySet(1): -1, PropertySet(1, 2): 1}
# # #
# # #     assert Union(Complement(PropertySet(1)), Union(PropertySet(2), Complement(PropertySet(3)))).cardinality() == \
# # #         {PropertySet(): 1, PropertySet(1, 3): -1, PropertySet(1, 2, 3): 1}
# # #
# # #     assert Intersection(Complement(PropertySet(1)), Intersection(PropertySet(2), Complement(PropertySet(3)))).cardinality() == \
# # #         {PropertySet(2): 1, PropertySet(1, 2): -1, PropertySet(2, 3): -1, PropertySet(1, 2, 3): 1}
# # #
