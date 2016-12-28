import pytest
from collections import namedtuple

from rxncon.util.utils import OrderedEnum


class Feelings(OrderedEnum):
    happy = 1
    sad   = 2
    tired = 3
    mixed = None

class Food(OrderedEnum):
    kebap  = 'kebap'
    burger = 'burger'
    pizza  = 'pizza'

EnumSortingTestCase = namedtuple('EnumSortingTestCase',  ['to_sort', 'expected_sorting'])
@pytest.fixture
def sorting_cases():
    return [
        EnumSortingTestCase([Feelings.tired, Feelings.sad, Feelings.happy],
                            [Feelings.happy, Feelings.sad, Feelings.tired]),

        EnumSortingTestCase([Feelings.tired, Feelings.sad, Feelings.happy, Feelings.mixed],
                            [Feelings.mixed, Feelings.happy, Feelings.sad, Feelings.tired]),

        EnumSortingTestCase([Food.pizza, Food.kebap, Food.burger],
                            [Food.burger, Food.kebap, Food.pizza]),
    ]

def test_OrderedEnum(sorting_cases):
    for case in sorting_cases:
        assert sorted(case.to_sort) == case.expected_sorting
