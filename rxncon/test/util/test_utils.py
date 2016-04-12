import pytest
from rxncon.util.utils import OrderedEnum


class Feelings(OrderedEnum):
    happy = 1
    sad = 2
    tired = 3

class Food(OrderedEnum):
    kebap = 'kebap'
    burger = 'burger'
    pizza = 'pizza'

@pytest.fixture
def to_sort():
    return [
        [Feelings.tired, Feelings.sad, Feelings.happy],
        [Food.pizza, Food.kebap, Food.burger]
    ]

@pytest.fixture
def expected_sorting():
    return [
        [Feelings.happy, Feelings.sad, Feelings.tired],
        [Food.burger, Food.kebap, Food.pizza]
    ]

def test_OrderedEnum(to_sort, expected_sorting):
    for i, element in enumerate(to_sort):
        assert sorted(element) == expected_sorting[i]
