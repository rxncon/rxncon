import pytest
from rxncon.util.utils import OrderedEnum
import rxncon.core.reaction as rxn


class Feelings(OrderedEnum):
    happy = 1
    sad = 2
    tired = 3
    mixed = None

class Food(OrderedEnum):
    kebap = 'kebap'
    burger = 'burger'
    pizza = 'pizza'

@pytest.fixture
def to_sort():
    return [
        [Feelings.tired, Feelings.sad, Feelings.happy],
        [Feelings.tired, Feelings.sad, Feelings.happy, Feelings.mixed],
        [Food.pizza, Food.kebap, Food.burger],
        [rxn.CovalentReactionModifier.ubiquitin, rxn.CovalentReactionModifier.phosphor,
         rxn.CovalentReactionModifier.guanosintriphosphat, rxn.CovalentReactionModifier.undefined,
         rxn.CovalentReactionModifier.truncated]
    ]

@pytest.fixture
def expected_sorting():
    return [
        [Feelings.happy, Feelings.sad, Feelings.tired],
        [Feelings.mixed, Feelings.happy, Feelings.sad, Feelings.tired],
        [Food.burger, Food.kebap, Food.pizza],
        [rxn.CovalentReactionModifier.undefined, rxn.CovalentReactionModifier.guanosintriphosphat,
         rxn.CovalentReactionModifier.phosphor, rxn.CovalentReactionModifier.truncated,
         rxn.CovalentReactionModifier.ubiquitin]
    ]

def test_OrderedEnum(to_sort, expected_sorting):
    for i, element in enumerate(to_sort):
        assert sorted(element) == expected_sorting[i]
