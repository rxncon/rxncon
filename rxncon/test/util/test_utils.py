import pytest
from collections import namedtuple

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

EnumSortingTestCase = namedtuple('EnumSortingTestCase',  ['to_sort', 'expected_sorting'])
@pytest.fixture
def the_case_sorting():
    return [
        EnumSortingTestCase([Feelings.tired, Feelings.sad, Feelings.happy],
                            [Feelings.happy, Feelings.sad, Feelings.tired],),

        EnumSortingTestCase([Feelings.tired, Feelings.sad, Feelings.happy, Feelings.mixed],
                            [Feelings.mixed, Feelings.happy, Feelings.sad, Feelings.tired],),

        EnumSortingTestCase([Food.pizza, Food.kebap, Food.burger],
                            [Food.burger, Food.kebap, Food.pizza],),

        EnumSortingTestCase([rxn.CovalentReactionModifier.ubiquitin, rxn.CovalentReactionModifier.phosphor,
                            rxn.CovalentReactionModifier.guanosintriphosphat, rxn.CovalentReactionModifier.undefined,
                            rxn.CovalentReactionModifier.truncated],

                            [rxn.CovalentReactionModifier.undefined, rxn.CovalentReactionModifier.guanosintriphosphat,
                             rxn.CovalentReactionModifier.phosphor, rxn.CovalentReactionModifier.truncated,
                             rxn.CovalentReactionModifier.ubiquitin]
                            )
    ]

def test_OrderedEnum(the_case_sorting):
    for the_case in the_case_sorting:
        assert sorted(the_case.to_sort) == the_case.expected_sorting
