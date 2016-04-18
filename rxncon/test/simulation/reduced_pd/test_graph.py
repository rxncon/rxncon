import pytest
from collections import namedtuple
import rxncon.simulation.reducedPD.graph as gra
from rxncon.simulation.rule_based.molecule_from_string import mol_def_from_string, rule_from_string

RuleTestCase = namedtuple('RuleTestCase', ['mol_def_strings', 'rule_strings', 'graph'])
@pytest.fixture
def simple_test_cases():
    return [
        RuleTestCase(
            ['A#', 'B#mod/B_[(Asite)]:u~p'],
            ['A# + B#mod/B_[(Asite)]:u -> A# + B#mod/B_[(Asite)]:p @ k_A_p+_B'],
            []
        )
    ]

def test_rule_generation(the_cases):
    for test_case in the_cases:
        assert is_graph_test_case_correct(test_case)

@pytest.fixture
def the_cases(simple_test_cases):
    return simple_test_cases

def is_graph_test_case_correct(test_case):
    g = gra.Graph()

    expected_mol_defs = {mol_def_from_string(x) for x in test_case.mol_def_strings}
    expected_rules    = {rule_from_string(expected_mol_defs, x) for x in test_case.rule_strings}
    rule = expected_rules[0]
    g.add_node(rule.left_hand_side)

