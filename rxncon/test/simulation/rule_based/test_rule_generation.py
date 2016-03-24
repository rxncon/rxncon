from collections import namedtuple
from typing import List
import pytest
from rxncon.input.quick.quick import Quick
from rxncon.simulation.rule_based.molecule_from_string import mol_def_from_string, rule_from_string
from rxncon.simulation.rule_based.rbm_from_rxncon import RuleBasedModelSupervisor
from rxncon.simulation.rule_based.rule_based_model import Rule


RuleTestCase = namedtuple('RuleTestCase', ['quick_string', 'mol_def_strings', 'rule_strings'])


def test_rule_generation(test_cases):
    for test_case in test_cases:
        assert is_rule_test_case_correct(test_case)


@pytest.fixture
def test_cases():
    return [
        RuleTestCase(
            'A_ppi_B',
            ['A#ass/A_[Bassoc]:B_[Aassoc]', 'B#ass/B_[Aassoc]:A_[Bassoc'],
            ['A#ass/A_[Bassoc]: + B#ass/B_[Aassoc]: <-> A#ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0']
        ),
        RuleTestCase(
            '''
            A_ppi_B; ! A-{p}
            C_p+_A''',
            ['A#ass/A_[Bassoc]:B_[Aassoc],mod/A_[(Csite)]:u~p', 'B#ass/B_[Aassoc]:A_[Bassoc]', 'C#'],
            ['A#ass/A_[Bassoc]:,mod/A_[(Csite)]:p + B#ass/B_[Aassoc]: <-> A#ass/A_[Bassoc]:B_[Aassoc]~0,mod/A_[(Csite)]:p.B#ass/B_[Aassoc]:A_[Bassoc]~0',
             'C# + A#mod/A_[(Csite)]:u <-> C# + A#mod/A_[(Csite)]:p']
        )
    ]


def is_rule_test_case_correct(test_case) -> bool:
    rxncon = Quick(test_case.quick_string).rxncon_system
    actual_mol_defs = set(RuleBasedModelSupervisor(rxncon).mol_defs.values())
    actual_rules    = RuleBasedModelSupervisor(rxncon).rules

    expected_mol_defs = {mol_def_from_string(x) for x in test_case.mol_def_strings}
    expected_rules    = [rule_from_string(expected_mol_defs, x) for x in test_case.rule_strings]

    return actual_mol_defs == expected_mol_defs and are_rule_lists_equivalent(actual_rules, expected_rules)


def are_rule_lists_equivalent(first_list: List[Rule], second_list: List[Rule]) -> bool:
    while first_list:
        first_rule = first_list.pop()
        for second_rule in second_list:
            if are_rules_equivalent(first_rule, second_rule):
                second_list.remove(second_rule)

    return len(second_list) == 0


def are_rules_equivalent(first_rule: Rule, second_rule: Rule) -> bool:
    # @todo we disregard the rates in this equivalence.
    return set(first_rule.left_hand_side) == set(second_rule.left_hand_side) and \
        set(first_rule.right_hand_side) == set(second_rule.right_hand_side) and \
        first_rule.arrow_type == second_rule.arrow_type