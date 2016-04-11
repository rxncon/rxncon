from rxncon.simulation.rule_based.rbm_from_rxncon import rules_from_reaction, rbm_from_rxncon_sys
from rxncon.input.quick.quick import Quick
from rxncon.syntax.rxncon_from_string import reaction_from_string
from rxncon.simulation.rule_based.molecule_from_string import rule_from_string, mol_def_from_string


def test_single_reaction():
    rxncon_sys = Quick('A_ppi_B').rxncon_system

    rules = rules_from_reaction(rxncon_sys, reaction_from_string('A_ppi_B'))
    assert len(rules) == 1

    the_rule = list(rules)[0]

    assert the_rule == rule_from_string(
        ['A#ass/A_[Bassoc]:B_[Aassoc]', 'B#ass/B_[Aassoc]:A_[Bassoc]'],
        'A#ass/A_[Bassoc]: + B#ass/B_[Aassoc]: <-> A#ass/A_[Bassoc]:B_[Aassoc].B#ass/B_[Aassoc]:A_[Bassoc] @ kf_A_ppi_B, kr_A_ppi_B'
    )