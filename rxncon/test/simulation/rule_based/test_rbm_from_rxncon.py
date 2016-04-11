from rxncon.simulation.rule_based.rbm_from_rxncon import rules_from_reaction, rbm_from_rxncon_sys
from rxncon.input.quick.quick import Quick
from rxncon.syntax.rxncon_from_string import reaction_from_string

def test_single_reaction():
    rxncon_sys = Quick('A_ppi_B').rxncon_system

    rules = rules_from_reaction(rxncon_sys, reaction_from_string('A_ppi_B'))

    for rule in rules:
        print(rule)