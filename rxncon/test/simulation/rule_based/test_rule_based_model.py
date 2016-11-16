from rxncon.input.quick.quick import Quick
from rxncon.simulation.rule_based.rule_based_model import rule_based_model_from_rxncon


def test_simple_system():
    rbm = rule_based_model_from_rxncon(Quick("""A_[b]_ppi+_B_[a]; ! A_[(r)]-{p}
                                             A_[b]_ppi-_B_[a]
                                             C_p+_A_[(r)]
                                             D_p-_A_[(r)]""").rxncon_system)

    print()
    for rule in rbm.rules:
        print(rule)


def test_2():
    rbm = rule_based_model_from_rxncon(Quick("""E_p+_C
                                             D_p+_C
                                             A_ppi+_C
                                             C_p+_B_[(r)]; ! C-{p}
                                             C_ub+_B_[(r)]; ! A--C; ! C-{P}""").rxncon_system)

    print()
    for rule in rbm.rules:
        print(rule)

