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


def test_calc_positive_solutions():
    rxncon_sys = Quick("""A_[b]_ppi+_B_[a]; x A_[(r)]-{p}
                          E_[x]_ppi+_B_[a]
                          C_p+_A_[(r)]
                          D_ub+_A_[(r)]""").rxncon_system

    for x in rule_based_model_from_rxncon(rxncon_sys).rules:
        print(x)



def test_2():
    rbm = rule_based_model_from_rxncon(Quick('''A_ppi+_C
                                             C_ppi+_D
                                             B_ppi+_E
                                             B_ppi+_F
                                             A_ppi+_B; x <comp1>
                                             <comp1>; OR <comp1C1>
                                             <comp1>; OR <comp2C1>
                                             <comp1C1>; AND A--C
                                             <comp1C1>; AND C--D
                                             <comp2C1>; AND B--F
                                             <comp2C1>; AND B--E''').rxncon_system)

    print()
    for rule in rbm.rules:
        print(rule)

def test_3():
    rbm = rule_based_model_from_rxncon(Quick('''A_ppi+_C
                                             C_ppi+_D
                                             B_ppi+_E
                                             A_ppi+_B; x <comp1>
                                             <comp1>; AND <comp1C1>
                                             <comp1>; AND <comp2C1>
                                             <comp1C1>; OR A--C
                                             <comp1C1>; OR C--D
                                             <comp2C1>; OR A--C
                                             <comp2C1>; OR B--E''').rxncon_system)

    print()
    for rule in rbm.rules:
        print(rule)


def test_mutually_exclusive_bindings():
    rbm = rule_based_model_from_rxncon(Quick('''C_[A]_ppi+_A_[x]
                                             D_[A]_ppi+_A_[x]
                                             B_p+_A; x C_[A]--A_[x]''').rxncon_system)
    print()
    for rule in rbm.rules:
        print(rule)