from rxncon.input.quick.quick import Quick
from rxncon.simulation.rule_based.rule_based_model import rule_based_model_from_rxncon, calc_positive_solutions
from rxncon.core.state import state_from_str

def test_simple_system():
    rbm = rule_based_model_from_rxncon(Quick("""A_[b]_ppi+_B_[a]; ! A_[(r)]-{p}
                                             A_[b]_ppi-_B_[a]
                                             C_p+_A_[(r)]
                                             D_p-_A_[(r)]""").rxncon_system)

    print()
    for rule in rbm.rules:
        print(rule)


def test_calc_positive_solutions():
    rxncon_sys = Quick("""A_[b]_ppi+_B_[a]; ! A_[(r)]-{p}
                          E_[x]_ppi+_B_[a]
                          C_p+_A_[(r)]
                          D_ub+_A_[(r)]""").rxncon_system

    soln = {state_from_str('A@1_[b]--B@2_[a]'): False}

    pos_solns = calc_positive_solutions(rxncon_sys, soln)

    for x in pos_solns:
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

