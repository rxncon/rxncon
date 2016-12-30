from itertools import combinations
import pytest

from rxncon.input.quick.quick import Quick
from rxncon.simulation.rule_based.rule_based_model import rule_based_model_from_rxncon, rule_from_str, complex_from_str


# Test the *_from_str functions.
def test_complex_from_str_equivalent():
    complex_classes = [
        ['A()', 'A()'],
        ['A(x!1).A(x!1)', 'A(x!2).A(x!2)'],
        ['A(x!1).A(y!1)', 'A(y!3).A(x!3)'],
        ['A(x!1,r~p).B(rr~0,y!1,cc!2).C(c!2)', 'B(rr~0,y!4,cc!1).A(x!4,r~p).C(c!1)']
    ]

    for complex_class in complex_classes:
        for first_complex, second_complex in combinations(complex_class, 2):
            assert complex_from_str(first_complex).is_equivalent_to(complex_from_str(second_complex))
            assert complex_from_str(second_complex).is_equivalent_to(complex_from_str(first_complex))


def test_complex_from_str_inequivalent():
    complex_classes = [
        ['A()', 'B()'],
        ['A(x!1,r~p).B(rr~0,y!2,cc!1).C(c!2)', 'B(rr~0,y!4,cc!1).A(x!4,r~p).C(c!1)']
    ]

    for complex_class in complex_classes:
        for first_complex, second_complex in combinations(complex_class, 2):
            assert not complex_from_str(first_complex).is_equivalent_to(complex_from_str(second_complex))
            assert not complex_from_str(second_complex).is_equivalent_to(complex_from_str(first_complex))


def test_invalid_complexes():
    # Two molecules, not connected by any bond.
    with pytest.raises(AssertionError):
        c = complex_from_str('A().B()')

    # Single dangling bond.
    with pytest.raises(AssertionError):
        c = complex_from_str('A(x!1)')

    # Two molecules with two dangling bonds.
    with pytest.raises(AssertionError):
        c = complex_from_str('A(x!1).B(y!2)')


def test_rule_from_str_equivalent():
    rule_classes = [
        ['A(x) + B(y) -> A(x!1).B(y!1) k', 'B(y) + A(x) -> A(x!1).B(y!1) k', 'B(y) + A(x) -> B(y!4).A(x!4) k'],
        ['A(x) + A(x,r~p) -> A(x!1,r~p).A(x!1) kk', 'A(x) + A(x,r~p) -> A(x!2).A(x!2,r~p) kk']
    ]

    for rule_class in rule_classes:
        for first_rule, second_rule in combinations(rule_class, 2):
            assert rule_from_str(first_rule).is_equivalent_to(rule_from_str(second_rule))
            assert rule_from_str(second_rule).is_equivalent_to(rule_from_str(first_rule))


def test_rule_from_str_inequivalent():
    rule_classes = [
        ['A(x) + A(y,r~p) -> A(x!1,r~p).A(y!1) kk', 'A(x) + A(y,r~p) -> A(y!1).A(x!1,r~p) kk']
    ]

    for rule_class in rule_classes:
        for first_rule, second_rule in combinations(rule_class, 2):
            assert rule_from_str(first_rule).is_equivalent_to(rule_from_str(second_rule))
            assert rule_from_str(second_rule).is_equivalent_to(rule_from_str(first_rule))



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

    rbm = rule_based_model_from_rxncon(rxncon_sys)

    for x in rbm.rules:
        print(x)

    print()
    for ic in rbm.initial_conditions:
        print(ic)


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

    print()
    for ic in rbm.initial_conditions:
        print(ic)

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