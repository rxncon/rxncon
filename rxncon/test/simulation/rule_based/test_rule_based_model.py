from itertools import combinations
import pytest

from rxncon.input.quick.quick import Quick
from rxncon.simulation.rule_based.rule_based_model import rule_based_model_from_rxncon, rule_from_str, complex_from_str, \
    initial_condition_from_str


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
        ['A(x) + A(x,r~p) -> A(x!1,r~p).A(x!1) kk', 'A(x) + A(x,r~p) -> A(x!2).A(x!2,r~p) kk'],
        ['C() + A(rR~0) -> C() + A(rR~p) k', 'A(rR~0) + C() -> A(rR~p) + C() k']
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


# Test simple rxncon systems.

def test_single_requirement_modification():
    rbm = rule_based_model_from_rxncon(Quick("""A_[b]_ppi+_B_[a]; ! A_[(r)]-{p}
                                             A_[b]_ppi-_B_[a]
                                             C_p+_A_[(r)]
                                             D_p-_A_[(r)]""").rxncon_system)

    expected_rules = [
        'A(rR~p,bD) + B(aD) -> A(rR~p,bD!1).B(aD!1) k',
        'A(bD!1).B(aD!1) -> A(bD) + B(aD) k',
        'C() + A(rR~0) -> C() + A(rR~p) k',
        'A(rR~p) + D() -> A(rR~0) + D() k'
    ]

    assert len(rbm.rules) == len(expected_rules)

    for actual_rule in rbm.rules:
        assert any(rule_from_str(rule).is_equivalent_to(actual_rule) for rule in expected_rules)

    expected_ics = [
        'A(bD,rR~0) NumA',
        'B(aD) NumB',
        'D() NumD',
        'C() NumC'
    ]

    for actual_ic in rbm.initial_conditions:
        assert any(initial_condition_from_str(ic).is_equivalent_to(actual_ic) for ic in expected_ics)


def test_single_inhibition_modification():
    rbm = rule_based_model_from_rxncon(Quick("""A_[b]_ppi+_B_[a]; x A_[(r)]-{p}
                                             E_[x]_ppi+_B_[a]
                                             C_p+_A_[(r)]
                                             D_ub+_A_[(r)]""").rxncon_system)

    expected_rules = [
        'A(rR~ub,bD) + B(aD) -> A(rR~ub,bD!1).B(aD!1) k',
        'A(rR~0,bD) + B(aD) -> A(rR~0,bD!1).B(aD!1) k',
        'B(aD) + E(xD) -> B(aD!1).E(xD!1) k',
        'A(rR~0) + C() -> A(rR~p) + C() k',
        'A(rR~0) + D() -> A(rR~ub) + D() k'
    ]

    assert len(rbm.rules) == len(expected_rules)

    for actual_rule in rbm.rules:
        assert any(rule_from_str(rule).is_equivalent_to(actual_rule) for rule in expected_rules)

    expected_ics = [
        'A(bD,rR~0) NumA',
        'B(aD) NumB',
        'D() NumD',
        'C() NumC',
        'E(xD) NumE'
    ]

    for actual_ic in rbm.initial_conditions:
        assert any(initial_condition_from_str(ic).is_equivalent_to(actual_ic) for ic in expected_ics)


def test_single_inhibition_interaction():
    rbm = rule_based_model_from_rxncon(Quick("""A_[b]_ppi+_B_[a]; x A_[(r)]-{p}
                                             E_[x]_ppi+_B_[a]
                                             C_p+_A_[(r)]
                                             D_ub+_A_[(r)]""").rxncon_system)

    expected_rules = [
        'A(rR~ub,bD) + B(aD) -> A(rR~ub,bD!1).B(aD!1) k',
        'A(rR~0,bD) + B(aD) -> A(rR~0,bD!1).B(aD!1) k',
        'B(aD) + E(xD) -> B(aD!1).E(xD!1) k',
        'A(rR~0) + C() -> A(rR~p) + C() k',
        'A(rR~0) + D() -> A(rR~ub) + D() k'
    ]

    assert len(rbm.rules) == len(expected_rules)

    for actual_rule in rbm.rules:
        assert any(rule_from_str(rule).is_equivalent_to(actual_rule) for rule in expected_rules)

    expected_ics = [
        'A(bD,rR~0) NumA',
        'B(aD) NumB',
        'D() NumD',
        'C() NumC',
        'E(xD) NumE'
    ]

    for actual_ic in rbm.initial_conditions:
        assert any(initial_condition_from_str(ic).is_equivalent_to(actual_ic) for ic in expected_ics)


def test_boolean_requirement_interaction():
    rbm = rule_based_model_from_rxncon(Quick('''A_ppi+_C
                                             C_ppi+_D
                                             B_ppi+_E
                                             B_ppi+_F
                                             A_ppi+_B; ! <comp1>
                                             <comp1>; OR <comp1C1>
                                             <comp1>; OR <comp2C1>
                                             <comp1C1>; AND A--C
                                             <comp1C1>; AND C--D
                                             <comp2C1>; AND B--F
                                             <comp2C1>; AND B--E''').rxncon_system)

    expected_rules = [
        'A(CD) + C(AD) -> A(CD!1).C(AD!1) k',
        'C(DD) + D(CD) -> C(DD!1).D(CD!1) k',
        'B(ED) + E(BD) -> B(ED!1).E(BD!1) k',
        'B(FD) + F(BD) -> B(FD!1).F(BD!1) k',
        'A(BD,CD) + B(AD,ED!2,FD!1).E(BD!2).F(BD!1) -> A(BD!3,CD).B(AD!3,ED!2,FD!1).E(BD!2).F(BD!1) k',  # B--F, B--E, NOT A--C
        'A(BD,CD!2).C(AD!2,DD) + B(AD,ED!3,FD!1).E(BD!3).F(BD!1) -> A(BD!4,CD!2).B(AD!4,ED!3,FD!1).C(AD!2,DD).E(BD!3).F(BD!1) k',  # B--F, B--E, A--C, NOT C--D
        'A(BD,CD!1).C(AD!1,DD!2).D(CD!2) + B(AD) -> A(BD!3,CD!1).B(AD!3).C(AD!1,DD!2).D(CD!2) k',  # A--C, C--D
    ]

    assert len(rbm.rules) == len(expected_rules)

    for actual_rule in rbm.rules:
        assert any(rule_from_str(rule).is_equivalent_to(actual_rule) for rule in expected_rules)


def test_boolean_inhibition_interaction():
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

    expected_rules = [
        'A(CD) + C(AD) -> A(CD!1).C(AD!1) k',
        'C(DD) + D(CD) -> C(DD!1).D(CD!1) k',
        'B(ED) + E(BD) -> B(ED!1).E(BD!1) k',
        'B(FD) + F(BD) -> B(FD!1).F(BD!1) k',
        'A(BD,CD) + B(AD,FD) -> A(BD!1,CD).B(AD!1,FD) k',
        'A(BD,CD) + B(AD,ED,FD!1).F(BD!1) -> A(BD!2,CD).B(AD!2,ED,FD!1).F(BD!1) k',
        'A(BD,CD!1).C(AD!1,DD) + B(AD,FD) -> A(BD!2,CD!1).B(AD!2,FD).C(AD!1,DD) k',
        'A(BD,CD!1).C(AD!1,DD) + B(AD,ED,FD!2).F(BD!2) -> A(BD!3,CD!1).B(AD!3,ED,FD!2).C(AD!1,DD).F(BD!2) k',
    ]

    assert len(rbm.rules) == len(expected_rules)

    for actual_rule in rbm.rules:
        assert any(rule_from_str(rule).is_equivalent_to(actual_rule) for rule in expected_rules)


def test_mutually_exclusive_bindings():
    rbm = rule_based_model_from_rxncon(Quick('''C_[A]_ppi+_A_[x]
                                             D_[A]_ppi+_A_[x]
                                             B_p+_A; x C_[A]--A_[x]''').rxncon_system)

    expected_rules = [
        'A(xD) + C(AD) -> A(xD!1).C(AD!1) k',
        'A(xD) + D(AD) -> A(xD!1).D(AD!1) k',
        'A(BR~0,xD) + B() -> A(BR~p,xD) + B() k',
        'A(BR~0,xD!1).D(AD!1) + B() -> A(BR~p,xD!1).D(AD!1) + B() k'
    ]

    assert len(rbm.rules) == len(expected_rules)

    for actual_rule in rbm.rules:
        assert any(rule_from_str(rule).is_equivalent_to(actual_rule) for rule in expected_rules)
