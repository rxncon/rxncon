import pytest
from itertools import combinations

from rxncon.input.quick.quick import Quick
from rxncon.simulation.rule_based.rule_based_model import complex_from_str, rule_from_str, rule_based_model_from_rxncon,\
    initial_condition_from_str, bond_complexes, components_microstate, BondComplex
from rxncon.core.state import state_from_str
from rxncon.venntastic.sets import ValueSet, Union, venn_from_str
import rxncon.core.reaction as reaction
from rxncon.core.spec import Spec, LocusResolution, MRNASpec, spec_from_str


# Test the *_from_str functions.
def test_complex_from_str_equivalent() -> None:
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


def test_complex_from_str_inequivalent() -> None:
    complex_classes = [
        ['A()', 'B()'],
        ['A(x!1,r~p).B(rr~0,y!2,cc!1).C(c!2)', 'B(rr~0,y!4,cc!1).A(x!4,r~p).C(c!1)']
    ]

    for complex_class in complex_classes:
        for first_complex, second_complex in combinations(complex_class, 2):
            assert not complex_from_str(first_complex).is_equivalent_to(complex_from_str(second_complex))
            assert not complex_from_str(second_complex).is_equivalent_to(complex_from_str(first_complex))


def test_invalid_complexes() -> None:
    # Two molecules, not connected by any bond.
    with pytest.raises(AssertionError):
        complex_from_str('A().B()')

    # Single dangling bond.
    with pytest.raises(AssertionError):
        complex_from_str('A(x!1)')

    # Two molecules with two dangling bonds.
    with pytest.raises(AssertionError):
        complex_from_str('A(x!1).B(y!2)')


def test_rule_from_str_equivalent() -> None:
    rule_classes = [
        ['A(x) + B(y) -> A(x!1).B(y!1) k', 'B(y) + A(x) -> A(x!1).B(y!1) k', 'B(y) + A(x) -> B(y!4).A(x!4) k'],
        ['A(x) + A(x,r~p) -> A(x!1,r~p).A(x!1) kk', 'A(x) + A(x,r~p) -> A(x!2).A(x!2,r~p) kk'],
        ['C() + A(rR~0) -> C() + A(rR~p) k', 'A(rR~0) + C() -> A(rR~p) + C() k']
    ]

    for rule_class in rule_classes:
        for first_rule, second_rule in combinations(rule_class, 2):
            assert rule_from_str(first_rule).is_equivalent_to(rule_from_str(second_rule))
            assert rule_from_str(second_rule).is_equivalent_to(rule_from_str(first_rule))


def test_rule_from_str_inequivalent() -> None:
    rule_classes = [
        ['A(x) + A(y,r~p) -> A(x!1,r~p).A(y!1) kk', 'A(x) + A(y,r~p) -> A(y!1).A(x!1,r~p) kk']
    ]

    for rule_class in rule_classes:
        for first_rule, second_rule in combinations(rule_class, 2):
            assert rule_from_str(first_rule).is_equivalent_to(rule_from_str(second_rule))
            assert rule_from_str(second_rule).is_equivalent_to(rule_from_str(first_rule))


def test_bond_complex_class() -> None:
    first_complex = BondComplex({spec_from_str('A@0')}, {state_from_str('A@0_[C]--C@2_[A]'): True}, set(), set())
    second_complex = BondComplex({spec_from_str('C@2')}, {state_from_str('A@0_[C]--C@2_[A]'): True}, set(), set())

    assert first_complex.dangling_bonds() == {state_from_str('A@0_[C]--C@2_[A]')}
    assert second_complex.dangling_bonds() == {state_from_str('A@0_[C]--C@2_[A]')}
    assert first_complex.can_connect_with(second_complex)

    first_and_second = first_complex.combined_with(second_complex)
    assert first_complex.dangling_bonds() == {state_from_str('A@0_[C]--C@2_[A]')}
    assert second_complex.dangling_bonds() == {state_from_str('A@0_[C]--C@2_[A]')}
    assert first_and_second.is_connected()
    assert not first_and_second.dangling_bonds()


def test_bond_complexes_calculation() -> None:
    states = Union(*(ValueSet(state_from_str(x))
                     for x in ('A@0_[C]--C@2_[A]', 'C@2_[(r)]-{p}',
                               'A@0_[C]--0', 'C@2_[A]--0', 'C@2_[(r)]-{0}')))

    assert len(bond_complexes(states)) == 3


def test_bond_complexes_ring_calculation() -> None:
    states = Union(*(ValueSet(state_from_str(x))
                     for x in ('A@0_[B]--B@2_[A]', 'B@2_[C]--C@3_[B]',
                               'C@3_[D]--D@4_[C]', 'D@4_[A]--A@0_[D]')))

    assert len(bond_complexes(states)) == 11


def test_bond_complexes_poly_calculation() -> None:
    states = Union(*(ValueSet(state_from_str(x))
                     for x in ('R@1_[Y]--Y@2_[R]', 'Y@2_[Y]--Y@3_[Y]',
                               'Y@3_[R]--R@4_[Y]', 'R@4_[Ras]--Ras@5_[R]')))

    assert len(bond_complexes(states)) == 5


def test_bond_complexes_calculation_unordered() -> None:
    states = Union(*(ValueSet(state_from_str(x))
                     for x in ('A@0_[C]--C@5_[A]', 'C@5_[B]--B@3_[C]', 'B@3_[D]--D@2_[B]')))

    assert len(bond_complexes(states)) == 4


# Test simple rxncon systems.
def test_single_requirement_modification() -> None:
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
        print(actual_rule)
        # assert any(rule_from_str(rule).is_equivalent_to(actual_rule) for rule in expected_rules)

    expected_ics = [
        'A(bD,rR~0) NumA',
        'B(aD) NumB',
        'D() NumD',
        'C() NumC'
    ]

    for actual_ic in rbm.initial_conditions:
        assert any(initial_condition_from_str(ic).is_equivalent_to(actual_ic) for ic in expected_ics)


def test_single_inhibition_modification() -> None:
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


def test_single_inhibition_interaction() -> None:
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


def test_boolean_requirement_interaction() -> None:
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
        'A(BD,CD) + B(AD,ED!2,FD!1).E(BD!2).F(BD!1) -> A(BD!3,CD).B(AD!3,ED!2,FD!1).E(BD!2).F(BD!1) k',
        'A(BD,CD!2).C(AD!2,DD) + B(AD,ED!3,FD!1).E(BD!3).F(BD!1) -> '
        'A(BD!4,CD!2).B(AD!4,ED!3,FD!1).C(AD!2,DD).E(BD!3).F(BD!1) k',
        'A(BD,CD!1).C(AD!1,DD!2).D(CD!2) + B(AD,FD) -> A(BD!3,CD!1).B(AD!3,FD).C(AD!1,DD!2).D(CD!2) k',
        'A(BD,CD!1).C(AD!1,DD!2).D(CD!2) + B(AD,FD!4).F(BD!4) -> '
        'A(BD!3,CD!1).B(AD!3,FD!4).F(BD!4).C(AD!1,DD!2).D(CD!2) k'
    ]

    assert len(rbm.rules) == len(expected_rules)

    for actual_rule in rbm.rules:
        assert any(rule_from_str(rule).is_equivalent_to(actual_rule) for rule in expected_rules)


def test_boolean_inhibition_interaction() -> None:
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


def test_mutually_exclusive_bindings() -> None:
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


def test_kplus_kminus() -> None:
    rbm = rule_based_model_from_rxncon(Quick("""A_[b]_ppi+_B_[a]; k+ A_[(r1)]-{p}
                                                A_[b]_ppi+_B_[a]; k- A_[(r2)]-{p}
                                                C_p+_A_[(r1)]
                                                D_p+_A_[(r2)]""").rxncon_system)

    expected_rules = [
        'A(r1R~p,r2R~p,bD) + B(aD) -> A(r1R~p,r2R~p,bD!1).B(aD!1) k',
        'A(r1R~p,r2R~0,bD) + B(aD) -> A(r1R~p,r2R~0,bD!1).B(aD!1) k',
        'A(r1R~0,r2R~p,bD) + B(aD) -> A(r1R~0,r2R~p,bD!1).B(aD!1) k',
        'A(r1R~0,r2R~0,bD) + B(aD) -> A(r1R~0,r2R~0,bD!1).B(aD!1) k',
        'A(r1R~0) + C() -> A(r1R~p) + C() k',
        'A(r2R~0) + D() -> A(r2R~p) + D() k'
    ]

    assert len(rbm.rules) == len(expected_rules)

    for actual_rule in rbm.rules:
        assert any(rule_from_str(rule).is_equivalent_to(actual_rule) for rule in expected_rules)


def test_self_regulation() -> None:
    rxn_system = Quick("""Rlm1_[MADS]_bind+_Rlm1Gene
                          PolII_trsc_Rlm1Gene
                          Ribo_trsl_Rlm1mRNA""").rxncon_system

    rbm = rule_based_model_from_rxncon(rxn_system)

    expected_rules = [
        'Rlm1(MADSD) + Rlm1Gene(Rlm1D) -> Rlm1(MADSD!1).Rlm1Gene(Rlm1D!1) k',
        'PolII() + Rlm1Gene() -> PolII() + Rlm1Gene() + Rlm1mRNA() k',
        'Ribo() + Rlm1mRNA() -> Ribo() + Rlm1(MADSD) + Rlm1mRNA() k'
    ]

    assert len(rbm.rules) == len(expected_rules)

    for actual_rule in rbm.rules:
        assert any(rule_from_str(rule).is_equivalent_to(actual_rule) for rule in expected_rules)


def test_trslprocat() -> None:
    procatdef = reaction.ReactionDef(
        'pro-cat-translation',
        '$x_trslprocat_$y',
        {
            '$x': (Spec, LocusResolution.component),
            '$y': (MRNASpec, LocusResolution.component)
        },
        '$x%# + $y%# -> $x%# + $y%# + $y.to_protein_component_spec().with_name_suffix(\'PRO\')%!$y.to_protein_component_spec().with_name_suffix(\'CAT\')%#'  # pylint: disable=line-too-long
        '$y.to_protein_component_spec().with_name_suffix(\'PRO\').with_domain(\'PROCAT\')%--$y.to_protein_component_spec().with_name_suffix(\'CAT\').with_domain(\'CATPRO\')%!0'
    )

    reaction.REACTION_DEFS = reaction.DEFAULT_REACTION_DEFS + [procatdef]

    rxn_system = Quick("""Ribo_trslprocat_Ssy5mRNA
                          A_p+_Ssy5CAT
                          B_deg_Ssy5PRO""").rxncon_system

    rbm = rule_based_model_from_rxncon(rxn_system)

    expected_rules = [
        'Ribo() + Ssy5mRNA() -> Ribo() + Ssy5CAT(AR~0,CATPROD!1).Ssy5PRO(PROCATD!1) + Ssy5mRNA() k',
        'A() + Ssy5CAT(AR~0) -> A() + Ssy5CAT(AR~p) k',
        'B() + Ssy5PRO() -> B() k'
    ]

    assert len(rbm.rules) == len(expected_rules)

    for actual_rule in rbm.rules:
        assert any(rule_from_str(rule).is_equivalent_to(actual_rule) for rule in expected_rules)

    reaction.REACTION_DEFS = reaction.DEFAULT_REACTION_DEFS


def test_kplus_overlaps_with_reaction() -> None:
    rxncon_system = Quick("""A_p+_B_[(s)]; ! <Aphos> ; k+ A_[(r1)]-{p} ; k+ A_[(r2)]-{p}
                          <Aphos>; OR A_[(r1)]-{p}
                          <Aphos>; OR A_[(r2)]-{p}
                          C_p+_A_[(r1)]
                          D_p+_A_[(r2)]""").rxncon_system

    rbm = rule_based_model_from_rxncon(rxncon_system)

    for rule in rbm.rules:
        print(rule)