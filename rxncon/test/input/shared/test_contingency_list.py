import pytest
from rxncon.input.shared.contingency_list import contingency_list_entry_from_strs as cle_from_str, contingencies_from_contingency_list_entries
from rxncon.venntastic.sets import venn_from_str, Set as VennSet, ValueSet, Intersection, Union, Complement
from rxncon.core.state import state_from_str, State
from rxncon.core.effector import Effector, StateEffector, AndEffector, OrEffector, NotEffector


def venn_from_effector(effector: Effector) -> VennSet[State]:
    if isinstance(effector, StateEffector):
        return ValueSet(effector.expr)
    elif isinstance(effector, AndEffector):
        return Intersection(*(venn_from_effector(x) for x in effector.exprs))
    elif isinstance(effector, OrEffector):
        return Union(*(venn_from_effector(x) for x in effector.exprs))
    elif isinstance(effector, NotEffector):
        return Complement(venn_from_effector(effector.expr))
    else:
        raise AssertionError


def test_nested_boolean() -> None:
    cles = [
        cle_from_str('<C1>', 'AND', 'A_[x]--B_[y]'),
        cle_from_str('<C1>', 'AND', 'A_[(r)]-{p}'),
        cle_from_str('<C2>', 'AND', 'B_[z]--D_[y]'),
        cle_from_str('<C2>', 'AND', 'B_[(r1)]-{p}'),
        cle_from_str('<C2>', 'AND', 'B_[(r2)]-{p}'),
        cle_from_str('<C1C2>', 'AND', '<C1>'),
        cle_from_str('<C1C2>', 'AND', '<C2>'),
        cle_from_str('A_[q]_ppi+_Q_[a]', '!', '<C1C2>')
    ]

    contingencies = contingencies_from_contingency_list_entries(cles)
    assert len(contingencies) == 1

    expected = venn_from_str('( A_[x]--B_[y] & A_[(r)]-{p} & B_[z]--D_[y] & B_[(r1)]-{p} & B_[(r2)]-{p} )', state_from_str)
    assert venn_from_effector(contingencies[0].effector).is_equivalent_to(expected)


def test_contingency_equality() -> None:
    cles1 = [
        cle_from_str('<C1>', 'AND', 'A_[x]--B_[y]'),
        cle_from_str('<C1>', 'AND', 'A_[(r)]-{p}'),
        cle_from_str('<C2>', 'NOT', 'B_[(r1)]-{p}'),
        cle_from_str('<C1C2>', 'OR', '<C1>'),
        cle_from_str('<C1C2>', 'OR', '<C2>'),
        cle_from_str('A_[q]_ppi+_Q_[a]', '!', '<C1C2>')
    ]

    cles2 = [
        cle_from_str('<C1>', 'AND', 'A_[x]--B_[y]'),
        cle_from_str('<C1>', 'AND', 'A_[(r)]-{p}'),
        cle_from_str('<C2>', 'NOT', 'B_[(r1)]-{p}'),
        cle_from_str('<C1C2>', 'OR', '<C1>'),
        cle_from_str('<C1C2>', 'OR', '<C2>'),
        cle_from_str('A_[q]_ppi+_Q_[a]', '!', '<C1C2>')
    ]
    contingencies1 = contingencies_from_contingency_list_entries(cles1)
    assert len(contingencies1) == 1
    contingencies2 = contingencies_from_contingency_list_entries(cles2)
    assert len(contingencies2) == 1
    assert contingencies1 == contingencies2


def test_nested_OR_NOT_boolean() -> None:
    cles = [
        cle_from_str('<C1>', 'AND', 'A_[x]--B_[y]'),
        cle_from_str('<C1>', 'AND', 'A_[(r)]-{p}'),
        cle_from_str('<C2>', 'NOT', 'B_[(r1)]-{p}'),
        cle_from_str('<C1C2>', 'OR', '<C1>'),
        cle_from_str('<C1C2>', 'OR', '<C2>'),
        cle_from_str('A_[q]_ppi+_Q_[a]', '!', '<C1C2>')
    ]

    contingencies = contingencies_from_contingency_list_entries(cles)
    assert len(contingencies) == 1

    expected = venn_from_str('( ( A_[x]--B_[y] & A_[(r)]-{p} ) | ~( B_[(r1)]-{p} ) )', state_from_str)
    assert venn_from_effector(contingencies[0].effector).is_equivalent_to(expected)


def test_simple_namespace() -> None:
    cles = [
        cle_from_str('<C0>', 'AND', 'B@1_[(r)]-{p}'),
        cle_from_str('<C0>', 'AND', 'A@3_[ab]--B@1_[ba]'),
        cle_from_str('<C1>', 'AND', 'A@2_[ac]--C@4_[ca]'),
        cle_from_str('<C1>', 'AND', 'A@2_[ad]--D@5_[da]'),
        cle_from_str('<C>', 'AND', '<C0>#A@1=A@3'),
        cle_from_str('<C>', 'AND', '<C1>#<C0>.A@3=A@2'),
        cle_from_str('<C>', 'AND', 'X@0_[xa]--A@1_[ax]'),
        cle_from_str('X_p+_A_[(r)]', '!', '<C>#A@1=A@1#X@0=X@0'),
    ]

    contingencies = contingencies_from_contingency_list_entries(cles)

    effector = contingencies[0].to_structured().effector

    # Note: the precise numbering is not important, but currently is reproducible.
    assert state_from_str('X@0_[xa]--A@1_[ax]') in effector.states
    assert state_from_str('B@2_[(r)]-{p}') in effector.states
    assert state_from_str('A@1_[ab]--B@2_[ba]') in effector.states
    assert state_from_str('A@1_[ac]--C@3_[ca]') in effector.states
    assert state_from_str('A@1_[ad]--D@4_[da]') in effector.states


def test_insulin_homodimer() -> None:
    cles = [
        cle_from_str('IR_ap+_IR_[JM(Y972)]', '!', '<IRstar>#IR@0=<IR-Ins>.IR@0#IR@2=<IR-Ins>.<IR01-Ins>.IR@1#IR@2=<IR-Ins>.IR@2#insulin@3=<IR-Ins>.<IR01-Ins>.insulin@2#insulin@4=<IR-Ins>.<IR01-Ins>.insulin@3'),
        cle_from_str('<IRstar>', 'AND', '<IR-phos>'),
        cle_from_str('<IRstar>', 'AND', '<IR-Ins>#<IR-phos>.IR@0=IR@0'),

        cle_from_str('<IR-phos>', 'AND', 'IR@0_[TK(Y1158)]-{P}'),
        cle_from_str('<IR-phos>', 'AND', 'IR@0_[TK(Y1162)]-{P}'),
        cle_from_str('<IR-phos>', 'AND', 'IR@0_[TK(Y1163)]-{P}'),

        cle_from_str('<IR-Ins>', 'AND', 'IR@0_[IRBD]--IR@2_[IRBD]'),
        cle_from_str('<IR-Ins>', 'AND', '<IR01-Ins>#IR@0=IR@0'),##IR@1=IR@1

        cle_from_str('<IR01-Ins>', 'OR', 'IR@0_[lig]--insulin@2_[IR]'),
        cle_from_str('<IR01-Ins>', 'OR', 'IR@1_[lig]--insulin@3_[IR]'),
    ]

    states = contingencies_from_contingency_list_entries(cles)[0].to_structured().effector.states

    assert state_from_str('IR@0_[TK(Y1158)]-{p}') in states
    assert state_from_str('IR@0_[TK(Y1162)]-{p}') in states
    assert state_from_str('IR@0_[TK(Y1163)]-{p}') in states
    assert state_from_str('IR@0_[IRBD]--IR@2_[IRBD]') in states
    assert state_from_str('IR@0_[lig]--insulin@3_[IR]') in states
    assert state_from_str('IR@2_[lig]--insulin@4_[IR]') in states


def test_contingency_list_entry_from_strs():
    cles = [
        cle_from_str('A_p+_B_[(a)]', '!', '<bool>#A@0=<bool1>.A@1#B@1=<bool2>.<bool3>.B@3#E@2=<bool2>.<bool3>.E@3'),
        cle_from_str('<bool>', 'OR', '<bool1>#<bool2>.<bool3>.A@2=A@1'),
        cle_from_str('<bool>', 'OR', '<bool2>#<bool1>.B@2=<bool3>.B@3'),
        cle_from_str('<bool1>', 'AND', 'A@1_[ab1]--B@2_[ba1]'),
        cle_from_str('<bool1>', 'AND', 'B@2_[be1]--E@3_[eb1]'),
        cle_from_str('<bool2>', 'AND', '<bool3>#<bool4>.A@2=A@2'),
        cle_from_str('<bool2>', 'AND', '<bool4>'),
        cle_from_str('<bool3>', 'OR', 'A@2_[ab2]--B@3_[ba2]'),
        cle_from_str('<bool3>', 'OR', 'A@2_[ae2]--E@3_[ea2]'),
        cle_from_str('<bool4>', 'OR', 'A@2_[(a1)]-{p}'),
        cle_from_str('<bool4>', 'OR', 'A@2_[(a2)]-{ub}'),
    ]

    effector = contingencies_from_contingency_list_entries(cles)[0].to_structured().effector
    states = effector.states

    assert state_from_str('A@0_[ab1]--B@1_[ba1]') in states
    assert state_from_str('B@1_[be1]--E@3_[eb1]') in states
    assert state_from_str('A@0_[ab2]--B@1_[ba2]') in states
    assert state_from_str('A@0_[ae2]--E@2_[ea2]') in states
    assert state_from_str('A@0_[(a1)]-{p}') in states
    assert state_from_str('A@0_[(a2)]-{ub}') in states


def test_single_structured_state_wrong_indexing() -> None:
    cles = [
        cle_from_str('A_p+_B_[(r)]', '!', 'B@0--A@1'),
    ]
    with pytest.raises(AssertionError):
        contingencies_from_contingency_list_entries(cles)[0].to_structured().effector.states


def test_single_structured_state_spec_not_in_reaction() -> None:
    cles = [
        cle_from_str('A_p+_B_[(r)]', '!', 'C@0-{p}'),
    ]
    with pytest.raises(AssertionError):
        contingencies_from_contingency_list_entries(cles)[0].to_structured().effector.states


def test_unary_and_contingency() -> None:
    cles = [
        cle_from_str('A_[b]_ppi+_B_[a]', '!', '<x>'),
        cle_from_str('<x>', 'AND', 'A-{p}')
    ]

    assert isinstance(contingencies_from_contingency_list_entries(cles)[0].effector, StateEffector)


def test_unary_or_contingency() -> None:
    cles = [
        cle_from_str('A_[b]_ppi+_B_[a]', '!', '<x>'),
        cle_from_str('<x>', 'OR', 'A-{p}')
    ]

    assert isinstance(contingencies_from_contingency_list_entries(cles)[0].effector, StateEffector)
