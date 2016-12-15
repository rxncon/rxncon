from rxncon.input.shared.contingency_list import contingency_list_entry_from_strs as cle_from_str, contingencies_from_contingency_list_entries
from rxncon.venntastic.sets import venn_from_str, Set as VennSet, ValueSet, Intersection, Union, Complement
from rxncon.core.state import state_from_str
from rxncon.core.effector import Effector, StateEffector, AndEffector, OrEffector, NotEffector

from rxncon.core.reaction import reaction_from_str

def venn_from_effector(effector: Effector) -> VennSet:
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


def test_nested_boolean():
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


def test_simple_namespace():
    cles = [
        cle_from_str('<C0>', 'AND', 'B@1_[(r)]-{p}'),
        cle_from_str('<C0>', 'AND', 'A@0_[ab]--B@1_[ba]'),
        cle_from_str('<C1>', 'AND', 'A@1_[ac]--C@2_[ca]'),
        cle_from_str('<C1>', 'AND', 'A@1_[ad]--D@3_[da]'),
        cle_from_str('<C>', 'AND', '<C0>'),
        cle_from_str('<C>', 'AND', '<C1>'),
        cle_from_str('<C>', 'AND', 'X@0_[xa]--A@1_[ax]'),
        cle_from_str('<C>', 'EQV', '<C0>.A@0,<C1>.A@1'),
        cle_from_str('<C>', 'EQV', 'A@1,<C0>.A@0'),
        cle_from_str('X_p+_A_[(r)]', '!', '<C>#1,A@1#0,X@0'),
    ]

    contingencies = contingencies_from_contingency_list_entries(cles)

    effector = contingencies[0].to_structured().effector
    # Note: the precise numbering is not important, but currently is reproducible.
    assert state_from_str('X@0_[xa]--A@1_[ax]') in effector.states
    assert state_from_str('B@2_[(r)]-{p}') in effector.states
    assert state_from_str('A@1_[ab]--B@2_[ba]') in effector.states
    assert state_from_str('A@1_[ac]--C@3_[ca]') in effector.states
    assert state_from_str('A@1_[ad]--D@4_[da]') in effector.states
