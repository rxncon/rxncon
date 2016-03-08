import pytest
import rxncon.venntastic.sets as venn
import rxncon.simulation.bBM.bipartite_boolean_model as bbm
import rxncon.syntax.rxncon_from_string as rfs

def test_Factor():
    factor= venn.Intersection(venn.PropertySet(bbm.Node(rfs.state_from_string("A--B"))),
                              venn.PropertySet(bbm.Node(rfs.state_from_string("A-{P}"))))
    bbm_factor= bbm.Factor(factor)
    assert bbm_factor.factor.is_equivalent_to(factor)

def test_target():
    # todo: single state A
    target_node = bbm.Node(rfs.state_from_string("A--B"))
    assert target_node.value == rfs.state_from_string("A--B")


def test_rule():
    factor = venn.Intersection(venn.PropertySet(bbm.Node(rfs.state_from_string("A--B"))),
                               venn.PropertySet(bbm.Node(rfs.state_from_string("A-{P}"))))
    rule = bbm.Rule(bbm.Node(rfs.state_from_string("A--B")), bbm.Factor(factor))
    assert rule.target == bbm.Node(rfs.state_from_string("A--B"))
    assert rule.factor.factor.is_equivalent_to(factor)

