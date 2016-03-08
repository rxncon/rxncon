import pytest
import rxncon.venntastic.sets as venn
import rxncon.simulation.bBM.bipartite_boolean_model as bbm
import rxncon.syntax.rxncon_from_string as rfs

def test_factor():
    factor= venn.Intersection(venn.PropertySet(rfs.state_from_string("A--B")),
                              venn.PropertySet(rfs.state_from_string("A-{P}")))
    bbm_factor= bbm.Factor(factor)
    assert bbm_factor.factor.is_equivalent_to(factor)

def test_target():
    target = bbm.Target("A")
    assert target.name=="A"

def test_rule():
    factor = venn.Intersection(venn.PropertySet(rfs.state_from_string("A--B")),
                               venn.PropertySet(rfs.state_from_string("A-{P}")))
    rule = bbm.Rule(bbm.Target("A"), bbm.Factor(factor))


    assert rule.target == bbm.Target("A")
    assert rule.factor.factor.is_equivalent_to(factor)
    #assert rule.factor==factor

