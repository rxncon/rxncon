import pytest
import rxncon.venntastic.sets as venn
import rxncon.simulation.bBM.bipartite_boolean_model as bbm
import rxncon.syntax.rxncon_from_string as rfs

def test_Factor():
    factor= venn.Intersection(venn.PropertySet(rfs.state_from_string("A--B")),
                              venn.PropertySet(rfs.state_from_string("A-{P}")))
    bbm_factor= bbm.Factor(factor)
    assert bbm_factor.is_equivalent_to(factor)
