import pytest

import rxncon.simulation.bBM.bipartite_boolean_model as bbm
import rxncon.simulation.bBM.bBM_boolnet_exporter as bbe

import rxncon.syntax.rxncon_from_string as rfs
import rxncon.venntastic.sets as venn


def test_generate_name():
    a_pplus_b = venn.PropertySet(bbm.Node(rfs.reaction_from_string("a_p+_b")))
    assert bbe.string_from_reaction(a_pplus_b.value) == "a_pplus_b"

    a_dash_dash_b = venn.PropertySet(rfs.state_from_string("A--B"))
    assert bbe.string_from_inter_protein_interaction_state(a_dash_dash_b.value) == "A__B"

    b_intra = venn.PropertySet(rfs.state_from_string("b_[n]--[m]"))

    assert bbe.string_from_intra_protein_interaction_state(b_intra.value) == "b_OpenBreaket_n_CloseBreaket__OpenBreaket_m_CloseBreaket"

