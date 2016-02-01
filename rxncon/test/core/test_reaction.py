import rxncon.core.reaction as rxn
import rxncon.core.state as sta
import rxncon.syntax.rxncon_from_string as fst


def test_states_from_reaction_ppi():
    reaction = fst.reaction_from_string('Fus3_[CD]_ppi_Msg5_[n]')
    states = rxn.states_from_reaction(reaction)

    assert states.source_state is None
    assert isinstance(states.product_state, sta.InterProteinInteractionState)
    assert str(states.product_state) == 'Fus3_[CD]--Msg5_[n]'


def test_states_from_reaction_p_plus():
    reaction = fst.reaction_from_string('Fus3_[KD]_P+_Sst2_[(S539)]')
    states = rxn.states_from_reaction(reaction)

    assert states.source_state is None
    assert isinstance(states.product_state, sta.CovalentModificationState)
    assert str(states.product_state) == 'Sst2_[(S539)]-{p}'


def test_states_from_reaction_p_minus():
    reaction = fst.reaction_from_string('Msg5_[PD]_P-_Slt2_[(Y192)]')
    states = rxn.states_from_reaction(reaction)

    assert states.product_state is None
    assert isinstance(states.source_state, sta.CovalentModificationState)
    assert str(states.source_state) == 'Slt2_[(Y192)]-{p}'


def test_states_from_reaction_gef():
    reaction = fst.reaction_from_string('Rom2_[DH]_GEF_Rho1_[GnP]')
    states = rxn.states_from_reaction(reaction)

    assert states.source_state is None
    assert isinstance(states.product_state, sta.CovalentModificationState)
    assert str(states.product_state) == 'Rho1_[GnP]-{p}'


def test_states_from_reaction_gap():
    reaction = fst.reaction_from_string('Lrg1_[GAP]_GAP_Rho1_[GnP]')
    states = rxn.states_from_reaction(reaction)

    assert states.product_state is None
    assert isinstance(states.source_state, sta.CovalentModificationState)
    assert str(states.source_state) == 'Rho1_[GnP]-{p}'


def test_states_from_reaction_ub_plus():
    reaction = fst.reaction_from_string('SCF_Ub+_Tec1')
    states = rxn.states_from_reaction(reaction)

    assert states.source_state is None
    assert isinstance(states.product_state, sta.CovalentModificationState)
    assert str(states.product_state) == 'Tec1-{ub}'


def test_states_from_reaction_ap():
    reaction = fst.reaction_from_string('Rck2_AP_Rck2_[Ser]')
    states = rxn.states_from_reaction(reaction)

    assert states.source_state is None
    assert isinstance(states.product_state, sta.CovalentModificationState)
    assert str(states.product_state) == 'Rck2_[Ser]-{p}'


def test_states_from_reaction_cut():
    reaction = fst.reaction_from_string('Yps1_CUT_Msb2_[HMH/CD]')
    states = rxn.states_from_reaction(reaction)

    assert states.source_state is None
    assert isinstance(states.product_state, sta.CovalentModificationState)
    assert str(states.product_state) == 'Msb2_[HMH/CD]-{truncated}'


def test_states_from_reaction_pt():
    reaction = fst.reaction_from_string('A_pt_B')
    states = rxn.states_from_reaction(reaction)

    assert states.source_state == fst.state_from_string('A-{p}')
    assert states.product_state == fst.state_from_string('B-{p}')


def test_states_from_reaction_synth():
    reaction = fst.reaction_from_string('A_syn_B')
    states = rxn.states_from_reaction(reaction)

    assert not states.source_state
    assert str(states.product_state) == 'B'
