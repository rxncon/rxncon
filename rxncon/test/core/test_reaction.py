import rxncon.core.reaction as rxn
import rxncon.core.state as sta
import rxncon.input.shared.from_string as fst


def test_states_from_reaction_ppi():
    reaction = fst.reaction_from_string('Fus3_[CD]_ppi_Msg5_[n]')
    states = rxn.states_from_reaction(reaction)

    assert states.source_state is None
    assert isinstance(states.product_state, sta.InteractionState)
    assert states.product_state.full_name == 'Fus3_[CD]--Msg5_[n]'


def test_states_from_reaction_p_plus():
    reaction = fst.reaction_from_string('Fus3_[KD]_P+_Sst2_[(S539)]')
    states = rxn.states_from_reaction(reaction)

    assert states.source_state is None
    assert isinstance(states.product_state, sta.CovalentModificationState)
    assert states.product_state.full_name == 'Sst2_[(S539)]-{P}'


def test_states_from_reaction_p_minus():
    reaction = fst.reaction_from_string('Msg5_[PD]_P-_Slt2_[(Y192)]')
    states = rxn.states_from_reaction(reaction)

    assert states.product_state is None
    assert isinstance(states.source_state, sta.CovalentModificationState)
    assert states.source_state.full_name == 'Slt2_[(Y192)]-{P}'


def test_states_from_reaction_gef():
    reaction = fst.reaction_from_string('Rom2_[DH]_GEF_Rho1_[GnP]')
    states = rxn.states_from_reaction(reaction)

    assert states.source_state is None
    assert isinstance(states.product_state, sta.CovalentModificationState)
    assert states.product_state.full_name == 'Rho1_[GnP]-{P}'


def test_states_from_reaction_gap():
    reaction = fst.reaction_from_string('Lrg1_[GAP]_GAP_Rho1_[GnP]')
    states = rxn.states_from_reaction(reaction)

    assert states.product_state is None
    assert isinstance(states.source_state, sta.CovalentModificationState)
    assert states.source_state.full_name == 'Rho1_[GnP]-{P}'


def test_states_from_reaction_ub_plus():
    reaction = fst.reaction_from_string('SCF_Ub+_Tec1')
    states = rxn.states_from_reaction(reaction)

    assert states.source_state is None
    assert isinstance(states.product_state, sta.CovalentModificationState)
    assert states.product_state.full_name == 'Tec1-{Ub}'


def test_states_from_reaction_ap():
    reaction = fst.reaction_from_string('Rck2_AP_Rck2_[Ser]')
    states = rxn.states_from_reaction(reaction)

    assert states.source_state is None
    assert isinstance(states.product_state, sta.CovalentModificationState)
    assert states.product_state.full_name == 'Rck2_[Ser]-{P}'


def test_states_from_reaction_cut():
    reaction = fst.reaction_from_string('Yps1_CUT_Msb2_[HMH/CD]')
    states = rxn.states_from_reaction(reaction)

    assert states.source_state is None
    assert isinstance(states.product_state, sta.CovalentModificationState)
    assert states.product_state.full_name == 'Msb2_[HMH/CD]-{Truncated}'

