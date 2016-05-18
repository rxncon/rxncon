from collections import namedtuple
import pytest
import rxncon.core.state as sta
import rxncon.syntax.rxncon_from_string as fst


# Component from string #
ComponentFromStringTestCase = namedtuple('ComponentFromStringTestCase', ['component_string', 'expected_name',
                                                                         'expected_domain', 'expected_subdomain',
                                                                         'expected_residue'])
def test_component_from_string(the_case_component_from_string):
    for the_case in the_case_component_from_string:
        is_component_correct(the_case)

def is_component_correct(the_case):
    component = fst.specification_from_string(the_case.component_string)
    assert str(component) == the_case.component_string
    assert component.name == the_case.expected_name
    assert component.domain == the_case.expected_domain
    assert component.subdomain == the_case.expected_subdomain
    assert component.residue == the_case.expected_residue


@pytest.fixture
def the_case_component_from_string():
    return [
        ComponentFromStringTestCase('Sln1', 'Sln1', None, None, None),
        ComponentFromStringTestCase('Pkc1_[HR1]', 'Pkc1', 'HR1', None, None),
        ComponentFromStringTestCase('Pbs2_[RSD2/PR]', 'Pbs2', 'RSD2', 'PR', None),
        ComponentFromStringTestCase('Sln1_[HK(H576)]', 'Sln1', 'HK', None, 'H576'),
        ComponentFromStringTestCase('A_[B/C(D)]', 'A', 'B', 'C', 'D'),
        ComponentFromStringTestCase('A_[(b)]', 'A', None, None, 'b')
    ]


# State from string
StateFromStringTestCase = namedtuple('StateFromStringTestCase', ['state_string', "expected_state", 'expected_first_component',
                                                                 'expected_second_component', 'expected_substrate',
                                                                 'expected_modifier'])


def _test_rxncon_from_string_state(the_case_rxncon_from_string_state):
    for the_case in the_case_rxncon_from_string_state:
        is_state_correct(the_case)


def is_state_correct(the_case):
    state = fst.state_from_string(the_case.state_string)
    assert isinstance(state, the_case.expected_state)
    assert str(state) == the_case.state_string
    if the_case.expected_first_component and the_case.expected_second_component:
        assert str(state.first_component) == the_case.expected_first_component
        assert str(state.second_component) == the_case.expected_second_component

    elif the_case.expected_substrate:
        assert str(state.substrate) == the_case.expected_substrate
        assert state.modifier == the_case.expected_modifier

    else:
        return NotImplementedError


@pytest.fixture
def the_case_rxncon_from_string_state():
    StateFromStringTestCase('Fus3_[CD]--Msg5_[n]', sta.InteractionState, 'Fus_[CD]', 'Msg5_[n]', None),
    StateFromStringTestCase('A_[n]--[m]', sta.SelfInteractionState, 'A_[n]', 'A_[m]', None),
    StateFromStringTestCase('Slt2_[(Y192)]-{p}', sta.CovalentModificationState, None, None, 'Slt2_[(Y192)]',
                            sta.StateModifier.phosphor)


def test_rxncon_from_string_state_intra_protein_interaction():
    state = fst.state_from_string('A_[n]--[m]')

    assert isinstance(state, sta.SelfInteractionState)
    assert str(state) == 'A_[n]--[m]'
    assert str(state.first_component) == 'A_[n]'
    assert str(state.second_component) == 'A_[m]'



def test_rxncon_from_string_state_phosphorylation():
    state = fst.state_from_string('Slt2_[(Y192)]-{p}')

    assert isinstance(state, sta.CovalentModificationState)
    assert str(state) == 'Slt2_[(Y192)]-{p}'
    assert str(state.substrate) == 'Slt2_[(Y192)]'
    assert state.modifier == sta.StateModifier.phosphor


# Reaction from string #
def test_rxncon_from_string_reaction_ppi():
    reaction = fst.reaction_from_string('Fus3_[CD]_ppi_Msg5_[n]')

    #assert reaction.classification_code == '2.1.1.1'

    assert reaction.subject.name == 'Fus3'
    assert reaction.subject.domain == 'CD'
    assert reaction.subject.subdomain is None
    assert reaction.subject.residue is None

    assert reaction.object.name == 'Msg5'
    assert reaction.object.domain == 'n'
    assert reaction.object.subdomain is None
    assert reaction.object.residue is None

    assert reaction.source is None
    assert str(reaction.product) == 'Fus3_[CD]--Msg5_[n]'


def test_rxncon_from_string_reaction_ipi():
    reaction = fst.reaction_from_string('A_[n]_ipi_A_[m]')

    #assert reaction.classification_code == '2.1.1.2'

    assert reaction.subject.name == 'A'
    assert reaction.subject.domain == 'n'
    assert reaction.subject.subdomain is None
    assert reaction.subject.residue is None

    assert reaction.object.name == 'A'
    assert reaction.object.domain == 'm'
    assert reaction.object.subdomain is None
    assert reaction.object.residue is None

    assert reaction.source is None
    assert str(reaction.product) == 'A_[n]--[m]'


def test_rxncon_from_string_reaction_i():
    reaction = fst.reaction_from_string('Pkc1_i_PS')

    #assert reaction.classification_code == '2.1.1.1'

    assert reaction.subject.name == 'Pkc1'
    assert reaction.subject.domain is None
    assert reaction.subject.subdomain is None
    assert reaction.subject.residue is None

    assert reaction.object.name == 'PS'
    assert reaction.object.domain is None
    assert reaction.object.subdomain is None
    assert reaction.object.residue is None

    assert reaction.source is None
    assert str(reaction.product) == 'Pkc1--PS'


def test_rxncon_from_string_reaction_bind():
    reaction = fst.reaction_from_string('Tec1_[n/TEA]_BIND_TCS')

    #assert reaction.classification_code == '2.1.1'

    assert reaction.subject.name == 'Tec1'
    assert reaction.subject.domain == 'n'
    assert reaction.subject.subdomain == 'TEA'
    assert reaction.subject.residue is None

    assert reaction.object.name == 'TCS'
    assert reaction.object.domain is None
    assert reaction.object.subdomain is None
    assert reaction.object.residue is None

    assert reaction.source is None
    assert str(reaction.product) == 'Tec1_[n/TEA]--TCS'


def test_rxncon_from_string_reaction_p_plus():
    reaction = fst.reaction_from_string('Fus3_[KD]_P+_Sst2_[(S539)]')

    #assert reaction.classification_code == '1.1.1'

    assert reaction.subject.name == 'Fus3'
    assert reaction.subject.domain == 'KD'
    assert reaction.subject.subdomain is None
    assert reaction.subject.residue is None

    assert reaction.object.name == 'Sst2'
    assert reaction.object.domain is None
    assert reaction.object.subdomain is None
    assert reaction.object.residue == 'S539'

    assert reaction.source is None
    assert str(reaction.product) == 'Sst2_[(S539)]-{p}'


def test_rxncon_from_string_reaction_p_minus():
    reaction = fst.reaction_from_string('Msg5_[PD]_P-_Slt2_[(Y192)]')

    #assert reaction.classification_code == '1.1.2'

    assert reaction.subject.name == 'Msg5'
    assert reaction.subject.domain == 'PD'
    assert reaction.subject.subdomain is None
    assert reaction.subject.residue is None

    assert reaction.object.name == 'Slt2'
    assert reaction.object.domain is None
    assert reaction.object.subdomain is None
    assert reaction.object.residue == 'Y192'

    assert str(reaction.source) == 'Slt2_[(Y192)]-{p}'
    assert reaction.product is None


def test_rxncon_from_string_reaction_gef():
    reaction = fst.reaction_from_string('Rom2_[DH]_GEF_Rho1_[GnP]')

    #assert reaction.classification_code == '1.1.1'

    assert reaction.subject.name == 'Rom2'
    assert reaction.subject.domain == 'DH'
    assert reaction.subject.subdomain is None
    assert reaction.subject.residue is None

    assert reaction.object.name == 'Rho1'
    assert reaction.object.domain == 'GnP'
    assert reaction.object.subdomain is None
    assert reaction.object.residue is None

    assert reaction.source is None
    assert str(reaction.product) == 'Rho1_[GnP]-{gtp}'


def test_rxncon_from_string_reaction_gap():
    reaction = fst.reaction_from_string('Lrg1_[GAP]_GAP_Rho1_[GnP]')

    #assert reaction.classification_code == '1.1.2'

    assert reaction.subject.name == 'Lrg1'
    assert reaction.subject.domain == 'GAP'
    assert reaction.subject.subdomain is None
    assert reaction.subject.residue is None

    assert reaction.object.name == 'Rho1'
    assert reaction.object.domain == 'GnP'
    assert reaction.object.subdomain is None
    assert reaction.object.residue is None

    assert str(reaction.source) == 'Rho1_[GnP]-{gtp}'
    assert reaction.product is None


def test_rxncon_from_string_reaction_ub_plus():
    reaction = fst.reaction_from_string('SCF_Ub+_Tec1')

    #assert reaction.classification_code == '1.1.1'

    assert reaction.subject.name == 'SCF'
    assert reaction.subject.domain is None
    assert reaction.subject.subdomain is None
    assert reaction.subject.residue is None

    assert reaction.object.name == 'Tec1'
    assert reaction.object.domain is None
    assert reaction.object.subdomain is None
    assert reaction.object.residue is None

    assert reaction.source is None
    assert str(reaction.product) == 'Tec1-{ub}'


def test_rxncon_from_string_reaction_ap():
    reaction = fst.reaction_from_string('Rck2_AP_Rck2_[Ser]')

    #assert reaction.classification_code == '1.1.1'

    assert reaction.subject.name == 'Rck2'
    assert reaction.subject.domain is None
    assert reaction.subject.subdomain is None
    assert reaction.subject.residue is None

    assert reaction.object.name == 'Rck2'
    assert reaction.object.domain == 'Ser'
    assert reaction.object.subdomain is None
    assert reaction.object.residue is None

    assert reaction.source is None
    assert str(reaction.product) == 'Rck2_[Ser]-{p}'


def test_rxncon_from_string_reaction_pt():
    reaction = fst.reaction_from_string('Sln1_[HK(H576)]_PT_Sln1_[RR(D1144)]')

    #assert reaction.classification_code == '1.1.3'

    assert reaction.subject.name == 'Sln1'
    assert reaction.subject.domain == 'HK'
    assert reaction.subject.subdomain is None
    assert reaction.subject.residue == 'H576'

    assert reaction.object.name == 'Sln1'
    assert reaction.object.domain == 'RR'
    assert reaction.object.subdomain is None
    assert reaction.object.residue == 'D1144'

    assert str(reaction.source) == 'Sln1_[HK(H576)]-{p}'
    assert str(reaction.product) == 'Sln1_[RR(D1144)]-{p}'


def test_rxncon_from_string_reaction_deg():
    #todo: does this makes sense??
    reaction = fst.reaction_from_string('Bar1_[PepD]_DEG_MFalpha_[(L6-K7)]')

    #assert reaction.classification_code == '3.2.2'

    assert reaction.subject.name == 'Bar1'
    assert reaction.subject.domain == 'PepD'
    assert reaction.subject.subdomain is None
    assert reaction.subject.residue is None

    assert reaction.object.name == 'MFalpha'
    assert reaction.object.domain is None
    assert reaction.object.subdomain is None
    assert reaction.object.residue == 'L6-K7'

    #todo: SynthesisDegedationStates are ComponentenStates
    #assert str(reaction.source) == 'MFalpha_[(L6-K7)]'
    assert str(reaction.source) == 'MFalpha'
    assert reaction.product is None


def test_rxncon_from_string_reaction_cut():
    reaction = fst.reaction_from_string('Yps1_CUT_Msb2_[HMH/CD]')

    #assert reaction.classification_code == '1.2.1'

    assert reaction.subject.name == 'Yps1'
    assert reaction.subject.domain is None
    assert reaction.subject.subdomain is None
    assert reaction.subject.residue is None

    assert reaction.object.name == 'Msb2'
    assert reaction.object.domain == 'HMH'
    assert reaction.object.subdomain == 'CD'
    assert reaction.object.residue is None

    assert reaction.source is None
    assert str(reaction.product) == 'Msb2_[HMH/CD]-{truncated}'


