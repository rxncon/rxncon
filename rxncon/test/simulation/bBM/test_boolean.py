from collections import namedtuple

from rxncon.venntastic.sets import venn_from_str
from rxncon.input.quick.quick import Quick
from rxncon.core.reaction import reaction_from_str
from rxncon.core.state import state_from_str
from rxncon.core.spec import mol_spec_from_str
from rxncon.simulation.bBM.boolean_model import boolean_model_from_rxncon, boolnet_str_from_boolean_model, ReactionTarget, \
    StateTarget, ComponentStateTarget


def target_from_str(target_str):
    try:
        return ReactionTarget(reaction_from_str(target_str))
    except SyntaxError:
        pass

    try:
        return StateTarget(state_from_str(target_str))
    except SyntaxError:
        pass

    try:
        return ComponentStateTarget(mol_spec_from_str(target_str))
    except SyntaxError:
        raise SyntaxError('Could not parse target str {}'.format(target_str))



def test_simple_system():
    boolean_model = boolean_model_from_rxncon(Quick("""A_[b]_ppi+_B_[a]; ! A_[(r)]-{p}
                                                    A_[b]_ppi-_B_[a]
                                                    C_p+_A_[(r)]
                                                    D_p-_A_[(r)]""").rxncon_system)

    A = '(( A_[b]--0 | A_[b]--B_[a] ) & ( A_[(r)]-{0} | A_[(r)]-{p} ))'
    B = '( B_[a]--0 | A_[b]--B_[a] )'
    C = 'C'
    D = 'D'

    expected_rules = {
        target_from_str('A_[b]_ppi+_B_[a]'): venn_from_str('{0} & {1} & A_[(r)]-{{p}}'.format(A, B), target_from_str),
        target_from_str('A_[b]_ppi-_B_[a]'): venn_from_str('{0} & {1}'.format(A, B), target_from_str),
        target_from_str('C_p+_A_[(r)]'):     venn_from_str('{0} & {1}'.format(A, C), target_from_str),
        target_from_str('D_p-_A_[(r)]'):     venn_from_str('{0} & {1}'.format(A, D), target_from_str),
        target_from_str('A_[(r)]-{p}'):      venn_from_str('{0} & (( C_p+_A_[(r)] & A_[(r)]-{{0}} ) | ( A_[(r)]-{{p}} & ~( D_p-_A_[(r)] & A_[(r)]-{{p}} )))'.format(A), target_from_str),
        target_from_str('A_[(r)]-{0}'):      venn_from_str('{0} & (( D_p-_A_[(r)] & A_[(r)]-{{p}} ) | ( A_[(r)]-{{0}} & ~( C_p+_A_[(r)] & A_[(r)]-{{0}} )))'.format(A), target_from_str),
        target_from_str('A_[b]--0'):         venn_from_str('{0} & (( A_[b]_ppi-_B_[a] & A_[b]--B_[a] ) | ( A_[b]--0 & ~( A_[b]_ppi+_B_[a] & A_[b]--0 & B_[a]--0 )))'.format(A), target_from_str),
        target_from_str('A_[b]--B_[a]'):     venn_from_str('{0} & {1} & (( A_[b]_ppi+_B_[a] & A_[b]--0 & B_[a]--0 ) | ( A_[b]--B_[a] & ~( A_[b]_ppi-_B_[a] & A_[b]--B_[a] )))'.format(A, B), target_from_str),
        target_from_str('B_[a]--0'):         venn_from_str('{0} & (( A_[b]_ppi-_B_[a] & A_[b]--B_[a] ) | ( B_[a]--0 & ~( A_[b]_ppi+_B_[a] & A_[b]--0 & B_[a]--0 )))'.format(B), target_from_str),
        target_from_str('C'):                venn_from_str('{0}'.format(C), target_from_str),
        target_from_str('D'):                venn_from_str('{0}'.format(D), target_from_str),
    }

    assert len(boolean_model.update_rules) == len(expected_rules)

    for update_rule in boolean_model.update_rules:
        assert update_rule.factor.is_equivalent_to(expected_rules[update_rule.target])


def test_insulin_system():
    rxncon_sys = Quick("""IR_[IRBD]_ppi+_IR_[IRBD]
    IR_[IRBD]_ppi-_IR_[IRBD]
    IR_[lig]_i+_insulin_[IR]; ! <IR-empty>
    IR_[lig]_i-_insulin_[IR]
    IR_p+_IR_[TK(Y1158)]; ! <IR0-IR1-Insulin2>
    IR_p+_IR_[TK(Y1162)]; ! <IR0-IR1-Insulin2>
    IR_p+_IR_[TK(Y1163)]; ! <IR0-IR1-Insulin2>
    IR_ap+_IR_[JM(Y972)]; ! <IR-IR-active>
    IR_[JMY972]_ppi+_IRS_[PTB]; ! IR_[JM(Y972)]-{P}; ! IRS_[PH]--Phospholipids_[IRS]
    IR_[JMY972]_ppi-_IRS_[PTB]
    IRS_[PH]_i+_Phospholipids_[IRS]
    IRS_[PH]_i-_Phospholipids_[IRS]
    IR_p+_IRS_[(Y)]; ! IR_[JMY972]--IRS_[PTB]; ! <IR-active>
    Grb2_[SOS]_ppi+_SOS_[Grb2]
    Grb2_[SOS]_ppi-_SOS_[Grb2]
    Grb2_[SH2]_ppi+_IRS_[Y]; ! IRS_[(Y)]-{P}
    Grb2_[SH2]_ppi-_IRS_[Y]
    IR_[JMY972]_ppi+_Shc_[PTB]; ! IR_[JM(Y972)]-{P}
    IR_[JMY972]_ppi-_Shc_[PTB]
    IR_p+_Shc_[(YY239240)]; ! IR_[JMY972]--Shc_[PTB]; ! <IR-active>
    IR_p+_Shc_[(Y317)]; ! IR_[JMY972]--Shc_[PTB]; ! <IR-active>
    Grb2_[SH2]_ppi+_Shc_[YY]; ! Shc_[(YY239240)]-{P}
    Grb2_[SH2]_ppi-_Shc_[YY]
    Grb2_[SH2]_ppi+_Shc_[Y]; ! Shc_[(Y317)]-{P}
    Grb2_[SH2]_ppi-_Shc_[Y]
    IRS_[Y]_ppi+_PI3K_[SH2]; ! IRS_[(Y)]-{P}
    IRS_[Y]_ppi-_PI3K_[SH2]

    <IR-empty>; AND IR@0_[IRBD]--IR@2_[IRBD]; AND IR@0_[lig]--0; AND IR@2_[lig]--0

    <IR-IR-active>; AND <IR-phos>; AND <IR0-IR1-Insulin2>
    <IR-phos>; AND IR_[TK(Y1158)]-{P}; AND IR_[TK(Y1162)]-{P}; AND IR_[TK(Y1163)]-{P}

    <IR0-IR1-Insulin2>; OR <IR0-insulin2>; OR <IR1-insulin2>
    <IR0-insulin2>; AND IR@0_[IRBD]--IR@1_[IRBD]; AND IR@0_[lig]--insulin@2_[IR]
    <IR1-insulin2>; AND IR@0_[IRBD]--IR@1_[IRBD]; AND IR@1_[lig]--insulin@2_[IR]

    <IR-active>; AND <IR-phos>; AND <IR0-IR2-Insulin3>
    <IR0-IR2-Insulin3>; OR <IR0-insulin3>; OR <IR2-insulin3>
    <IR0-insulin3>; AND IR@0_[IRBD]--IR@2_[IRBD]; AND IR@0_[lig]--insulin@3_[IR]
    <IR2-insulin3>; AND IR@0_[IRBD]--IR@2_[IRBD]; AND IR@2_[lig]--insulin@3_[IR]
    """)
    boolmodel = boolean_model_from_rxncon(rxncon_sys.rxncon_system)
    print('hallo')

BooleanRuleTestCase = namedtuple('BooleanRuleTestCase', ['boolean_string', 'boolean_factor'])

