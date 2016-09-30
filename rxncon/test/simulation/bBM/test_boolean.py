import pytest
import rxncon.input.quick.quick as quick
import rxncon.simulation.bBM.boolean_model as bm


def test_simple_system():
    rxncon_sys = quick.Quick("""A_[b]_ppi+_B_[a]; ! A_[(r)]-{p}
                                A_[b]_ppi-_B_[a]
                                C_p+_A_[(r)]
                                D_p-_A_[(r)]""")
    boolmodel = bm.boolean_model_from_rxncon(rxncon_sys.rxncon_system)
    for rule in boolmodel.update_rules:
        print(rule)
    boolmodel

def test_insulin_system():
    rxncon_sys = quick.Quick("""IR_[IRBD]_ppi+_IR_[IRBD]
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
    boolmodel = bm.boolean_model_from_rxncon(rxncon_sys.rxncon_system)

def test_rule_from_str():
    test = bm.boolean_rule_from_str('a--b, (a-{p} | (a_ppi+_b &c--d))')
    test