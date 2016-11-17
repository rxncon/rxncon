from rxncon.input.quick.quick import Quick
from rxncon.simulation.rule_based.rule_based_model import rule_based_model_from_rxncon, calc_positive_solutions
from rxncon.core.state import state_from_str

def test_simple_system():
    rbm = rule_based_model_from_rxncon(Quick("""A_[b]_ppi+_B_[a]; ! A_[(r)]-{p}
                                             A_[b]_ppi-_B_[a]
                                             C_p+_A_[(r)]
                                             D_p-_A_[(r)]""").rxncon_system)

    print()
    for rule in rbm.rules:
        print(rule)


def test_calc_positive_solutions():
    rxncon_sys = Quick("""A_[b]_ppi+_B_[a]; ! A_[(r)]-{p}
                          E_[x]_ppi+_B_[a]
                          C_p+_A_[(r)]
                          D_ub+_A_[(r)]""").rxncon_system

    for x in rule_based_model_from_rxncon(rxncon_sys).rules:
        print(x)



def test_2():
    rbm = rule_based_model_from_rxncon(Quick('''A_ppi+_C
                                             C_ppi+_D
                                             B_ppi+_E
                                             B_ppi+_F
                                             A_ppi+_B; x <comp1>
                                             <comp1>; OR <comp1C1>
                                             <comp1>; OR <comp2C1>
                                             <comp1C1>; AND A--C
                                             <comp1C1>; AND C--D
                                             <comp2C1>; AND B--F
                                             <comp2C1>; AND B--E''').rxncon_system)

    print()
    for rule in rbm.rules:
        print(rule)

def test_3():
    rbm = rule_based_model_from_rxncon(Quick('''A_ppi+_C
                                             C_ppi+_D
                                             B_ppi+_E
                                             A_ppi+_B; x <comp1>
                                             <comp1>; AND <comp1C1>
                                             <comp1>; AND <comp2C1>
                                             <comp1C1>; OR A--C
                                             <comp1C1>; OR C--D
                                             <comp2C1>; OR A--C
                                             <comp2C1>; OR B--E''').rxncon_system)

    print()
    for rule in rbm.rules:
        print(rule)

def test_insulin():
    rbm = rule_based_model_from_rxncon(Quick("""<IR-empty>; AND IR@0_[IRBD]--IR@2_[IRBD]; AND IR@0_[lig]--0; AND IR@2_[lig]--0
                                                        <IR-IR-active>; AND <IR-phos>; AND <IR0-IR1-Insulin2>
                                                        <IR-phos>; AND IR@0_[TK(Y1158)]-{P}; AND IR@0_[TK(Y1162)]-{P}; AND IR@0_[TK(Y1163)]-{P}
                                                        <IR0-IR1-Insulin2>; OR <IR0-insulin2>; OR <IR1-insulin2>
                                                        <IR0-insulin2>; AND IR@0_[IRBD]--IR@1_[IRBD]; AND IR@0_[lig]--insulin@2_[IR]
                                                        <IR1-insulin2>; AND IR@0_[IRBD]--IR@1_[IRBD]; AND IR@1_[lig]--insulin@3_[IR]
                                                        <IR-active>; AND <IR-phos>; AND <IR0-IR1-Insulin2>
                                                        IR_p+_IR_[TK(Y1158)]; ! <IR0-IR1-Insulin2>
                                                        IR_[IRBD]_ppi-_IR_[IRBD]
                                                        IR_[lig]_i+_insulin_[IR]; ! <IR-empty>
                                                        IR_p+_IR_[TK(Y1162)]; ! <IR0-IR1-Insulin2>
                                                        IR_p+_IR_[TK(Y1163)]; ! <IR0-IR1-Insulin2>
                                                        IR_[IRBD]_ppi+_IR_[IRBD]
                                                        IR_[lig]_i-_insulin_[IR]
                                                        IR_ap+_IR_[JM(Y972)]; ! <IR-IR-active>
                                                        IR_[JMY972]_ppi+_IRS_[PTB]; ! IR@0_[JM(Y972)]-{P}; ! IRS@1_[PH]--Phospholipids@2_[IRS]
                                                        IR_[JMY972]_ppi-_IRS_[PTB]
                                                        IRS_[PH]_i+_Phospholipids_[IRS]
                                                        IRS_[PH]_i-_Phospholipids_[IRS]
                                                        IR_p+_IRS_[(Y)]; ! IR@0_[JMY972]--IRS@1_[PTB]; ! <IR-active>
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
                                                        IRS_[Y]_ppi-_PI3K_[SH2]""").rxncon_system)




    print()
    for x in rbm.rules:
        print()
        print(x)