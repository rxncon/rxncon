from rxncon.input.quick.quick import Quick
from rxncon.simulation.boolean.boolean_model import boolean_model_from_rxncon
from rxncon.test.simulation.boolean.utils import target_from_str
from rxncon.venntastic.sets import venn_from_str


def test_insulin_no_smoothing() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""IR_[IRBD]_ppi+_IR_[IRBD]
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
                                                    """).rxncon_system)

    # Component expressions
    IR            = '(( IR_[IRBD]--0 | IR_[IRBD]--IR_[IRBD] ) & ( IR_[lig]--0 | IR_[lig]--insulin_[IR] ) & ' \
                    ' ( IR_[TK(Y1158)]-{0} | IR_[TK(Y1158)]-{p} ) & ( IR_[TK(Y1162)]-{0} | IR_[TK(Y1162)]-{p} ) & ' \
                    ' ( IR_[TK(Y1163)]-{0} | IR_[TK(Y1163)]-{p} ) & ( IR_[JM(Y972)]-{0} | IR_[JM(Y972)]-{p} ) & ' \
                    ' ( IR_[JMY972]--0 | IR_[JMY972]--IRS_[PTB] | IR_[JMY972]--Shc_[PTB] ))'
    insulin       = '( insulin_[IR]--0 | IR_[lig]--insulin_[IR] )'
    IRS           = '(( IRS_[PTB]--0 | IR_[JMY972]--IRS_[PTB] ) & ( IRS_[PH]--0 | IRS_[PH]--Phospholipids_[IRS] ) & ' \
                    ' ( IRS_[(Y)]-{0} | IRS_[(Y)]-{p} ) & ( IRS_[Y]--0 | Grb2_[SH2]--IRS_[Y] | IRS_[Y]--PI3K_[SH2] ))'
    Phospholipids = '( Phospholipids_[IRS]--0 | IRS_[PH]--Phospholipids_[IRS] )'
    Grb2          = '(( Grb2_[SOS]--0 | Grb2_[SOS]--SOS_[Grb2] ) & ( Grb2_[SH2]--0 | Grb2_[SH2]--IRS_[Y] | Grb2_[SH2]--Shc_[Y] | Grb2_[SH2]--Shc_[YY] ))'
    SOS           = '( SOS_[Grb2]--0 | Grb2_[SOS]--SOS_[Grb2] )'
    Shc           = '(( Shc_[PTB]--0 | IR_[JMY972]--Shc_[PTB] ) & ( Shc_[(YY239240)]-{0} | Shc_[(YY239240)]-{p} ) & ' \
                    ' ( Shc_[(Y317)]-{0} | Shc_[(Y317)]-{p} ) & ( Shc_[Y]--0 | Grb2_[SH2]--Shc_[Y] ) & ( Shc_[YY]--0 | Grb2_[SH2]--Shc_[YY] ))'
    PI3K          = '( PI3K_[SH2]--0 | IRS_[Y]--PI3K_[SH2] )'

    expected_rules = {
        # Reaction updates
        'IR_[IRBD]_ppi+_IR_[IRBD]':        '{0}'.format(IR),
        'IR_[IRBD]_ppi-_IR_[IRBD]':        '{0}'.format(IR),
        'IR_[lig]_i+_insulin_[IR]':        '{0} & {1} & IR_[IRBD]--IR_[IRBD] & IR_[lig]--0'.format(IR, insulin),
        'IR_[lig]_i-_insulin_[IR]':        '{0} & {1}'.format(IR, insulin),
        'IR_p+_IR_[TK(Y1158)]':            '{0} & IR_[IRBD]--IR_[IRBD] & IR_[lig]--insulin_[IR]'.format(IR),
        'IR_p+_IR_[TK(Y1162)]':            '{0} & IR_[IRBD]--IR_[IRBD] & IR_[lig]--insulin_[IR]'.format(IR),
        'IR_p+_IR_[TK(Y1163)]':            '{0} & IR_[IRBD]--IR_[IRBD] & IR_[lig]--insulin_[IR]'.format(IR),
        'IR_ap+_IR_[JM(Y972)]':            '{0} & IR_[TK(Y1158)]-{{p}} & IR_[TK(Y1162)]-{{p}} & IR_[TK(Y1163)]-{{p}} & IR_[IRBD]--IR_[IRBD] & IR_[lig]--insulin_[IR]'.format(IR),
        'IR_[JMY972]_ppi+_IRS_[PTB]':      '{0} & {1} & IR_[JM(Y972)]-{{p}} & IRS_[PH]--Phospholipids_[IRS]'.format(IR, IRS),
        'IR_[JMY972]_ppi-_IRS_[PTB]':      '{0} & {1}'.format(IR, IRS),
        'IRS_[PH]_i+_Phospholipids_[IRS]': '{0} & {1}'.format(IRS, Phospholipids),
        'IRS_[PH]_i-_Phospholipids_[IRS]': '{0} & {1}'.format(IRS, Phospholipids),
        'IR_p+_IRS_[(Y)]':                 '{0} & {1} & IR_[JMY972]--IRS_[PTB] & IR_[TK(Y1158)]-{{p}} & IR_[TK(Y1162)]-{{p}} & IR_[TK(Y1163)]-{{p}} & IR_[IRBD]--IR_[IRBD] & IR_[lig]--insulin_[IR]'.format(IR, IRS),
        'Grb2_[SOS]_ppi+_SOS_[Grb2]':      '{0} & {1}'.format(Grb2, SOS),
        'Grb2_[SOS]_ppi-_SOS_[Grb2]':      '{0} & {1}'.format(Grb2, SOS),
        'Grb2_[SH2]_ppi+_IRS_[Y]':         '{0} & {1} & IRS_[(Y)]-{{p}}'.format(Grb2, IRS),
        'Grb2_[SH2]_ppi-_IRS_[Y]':         '{0} & {1}'.format(Grb2, IRS),
        'IR_[JMY972]_ppi+_Shc_[PTB]':      '{0} & {1} & IR_[JM(Y972)]-{{p}}'.format(IR, Shc),
        'IR_[JMY972]_ppi-_Shc_[PTB]':      '{0} & {1}'.format(IR, Shc),
        'IR_p+_Shc_[(YY239240)]':          '{0} & {1} & IR_[JMY972]--Shc_[PTB] & IR_[TK(Y1158)]-{{p}} & IR_[TK(Y1162)]-{{p}} & IR_[TK(Y1163)]-{{p}} & IR_[IRBD]--IR_[IRBD] & IR_[lig]--insulin_[IR]'.format(IR, Shc),
        'IR_p+_Shc_[(Y317)]':              '{0} & {1} & IR_[JMY972]--Shc_[PTB] & IR_[TK(Y1158)]-{{p}} & IR_[TK(Y1162)]-{{p}} & IR_[TK(Y1163)]-{{p}} & IR_[IRBD]--IR_[IRBD] & IR_[lig]--insulin_[IR]'.format(IR, Shc),
        'Grb2_[SH2]_ppi+_Shc_[YY]':        '{0} & {1} & Shc_[(YY239240)]-{{p}}'.format(Grb2, Shc),
        'Grb2_[SH2]_ppi-_Shc_[YY]':        '{0} & {1}'.format(Grb2, Shc),
        'Grb2_[SH2]_ppi+_Shc_[Y]':         '{0} & {1} & Shc_[(Y317)]-{{p}}'.format(Grb2, Shc),
        'Grb2_[SH2]_ppi-_Shc_[Y]':         '{0} & {1}'.format(Grb2, Shc),
        'IRS_[Y]_ppi+_PI3K_[SH2]':         '{0} & {1} & IRS_[(Y)]-{{p}}'.format(IRS, PI3K),
        'IRS_[Y]_ppi-_PI3K_[SH2]':         '{0} & {1}'.format(IRS, PI3K),
        # State updates
        'IR_[IRBD]--0':                    '{0} & (( IR_[IRBD]_ppi-_IR_[IRBD] & IR_[IRBD]--IR_[IRBD] ) | ( IR_[IRBD]--0 & ~( IR_[IRBD]_ppi+_IR_[IRBD] & IR_[IRBD]--0 )))'.format(IR),
        'IR_[IRBD]--IR_[IRBD]':            '{0} & (( IR_[IRBD]_ppi+_IR_[IRBD] & IR_[IRBD]--0 ) | ( IR_[IRBD]--IR_[IRBD] & ~( IR_[IRBD]_ppi-_IR_[IRBD] & IR_[IRBD]--IR_[IRBD] )))'.format(IR),
        'IR_[lig]--0':                     '{0} & (( IR_[lig]_i-_insulin_[IR] & IR_[lig]--insulin_[IR] ) | ( IR_[lig]--0 & ~( IR_[lig]_i+_insulin_[IR] & IR_[lig]--0 & insulin_[IR]--0 )))'.format(IR),
        'IR_[lig]--insulin_[IR]':          '{0} & {1} & (( IR_[lig]_i+_insulin_[IR] & IR_[lig]--0 & insulin_[IR]--0 ) | ( IR_[lig]--insulin_[IR] & ~( IR_[lig]_i-_insulin_[IR] & IR_[lig]--insulin_[IR] )))'.format(IR, insulin),
        'insulin_[IR]--0':                 '{0} & (( IR_[lig]_i-_insulin_[IR] & IR_[lig]--insulin_[IR] ) | ( insulin_[IR]--0 & ~( IR_[lig]_i+_insulin_[IR] & IR_[lig]--0 & insulin_[IR]--0 )))'.format(insulin),
        'IR_[TK(Y1158)]-{0}':              '{0} & ( IR_[TK(Y1158)]-{{0}} & ~( IR_p+_IR_[TK(Y1158)] & IR_[TK(Y1158)]-{{0}} ))'.format(IR),
        'IR_[TK(Y1162)]-{0}':              '{0} & ( IR_[TK(Y1162)]-{{0}} & ~( IR_p+_IR_[TK(Y1162)] & IR_[TK(Y1162)]-{{0}} ))'.format(IR),
        'IR_[TK(Y1163)]-{0}':              '{0} & ( IR_[TK(Y1163)]-{{0}} & ~( IR_p+_IR_[TK(Y1163)] & IR_[TK(Y1163)]-{{0}} ))'.format(IR),
        'IR_[JM(Y972)]-{0}':               '{0} & ( IR_[JM(Y972)]-{{0}} & ~( IR_ap+_IR_[JM(Y972)] & IR_[JM(Y972)]-{{0}} ))'.format(IR),
        'IR_[TK(Y1158)]-{p}':              '{0} & ( IR_[TK(Y1158)]-{{p}} | ( IR_p+_IR_[TK(Y1158)] & IR_[TK(Y1158)]-{{0}} ))'.format(IR),
        'IR_[TK(Y1162)]-{p}':              '{0} & ( IR_[TK(Y1162)]-{{p}} | ( IR_p+_IR_[TK(Y1162)] & IR_[TK(Y1162)]-{{0}} ))'.format(IR),
        'IR_[TK(Y1163)]-{p}':              '{0} & ( IR_[TK(Y1163)]-{{p}} | ( IR_p+_IR_[TK(Y1163)] & IR_[TK(Y1163)]-{{0}} ))'.format(IR),
        'IR_[JM(Y972)]-{p}':               '{0} & ( IR_[JM(Y972)]-{{p}} | ( IR_ap+_IR_[JM(Y972)] & IR_[JM(Y972)]-{{0}} ))'.format(IR),
        'IR_[JMY972]--0':                  '{0} & (( IR_[JMY972]_ppi-_IRS_[PTB] & IR_[JMY972]--IRS_[PTB] ) | ( IR_[JMY972]_ppi-_Shc_[PTB] & IR_[JMY972]--Shc_[PTB] ) | ( IR_[JMY972]--0 & ~( IR_[JMY972]_ppi+_IRS_[PTB] & IR_[JMY972]--0 & IRS_[PTB]--0 ) & ~( IR_[JMY972]_ppi+_Shc_[PTB] & IR_[JMY972]--0 & Shc_[PTB]--0 )))'.format(IR),
        'IRS_[PTB]--IR_[JMY972]':          '{0} & {1} & (( IR_[JMY972]_ppi+_IRS_[PTB] & IR_[JMY972]--0 & IRS_[PTB]--0 ) | ( IR_[JMY972]--IRS_[PTB] & ~( IR_[JMY972]_ppi-_IRS_[PTB] & IR_[JMY972]--IRS_[PTB] )))'.format(IR, IRS),
        'IRS_[PTB]--0':                    '{0} & (( IR_[JMY972]_ppi-_IRS_[PTB] & IR_[JMY972]--IRS_[PTB] ) | ( IRS_[PTB]--0 & ~( IR_[JMY972]_ppi+_IRS_[PTB] & IR_[JMY972]--0 & IRS_[PTB]--0 )))'.format(IRS),
        'IR_[JMY972]--Shc_[PTB]':          '{0} & {1} & (( IR_[JMY972]_ppi+_Shc_[PTB] & IR_[JMY972]--0 & Shc_[PTB]--0 ) | ( IR_[JMY972]--Shc_[PTB] & ~( IR_[JMY972]_ppi-_Shc_[PTB] & IR_[JMY972]--Shc_[PTB] )))'.format(IR, Shc),
        'Shc_[PTB]--0':                    '{0} & (( IR_[JMY972]_ppi-_Shc_[PTB] & IR_[JMY972]--Shc_[PTB] ) | ( Shc_[PTB]--0 & ~( IR_[JMY972]_ppi+_Shc_[PTB] & IR_[JMY972]--0 & Shc_[PTB]--0 )))'.format(Shc),
        'IRS_[PH]--0':                     '{0} & (( IRS_[PH]_i-_Phospholipids_[IRS] & IRS_[PH]--Phospholipids_[IRS] ) | ( IRS_[PH]--0 & ~( IRS_[PH]_i+_Phospholipids_[IRS] & IRS_[PH]--0 & Phospholipids_[IRS]--0 )))'.format(IRS),
        'IRS_[PH]--Phospholipids_[IRS]':   '{0} & {1} & (( IRS_[PH]_i+_Phospholipids_[IRS] & IRS_[PH]--0 & Phospholipids_[IRS]--0 ) | ( IRS_[PH]--Phospholipids_[IRS] & ~( IRS_[PH]_i-_Phospholipids_[IRS] & IRS_[PH]--Phospholipids_[IRS] )))'.format(IRS, Phospholipids),
        'Phospholipids_[IRS]--0':          '{0} & (( IRS_[PH]_i-_Phospholipids_[IRS] & IRS_[PH]--Phospholipids_[IRS] ) | ( Phospholipids_[IRS]--0 & ~( IRS_[PH]_i+_Phospholipids_[IRS] & IRS_[PH]--0 & Phospholipids_[IRS]--0 )))'.format(Phospholipids),
        'IRS_[(Y)]-{0}':                   '{0} & ( IRS_[(Y)]-{{0}} & ~( IR_p+_IRS_[(Y)] & IRS_[(Y)]-{{0}} ))'.format(IRS),
        'IRS_[(Y)]-{p}':                   '{0} & ( IRS_[(Y)]-{{p}} | ( IR_p+_IRS_[(Y)] & IRS_[(Y)]-{{0}} ))'.format(IRS),
        'Grb2_[SOS]--0':                   '{0} & (( Grb2_[SOS]_ppi-_SOS_[Grb2] & Grb2_[SOS]--SOS_[Grb2] ) | ( Grb2_[SOS]--0 & ~( Grb2_[SOS]_ppi+_SOS_[Grb2] & Grb2_[SOS]--0 & SOS_[Grb2]--0 )))'.format(Grb2),
        'Grb2_[SOS]--SOS_[Grb2]':          '{0} & {1} & (( Grb2_[SOS]_ppi+_SOS_[Grb2] & Grb2_[SOS]--0 & SOS_[Grb2]--0 ) | ( Grb2_[SOS]--SOS_[Grb2] & ~( Grb2_[SOS]_ppi-_SOS_[Grb2] & Grb2_[SOS]--SOS_[Grb2] )))'.format(Grb2, SOS),
        'SOS_[Grb2]--0':                   '{0} & (( Grb2_[SOS]_ppi-_SOS_[Grb2] & Grb2_[SOS]--SOS_[Grb2] ) | ( SOS_[Grb2]--0 & ~( Grb2_[SOS]_ppi+_SOS_[Grb2] & Grb2_[SOS]--0 & SOS_[Grb2]--0 )))'.format(SOS),
        'IRS_[Y]--0':                      '{0} & (( IRS_[Y]_ppi-_PI3K_[SH2] & IRS_[Y]--PI3K_[SH2] ) | ( Grb2_[SH2]_ppi-_IRS_[Y] & Grb2_[SH2]--IRS_[Y] ) | ( IRS_[Y]--0 & ~( IRS_[Y]_ppi+_PI3K_[SH2] & IRS_[Y]--0 & PI3K_[SH2]--0 ) & ~( Grb2_[SH2]_ppi+_IRS_[Y] & IRS_[Y]--0 & Grb2_[SH2]--0 )))'.format(IRS),
        'IRS_[Y]--PI3K_[SH2]':             '{0} & {1} & (( IRS_[Y]_ppi+_PI3K_[SH2] & IRS_[Y]--0 & PI3K_[SH2]--0 ) | ( IRS_[Y]--PI3K_[SH2] & ~( IRS_[Y]_ppi-_PI3K_[SH2] & IRS_[Y]--PI3K_[SH2] )))'.format(IRS, PI3K),
        'PI3K_[SH2]--0':                   '{0} & (( IRS_[Y]_ppi-_PI3K_[SH2] & IRS_[Y]--PI3K_[SH2] ) | ( PI3K_[SH2]--0 & ~( IRS_[Y]_ppi+_PI3K_[SH2] & IRS_[Y]--0 & PI3K_[SH2]--0 )))'.format(PI3K),
        'Grb2_[SH2]--IRS_[Y]':             '{0} & {1} & (( Grb2_[SH2]_ppi+_IRS_[Y] & Grb2_[SH2]--0 & IRS_[Y]--0 ) | ( Grb2_[SH2]--IRS_[Y] & ~( Grb2_[SH2]_ppi-_IRS_[Y] & Grb2_[SH2]--IRS_[Y] )))'.format(Grb2, IRS),
        'Grb2_[SH2]--0':                   '{0} & (( Grb2_[SH2]_ppi-_IRS_[Y] & Grb2_[SH2]--IRS_[Y] ) | ( Grb2_[SH2]_ppi-_Shc_[YY] & Grb2_[SH2]--Shc_[YY] ) | ( Grb2_[SH2]_ppi-_Shc_[Y] & Grb2_[SH2]--Shc_[Y] ) | ( Grb2_[SH2]--0 & ~( Grb2_[SH2]_ppi+_IRS_[Y] & Grb2_[SH2]--0 & IRS_[Y]--0 ) & ~( Grb2_[SH2]_ppi+_Shc_[YY] & Grb2_[SH2]--0 & Shc_[YY]--0 ) & ~( Grb2_[SH2]_ppi+_Shc_[Y] & Grb2_[SH2]--0 & Shc_[Y]--0 )))'.format(Grb2),
        'Grb2_[SH2]--Shc_[YY]':            '{0} & {1} & (( Grb2_[SH2]_ppi+_Shc_[YY] & Grb2_[SH2]--0 & Shc_[YY]--0 ) | ( Grb2_[SH2]--Shc_[YY] & ~( Grb2_[SH2]_ppi-_Shc_[YY] & Grb2_[SH2]--Shc_[YY] )))'.format(Grb2, Shc),
        'Shc_[YY]--0':                     '{0} & (( Grb2_[SH2]_ppi-_Shc_[YY] & Grb2_[SH2]--Shc_[YY] ) | ( Shc_[YY]--0 & ~( Grb2_[SH2]_ppi+_Shc_[YY] & Grb2_[SH2]--0 & Shc_[YY]--0 )))'.format(Shc),
        'Grb2_[SH2]--Shc_[Y]':             '{0} & {1} & (( Grb2_[SH2]_ppi+_Shc_[Y] & Grb2_[SH2]--0 & Shc_[Y]--0 ) | ( Grb2_[SH2]--Shc_[Y] & ~( Grb2_[SH2]_ppi-_Shc_[Y] & Grb2_[SH2]--Shc_[Y] )))'.format(Grb2, Shc),
        'Shc_[Y]--0':                      '{0} & (( Grb2_[SH2]_ppi-_Shc_[Y] & Grb2_[SH2]--Shc_[Y] ) | ( Shc_[Y]--0 & ~( Grb2_[SH2]_ppi+_Shc_[Y] & Grb2_[SH2]--0 & Shc_[Y]--0 )))'.format(Shc),
        'Shc_[(YY239240)]-{0}':            '{0} & ( Shc_[(YY239240)]-{{0}} & ~( IR_p+_Shc_[(YY239240)] & Shc_[(YY239240)]-{{0}} ))'.format(Shc),
        'Shc_[(YY239240)]-{p}':            '{0} & ( Shc_[(YY239240)]-{{p}} | ( IR_p+_Shc_[(YY239240)] & Shc_[(YY239240)]-{{0}} ))'.format(Shc),
        'Shc_[(Y317)]-{0}':                '{0} & ( Shc_[(Y317)]-{{0}} & ~( IR_p+_Shc_[(Y317)] & Shc_[(Y317)]-{{0}} ))'.format(Shc),
        'Shc_[(Y317)]-{p}':                '{0} & ( Shc_[(Y317)]-{{p}} | ( IR_p+_Shc_[(Y317)] & Shc_[(Y317)]-{{0}} ))'.format(Shc),
    }

    assert len(boolean_model.update_rules) == len(expected_rules)

    for update_rule in boolean_model.update_rules:
        assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)], target_from_str))
