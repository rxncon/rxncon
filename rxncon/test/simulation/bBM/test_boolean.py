import pytest
from collections import namedtuple
import rxncon.input.quick.quick as quick
import rxncon.simulation.bBM.boolean_model as bm
from rxncon.core.reaction import reaction_from_string
from rxncon.core.state import state_from_string
from rxncon.venntastic.sets import Union, Intersection, ValueSet, Complement


def test_simple_system():
    rxncon_sys = quick.Quick("""A_[b]_ppi+_B_[a]; ! A_[(r)]-{p}
                                A_[b]_ppi-_B_[a]
                                C_p+_A_[(r)]
                                D_p-_A_[(r)]""")
    boolmodel = bm.boolean_model_from_rxncon(rxncon_sys.rxncon_system)

    print(bm.boolnet_str_from_boolean_model(boolmodel))

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
    print('hallo')

BooleanRuleTestCase = namedtuple('BooleanRuleTestCase', ['boolean_string', 'boolean_factor'])

@pytest.fixture
def boolean_test_cases():
    return [
        BooleanRuleTestCase('<a--a | <b_ppi+_a | <<c--a&d--e> | <e_ppi-_a&f--r> >>>',
                            bm.Union(bm.ValueSet(bm.StateTarget(state_from_string('a--a'))), bm.Union(bm.ValueSet(bm.ReactionTarget(reaction_from_string('b_[a]_ppi+_a_[b]'))),
                                                                                                      Union(bm.Intersection(bm.ValueSet(bm.StateTarget(state_from_string('c--a'))), bm.ValueSet(bm.StateTarget(state_from_string('d--e')))),
                                                                                                            bm.Intersection(bm.ValueSet(bm.ReactionTarget(reaction_from_string('e_[a]_ppi-_a_[e]'))), bm.ValueSet(bm.StateTarget(state_from_string('f--r'))))
                                                                                                            )
                                                                                                      )
                                     )
                            ),
        BooleanRuleTestCase('<a-{p} | <a_ppi+_b &c--d>>',
                           Union(ValueSet(bm.StateTarget(state_from_string('a-{p}'))),
                                 Intersection(ValueSet(bm.ReactionTarget(reaction_from_string('a_ppi+_b'))), ValueSet(bm.StateTarget(state_from_string('c--d')))))
                           ),


        BooleanRuleTestCase('<a_syn_b | < < <<a--b & a_[b]--0> & <a--b & b_[a]--0 >> & <! c_deg_a & <a_ppi+_b & <a_[b]--0 & b_[a]--0>>> > | <a--b & <! c_deg_a & ! a_ppi-_b>> > >',
                              Union(ValueSet(bm.ReactionTarget(reaction_from_string('a_syn_b'))),
                                    Union(
                                            Intersection(
                                                        Intersection(
                                                                    Intersection(ValueSet(bm.StateTarget(state_from_string('a--b'))),
                                                                            ValueSet(bm.StateTarget(state_from_string('a_[b]--0')))
                                                                            ),
                                                                    Intersection(ValueSet(bm.StateTarget(state_from_string('a--b'))),
                                                                            ValueSet(bm.StateTarget(state_from_string('b_[a]--0')))
                                                                                )
                                                                    ),
                                                        Intersection(Complement(ValueSet(bm.ReactionTarget(reaction_from_string('c_deg_a')))),
                                                                     Intersection(ValueSet(bm.ReactionTarget(reaction_from_string('a_ppi+_b'))),
                                                                                  Intersection(ValueSet(bm.StateTarget(state_from_string('a_[b]--0'))),
                                                                                               ValueSet(bm.StateTarget(state_from_string('b_[a]--0')))
                                                                                               )
                                                                                  )
                                                                     )
                                                        ),
                                            Intersection(ValueSet(bm.StateTarget(state_from_string('a--b'))),
                                                         Intersection(Complement(ValueSet(bm.ReactionTarget(reaction_from_string('c_deg_a')))),
                                                                      Complement(ValueSet(bm.ReactionTarget(reaction_from_string('a_ppi-_b'))))
                                                                      )
                                                         )
                                        )
                                    )
                              )

    ]
def test_rule_from_str(boolean_test_cases):
    for the_case in boolean_test_cases:
        actual_boolean_factor = bm.rxncon_bool_str_to_venn(the_case.boolean_string)
        actual_boolean_factor.is_equivalent_to(the_case.boolean_factor)


def test_create_eda_compartible_str():
    venn = bm.rxncon_bool_str_to_venn('<a--a | <b_ppi+_a | <<c--a&d--e> | <e_ppi-_a&f--r> >>>')
