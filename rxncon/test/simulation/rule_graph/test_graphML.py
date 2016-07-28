import pytest
from collections import namedtuple
import os
import time
import tempfile

import rxncon.input.quick.quick as qui
import rxncon.simulation.rule_graph.regulatory_graph as reg
import rxncon.simulation.rule_graph.graphML as gml

RuleTestCase = namedtuple('RuleTestCase', ['quick_string', 'expected_xgmml'])

def test_graph_generation(the_cases):
   for test_case in the_cases:
       assert is_xgmml_export_test_case_correct(test_case)

def test_graph_writing(the_cases):
    for test_case in the_cases:
        assert can_xgmml_be_written_to_file(test_case)


@pytest.fixture
def the_cases(case_and_expected_graph):
    return case_and_expected_graph


@pytest.fixture
def case_and_expected_graph():
    return [
        # RuleTestCase('''A_[b]_ppi_B_[a]; ! A-{p}
        #                 C_p+_A_[(x)]''',
        #              '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
        #                 <graph directed="1"  xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns="http://www.cs.rpi.edu/XGMML">
        #                 <att name="selected" value="1" type="boolean" />
        #                 <att name="name" value="test_graph" type="string"/>
        #                 <att name="shared name" value="test_graph" type="string"/>
        #
        #                 <node id="A--B" label="A--B"><att name="type" value="state" /></node>
        #                 <node id="C_p+_A" label="C_p+_A"><att name="type" value="reaction" /></node>
        #                 <node id="A_ppi_B" label="A_ppi_B"><att name="type" value="reaction" /></node>
        #                 <node id="A-{p}" label="A-{p}"><att name="type" value="state" /></node>
        #                 <edge source="C_p+_A" target="A-{p}"><att name="interaction" value="produce"/></edge>
        #                 <edge source="A_ppi_B" target="A--B"><att name="interaction" value="produce"/></edge>
        #                 <edge source="A-{p}" target="A_ppi_B"><att name="interaction" value="!"/></edge>
        #                 </graph>'''),

        RuleTestCase('''A_[b]_ppi_B_[a]; ! <comp>; ! C-{p}
                        <comp>; AND A-{p}; AND A--C; AND A--D
                        A_[d]_ppi_D_[a]
                        A_[c]_ppi_C_[a]
                        C_p+_A_[(x)]
                        D_p+_C_[(x)]''',
                     '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
                        <graph directed="1"  xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns="http://www.cs.rpi.edu/XGMML">
                        <att name="selected" value="1" type="boolean" />
                        <att name="name" value="test_graph" type="string"/>
                        <att name="shared name" value="test_graph" type="string"/>

                        <node id="A_[c]--C_[a]" label="A_[c]--C_[a]"><att name="type" value="state" /></node>
                        <node id="A_[c]_ppi_C_[a]" label="A_[c]_ppi_C_[a]"><att name="type" value="reaction" /></node>
                        <node id="C_[(x)]-{p}" label="C_[(x)]-{p}"><att name="type" value="state" /></node>
                        <node id="D_p+_C_[(x)]" label="D_p+_C_[(x)]"><att name="type" value="reaction" /></node>
                        <node id="A_[b]_ppi_B_[a]" label="A_[b]_ppi_B_[a]"><att name="type" value="reaction" /></node>
                        <node id="comp" label="comp"><att name="type" value="boolean_and" /></node>
                        <node id="C_p+_A_[(x)]" label="C_p+_A_[(x)]"><att name="type" value="reaction" /></node>
                        <node id="A_[(x)]-{p}" label="A_[(x)]-{p}"><att name="type" value="state" /></node>
                        <node id="A_[b]--B_[a]" label="A_[b]--B_[a]"><att name="type" value="state" /></node>
                        <node id="A_[d]--D_[a]" label="A_[d]--D_[a]"><att name="type" value="state" /></node>
                        <node id="A_[d]_ppi_D_[a]" label="A_[d]_ppi_D_[a]"><att name="type" value="reaction" /></node>
                        <edge source="A_[c]--C_[a]" target="comp"><att name="interaction" value="AND"/></edge>
                        <edge source="A_[c]_ppi_C_[a]" target="A_[c]--C_[a]"><att name="interaction" value="produce"/></edge>
                        <edge source="C_[(x)]-{p}" target="A_[b]_ppi_B_[a]"><att name="interaction" value="!"/></edge>
                        <edge source="D_p+_C_[(x)]" target="C_[(x)]-{p}"><att name="interaction" value="produce"/></edge>
                        <edge source="A_[b]_ppi_B_[a]" target="A_[b]--B_[a]"><att name="interaction" value="produce"/></edge>
                        <edge source="comp" target="A_[b]_ppi_B_[a]"><att name="interaction" value="!"/></edge>
                        <edge source="C_p+_A_[(x)]" target="A_[(x)]-{p}"><att name="interaction" value="produce"/></edge>
                        <edge source="A_[(x)]-{p}" target="comp"><att name="interaction" value="AND"/></edge>
                        <edge source="A_[d]--D_[a]" target="comp"><att name="interaction" value="AND"/></edge>
                        <edge source="A_[d]_ppi_D_[a]" target="A_[d]--D_[a]"><att name="interaction" value="produce"/></edge>
                        </graph> '''),
        # RuleTestCase('''A_ppi_B; ! <comp>
        #                 <comp>; AND <comp1>; AND <comp2>
        #                 <comp1>; OR <comp3>; OR A--C
        #                 <comp2>; AND A--D; AND A--E
        #                 <comp3>; AND A--F; AND A--G
        #                 A_ppi_C
        #                 A_ppi_D
        #                 A_ppi_E
        #                 A_ppi_F
        #                 A_ppi_G''',
        #              '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
        #                     <graph directed="1"  xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns="http://www.cs.rpi.edu/XGMML">
        #                     <att name="selected" value="1" type="boolean" />
        #                     <att name="name" value="test_graph" type="string"/>
        #                     <att name="shared name" value="test_graph" type="string"/>
        #
        #                 <node id="A--D" label="A--D"><att name="type" value="state" /></node>
        #                 <node id="A--E" label="A--E"><att name="type" value="state" /></node>
        #                 <node id="A_ppi_D" label="A_ppi_D"><att name="type" value="reaction" /></node>
        #                 <node id="A--C" label="A--C"><att name="type" value="state" /></node>
        #                 <node id="A--G" label="A--G"><att name="type" value="state" /></node>
        #                 <node id="comp" label="comp"><att name="type" value="boolean_and" /></node>
        #                 <node id="comp3" label="comp3"><att name="type" value="boolean_and" /></node>
        #                 <node id="comp1" label="comp1"><att name="type" value="boolean_or" /></node>
        #                 <node id="A_ppi_C" label="A_ppi_C"><att name="type" value="reaction" /></node>
        #                 <node id="comp2" label="comp2"><att name="type" value="boolean_and" /></node>
        #                 <node id="A_ppi_G" label="A_ppi_G"><att name="type" value="reaction" /></node>
        #                 <node id="A_ppi_B" label="A_ppi_B"><att name="type" value="reaction" /></node>
        #                 <node id="A--F" label="A--F"><att name="type" value="state" /></node>
        #                 <node id="A_ppi_E" label="A_ppi_E"><att name="type" value="reaction" /></node>
        #                 <node id="A_ppi_F" label="A_ppi_F"><att name="type" value="reaction" /></node>
        #                 <node id="A--B" label="A--B"><att name="type" value="state" /></node>
        #                 <edge source="A--D" target="comp2"><att name="interaction" value="AND"/></edge>
        #                 <edge source="A--E" target="comp2"><att name="interaction" value="AND"/></edge>
        #                 <edge source="A_ppi_D" target="A--D"><att name="interaction" value="produce"/></edge>
        #                 <edge source="A--C" target="comp1"><att name="interaction" value="OR"/></edge>
        #                 <edge source="A--G" target="comp3"><att name="interaction" value="AND"/></edge>
        #                 <edge source="comp" target="A_ppi_B"><att name="interaction" value="!"/></edge>
        #                 <edge source="comp3" target="comp1"><att name="interaction" value="OR"/></edge>
        #                 <edge source="comp1" target="comp"><att name="interaction" value="AND"/></edge>
        #                 <edge source="A_ppi_C" target="A--C"><att name="interaction" value="produce"/></edge>
        #                 <edge source="comp2" target="comp"><att name="interaction" value="AND"/></edge>
        #                 <edge source="A_ppi_G" target="A--G"><att name="interaction" value="produce"/></edge>
        #                 <edge source="A_ppi_B" target="A--B"><att name="interaction" value="produce"/></edge>
        #                 <edge source="A--F" target="comp3"><att name="interaction" value="AND"/></edge>
        #                 <edge source="A_ppi_E" target="A--E"><att name="interaction" value="produce"/></edge>
        #                 <edge source="A_ppi_F" target="A--F"><att name="interaction" value="produce"/></edge>
        #                 </graph>''')
    ]

def is_xgmml_export_test_case_correct(test_case):
    actual_system = qui.Quick(test_case.quick_string)
    reg_system = reg.RegulatoryGraph(actual_system.rxncon_system)
    actual_graph = reg_system.to_graph()
    gml_system = gml.XGMML(actual_graph, "test_graph")

    expected_xgmml_lines= [line.strip() for line in test_case.expected_xgmml.split("\n")]
    produced_xgmml_lines = [line.strip() for line in gml_system.to_string().split("\n")]

    if not sorted(produced_xgmml_lines) == sorted(expected_xgmml_lines):
        for i, line in enumerate(produced_xgmml_lines):
            if line not in expected_xgmml_lines:
                print("generated line: {}".format(line))
                print("expected line: {}".format(expected_xgmml_lines[i]))

    return sorted(produced_xgmml_lines) == sorted(expected_xgmml_lines)

def can_xgmml_be_written_to_file(test_case):
    actual_system = qui.Quick(test_case.quick_string)
    reg_system = reg.RegulatoryGraph(actual_system.rxncon_system)
    actual_graph = reg_system.to_graph()
    name = "test{0}".format(time.time())
    gml_system = gml.XGMML(actual_graph, "name")

    path = "{0}/{1}.xgmml".format(tempfile.gettempdir(), name)
    gml_system.to_file(path)
    file_exists = os.path.exists(path)
    os.remove(path)
    return file_exists

def test_insulin_reg_graph():
    quick_system = qui.Quick("""IR_[IRBD]_ppi_IR_[IRBD]
IR_[lig]_i_insulin_[IR]; ! <IR-empty>
IR_p+_IR_[TK(Y1158)]; ! <IR0-IR1-Insulin2>
IR_p+_IR_[TK(Y1162)]; ! <IR0-IR1-Insulin2>
IR_p+_IR_[TK(Y1163)]; ! <IR0-IR1-Insulin2>
IR_ap_IR_[JM(Y972)]; ! <IR-IR-active>
IR_[JMY972]_ppi_IRS_[PTB]; ! IR_[JM(Y972)]-{P}; ! IRS_[PH]--PM_[phospholipids]
IRS_[PH]_i_PM_[phospholipids]
IR_p+_IRS_[(Y)]; ! IR_[JMY972]--IRS_[PTB]; ! <IR-active>
Grb2_[SOS]_ppi_SOS_[Grb2]
Grb2_[SH2]_ppi_IRS_[Y]; ! IRS_[(Y)]-{P}
IR_[JMY972]_ppi_Shc_[PTB]; ! IR_[JM(Y972)]-{P}
IR_p+_Shc_[(YY239240)]; ! IR_[JMY972]--Shc_[PTB]; ! <IR-active>
IR_p+_Shc_[(Y317)]; ! IR_[JMY972]--Shc_[PTB]; ! <IR-active>
Grb2_[SH2]_ppi_Shc_[YY]; ! Shc_[(YY239240)]-{P}
Grb2_[SH2]_ppi_Shc_[Y]; ! Shc_[(Y317)]-{P}
IRS_[Y]_ppi_PI3K_[SH2]; ! IRS_[(Y)]-{P}

<IR-empty>; AND IR@0--IR@2; AND IR@0_[lig]--0; AND IR@2_[lig]--0

<IR-IR-active>; AND <IR-phos>; AND <IR0-IR1-Insulin2>
<IR-phos>; AND IR_[TK(Y1158)]-{P}; AND IR_[TK(Y1162)]-{P}; AND IR_[TK(Y1163)]-{P}

<IR0-IR1-Insulin2>; OR <IR0-insulin2>; OR <IR1-insulin2>
<IR0-insulin2>; AND IR@0--IR@1; AND IR@0_[lig]--insulin@2
<IR1-insulin2>; AND IR@0--IR@1; AND IR@1_[lig]--insulin@2

<IR-active>; AND <IR-phos>; AND <IR0-IR2-Insulin3>
<IR0-IR2-Insulin3>; OR <IR0-insulin3>; OR <IR2-insulin3>
<IR0-insulin3>; AND IR@0--IR@2; AND IR@0_[lig]--insulin@3
<IR2-insulin3>; AND IR@0--IR@2; AND IR@2_[lig]--insulin@3
""")
    reg_graph = reg.RegulatoryGraph(quick_system.rxncon_system)
    actual_graph = reg_graph.to_graph()
    gml_system = gml.XGMML(actual_graph, "name")
    gml_system.to_file('/home/thiemese/data/reg_graph.xgmml')


"""[PI3KAkt signal]; ! IRS_[Y]--PI3K_[SH2]

[RasErk signalling]; ! <Grb2-SOS>
<Grb2-SOS>; OR <Grb2-Shc-Sos>; OR <Grb2-IRS-Sos>
<Grb2-Shc-Sos>; AND Grb2_[SH2]--Shc_[YY]; AND Grb2--SOS
<Grb2-IRS-Sos>; AND Grb2_[SH2]--IRS_[Y]; AND Grb2--SOS"""