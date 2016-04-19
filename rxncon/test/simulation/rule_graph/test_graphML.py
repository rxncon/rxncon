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
        RuleTestCase('''A_ppi_B; ! A-{p}
                        C_p+_A''',
                     '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
                        <graph directed="1"  xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns="http://www.cs.rpi.edu/XGMML">
                        <att name="selected" value="1" type="boolean" />
                        <att name="name" value="test_graph" type="string"/>
                        <att name="shared name" value="test_graph" type="string"/>

                        <node id="A--B" label="A--B"><att name="type" value="state" /></node>
                        <node id="C_p+_A" label="C_p+_A"><att name="type" value="reaction" /></node>
                        <node id="A_ppi_B" label="A_ppi_B"><att name="type" value="reaction" /></node>
                        <node id="A-{p}" label="A-{p}"><att name="type" value="state" /></node>
                        <edge source="C_p+_A" target="A-{p}"><att name="interaction" value="produce"/></edge>
                        <edge source="A_ppi_B" target="A--B"><att name="interaction" value="produce"/></edge>
                        <edge source="A-{p}" target="A_ppi_B"><att name="interaction" value="!"/></edge>
                        </graph>'''),

        RuleTestCase('''A_ppi_B; ! <comp>; ! C-{p}
                        <comp>; AND A-{p}; AND A--C
                        A_ppi_C
                        C_p+_A
                        D_p+_C''',
                     '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
                        <graph directed="1"  xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns="http://www.cs.rpi.edu/XGMML">
                        <att name="selected" value="1" type="boolean" />
                        <att name="name" value="test_graph" type="string"/>
                        <att name="shared name" value="test_graph" type="string"/>

                        <node id="A--C" label="A--C"><att name="type" value="state" /></node>
                        <node id="A_ppi_C" label="A_ppi_C"><att name="type" value="reaction" /></node>
                        <node id="C-{p}" label="C-{p}"><att name="type" value="state" /></node>
                        <node id="D_p+_C" label="D_p+_C"><att name="type" value="reaction" /></node>
                        <node id="A_ppi_B" label="A_ppi_B"><att name="type" value="reaction" /></node>
                        <node id="comp" label="comp"><att name="type" value="boolean_and" /></node>
                        <node id="C_p+_A" label="C_p+_A"><att name="type" value="reaction" /></node>
                        <node id="A-{p}" label="A-{p}"><att name="type" value="state" /></node>
                        <node id="A--B" label="A--B"><att name="type" value="state" /></node>
                        <edge source="A--C" target="comp"><att name="interaction" value="AND"/></edge>
                        <edge source="A_ppi_C" target="A--C"><att name="interaction" value="produce"/></edge>
                        <edge source="C-{p}" target="A_ppi_B"><att name="interaction" value="!"/></edge>
                        <edge source="D_p+_C" target="C-{p}"><att name="interaction" value="produce"/></edge>
                        <edge source="A_ppi_B" target="A--B"><att name="interaction" value="produce"/></edge>
                        <edge source="comp" target="A_ppi_B"><att name="interaction" value="!"/></edge>
                        <edge source="C_p+_A" target="A-{p}"><att name="interaction" value="produce"/></edge>
                        <edge source="A-{p}" target="comp"><att name="interaction" value="AND"/></edge>
                        </graph> '''),
        RuleTestCase('''A_ppi_B; ! <comp>
                        <comp>; AND <comp1>; AND <comp2>
                        <comp1>; OR <comp3>; OR A--C
                        <comp2>; AND A--D; AND A--E
                        <comp3>; AND A--F; AND A--G
                        A_ppi_C
                        A_ppi_D
                        A_ppi_E
                        A_ppi_F
                        A_ppi_G''',
                     '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
                            <graph directed="1"  xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns="http://www.cs.rpi.edu/XGMML">
                            <att name="selected" value="1" type="boolean" />
                            <att name="name" value="test_graph" type="string"/>
                            <att name="shared name" value="test_graph" type="string"/>

                        <node id="A--D" label="A--D"><att name="type" value="state" /></node>
                        <node id="A--E" label="A--E"><att name="type" value="state" /></node>
                        <node id="A_ppi_D" label="A_ppi_D"><att name="type" value="reaction" /></node>
                        <node id="A--C" label="A--C"><att name="type" value="state" /></node>
                        <node id="A--G" label="A--G"><att name="type" value="state" /></node>
                        <node id="comp" label="comp"><att name="type" value="boolean_and" /></node>
                        <node id="comp3" label="comp3"><att name="type" value="boolean_and" /></node>
                        <node id="comp1" label="comp1"><att name="type" value="boolean_or" /></node>
                        <node id="A_ppi_C" label="A_ppi_C"><att name="type" value="reaction" /></node>
                        <node id="comp2" label="comp2"><att name="type" value="boolean_and" /></node>
                        <node id="A_ppi_G" label="A_ppi_G"><att name="type" value="reaction" /></node>
                        <node id="A_ppi_B" label="A_ppi_B"><att name="type" value="reaction" /></node>
                        <node id="A--F" label="A--F"><att name="type" value="state" /></node>
                        <node id="A_ppi_E" label="A_ppi_E"><att name="type" value="reaction" /></node>
                        <node id="A_ppi_F" label="A_ppi_F"><att name="type" value="reaction" /></node>
                        <node id="A--B" label="A--B"><att name="type" value="state" /></node>
                        <edge source="A--D" target="comp2"><att name="interaction" value="AND"/></edge>
                        <edge source="A--E" target="comp2"><att name="interaction" value="AND"/></edge>
                        <edge source="A_ppi_D" target="A--D"><att name="interaction" value="produce"/></edge>
                        <edge source="A--C" target="comp1"><att name="interaction" value="OR"/></edge>
                        <edge source="A--G" target="comp3"><att name="interaction" value="AND"/></edge>
                        <edge source="comp" target="A_ppi_B"><att name="interaction" value="!"/></edge>
                        <edge source="comp3" target="comp1"><att name="interaction" value="OR"/></edge>
                        <edge source="comp1" target="comp"><att name="interaction" value="AND"/></edge>
                        <edge source="A_ppi_C" target="A--C"><att name="interaction" value="produce"/></edge>
                        <edge source="comp2" target="comp"><att name="interaction" value="AND"/></edge>
                        <edge source="A_ppi_G" target="A--G"><att name="interaction" value="produce"/></edge>
                        <edge source="A_ppi_B" target="A--B"><att name="interaction" value="produce"/></edge>
                        <edge source="A--F" target="comp3"><att name="interaction" value="AND"/></edge>
                        <edge source="A_ppi_E" target="A--E"><att name="interaction" value="produce"/></edge>
                        <edge source="A_ppi_F" target="A--F"><att name="interaction" value="produce"/></edge>
                        </graph>''')
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

