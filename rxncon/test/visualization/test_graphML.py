import os
import tempfile
import time
from xml.dom import minidom

import pytest
from networkx import DiGraph

from rxncon.input.quick.quick import Quick
from rxncon.visualization.graphML import map_layout2xgmml, XGMML, _get_labels_and_coordinates_dict
from rxncon.visualization.regulatory_graph import SpeciesReactionGraph


NODE_LAYOUT_MANIPULATION = os.path.join(os.path.dirname(__file__), 'example_node_layout.xgmml')
EXISTING_FILE = os.path.join(os.path.dirname(__file__), 'test_file_existence.xgmml')
NOT_A_DIRECTORY = os.path.join(os.path.dirname(__file__), 'NOT_A_DIR/test_file_existence.xgmml')


def test_file_exists() -> None:
    graph = DiGraph()
    graph.add_node('A', str_value='A', int_value=1, float_value=1.0)

    gml_system = XGMML(graph, "name")
    with pytest.raises(FileExistsError):
        gml_system.to_file(EXISTING_FILE)

    with pytest.raises(NotADirectoryError):
        gml_system.to_file(NOT_A_DIRECTORY)


def test_write_output() -> None:
    digraph = _digraph_from_quick_str('''A_[b]_ppi+_B_[a]; ! A_[(x)]-{p}
                                      C_p+_A_[(x)]''')

    name = "test{0}".format(time.time())
    gml_system = XGMML(digraph, "name")

    path = "{0}/{1}.xgmml".format(tempfile.gettempdir(), name)
    gml_system.to_file(path)
    file_exists = os.path.exists(path)
    file_size = os.stat(path).st_size
    os.remove(path)
    assert file_exists and file_size > 0


def test_graph_from_simple_system() -> None:
    """
    Testing a simple system of 2 reactions and 1 contingency.

    Returns:
        None

    Raises:
        AsssertionError: If the created xgmml differes from the expected.

    """
    digraph = _digraph_from_quick_str('''A_[b]_ppi+_B_[a]; ! A_[(x)]-{p}
                                      C_p+_A_[(x)]''')

    expected_xgmml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
                        <graph directed="1"  xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns="http://www.cs.rpi.edu/XGMML">
                        <att name="selected" value="1" type="boolean" />
                        <att name="name" value="test_graph" type="string"/>
                        <att name="shared name" value="test_graph" type="string"/>

                        <node id="A_[b]--0" label="A_[b]--0"><att name="rxnconID" value="A_[b]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="B_[a]--0" label="B_[a]--0"><att name="rxnconID" value="B_[a]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[b]--B_[a]" label="A_[b]--B_[a]"><att name="rxnconID" value="A_[b]--B_[a]" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[b]_ppi+_B_[a]" label="A_[b]_ppi+_B_[a]"><att name="rxnconID" value="A_[b]_ppi+_B_[a]" type="string"/><att name="type" value="reaction" type="string"/></node>
                        <node id="A_[(x)]-{p}" label="A_[(x)]-{p}"><att name="rxnconID" value="A_[(x)]-{p}" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[(x)]-{0}" label="A_[(x)]-{0}"><att name="rxnconID" value="A_[(x)]-{0}" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="C_p+_A_[(x)]" label="C_p+_A_[(x)]"><att name="rxnconID" value="C_p+_A_[(x)]" type="string"/><att name="type" value="reaction" type="string"/></node>
                        <edge source="A_[b]--0" target="A_[b]_ppi+_B_[a]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="B_[a]--0" target="A_[b]_ppi+_B_[a]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="A_[b]_ppi+_B_[a]" target="A_[b]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="A_[b]_ppi+_B_[a]" target="B_[a]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="A_[b]_ppi+_B_[a]" target="A_[b]--B_[a]"><att name="interaction" value="produce" type="string"/></edge>
                        <edge source="A_[(x)]-{p}" target="A_[b]_ppi+_B_[a]"><att name="interaction" value="!" type="string"/></edge>
                        <edge source="A_[(x)]-{0}" target="C_p+_A_[(x)]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="C_p+_A_[(x)]" target="A_[(x)]-{p}"><att name="interaction" value="produce" type="string"/></edge>
                        <edge source="C_p+_A_[(x)]" target="A_[(x)]-{0}"><att name="interaction" value="consume" type="string"/></edge>
                        </graph>"""

    assert is_digraph_to_xgmml_correct(digraph, expected_xgmml)


def test_graph_from_simple_boolean() -> None:
    """
    Testing if a DiGraph can be written as xgmml and writen to file.

    Note:
        A system of 5 reactions, 1 boolean contingency, 1 state contingency.

    Returns:
        None

    Raises:
        AsssertionError: If the created xgmml differes from the expected.

    """
    digraph = _digraph_from_quick_str('''A_[b]_ppi+_B_[a]; ! <comp>; ! C_[(x)]-{p}
                                      <comp>; AND A_[(x)]-{p}; AND A_[c]--C_[a]; AND A_[d]--D_[a]
                                      A_[d]_ppi+_D_[a]
                                      A_[c]_ppi+_C_[a]
                                      C_p+_A_[(x)]
                                      D_p+_C_[(x)]''')

    expected_xgmml = '''
                        <?xml version="1.0" encoding="UTF-8" standalone="yes"?>
                        <graph directed="1"  xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns="http://www.cs.rpi.edu/XGMML">
                        <att name="selected" value="1" type="boolean" />
                        <att name="name" value="test_graph" type="string"/>
                        <att name="shared name" value="test_graph" type="string"/>

                        <node id="A_[d]--0" label="A_[d]--0"><att name="rxnconID" value="A_[d]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[b]--B_[a]" label="A_[b]--B_[a]"><att name="rxnconID" value="A_[b]--B_[a]" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="D_p+_C_[(x)]" label="D_p+_C_[(x)]"><att name="rxnconID" value="D_p+_C_[(x)]" type="string"/><att name="type" value="reaction" type="string"/></node>
                        <node id="A_[c]--0" label="A_[c]--0"><att name="rxnconID" value="A_[c]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[b]--0" label="A_[b]--0"><att name="rxnconID" value="A_[b]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[b]_ppi+_B_[a]" label="A_[b]_ppi+_B_[a]"><att name="rxnconID" value="A_[b]_ppi+_B_[a]" type="string"/><att name="type" value="reaction" type="string"/></node>
                        <node id="A_[c]--C_[a]" label="A_[c]--C_[a]"><att name="rxnconID" value="A_[c]--C_[a]" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="comp" label="comp"><att name="rxnconID" value="comp" type="string"/><att name="type" value="boolean_and" type="string"/></node>
                        <node id="A_[d]_ppi+_D_[a]" label="A_[d]_ppi+_D_[a]"><att name="rxnconID" value="A_[d]_ppi+_D_[a]" type="string"/><att name="type" value="reaction" type="string"/></node>
                        <node id="C_[a]--0" label="C_[a]--0"><att name="rxnconID" value="C_[a]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="C_[(x)]-{0}" label="C_[(x)]-{0}"><att name="rxnconID" value="C_[(x)]-{0}" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[c]_ppi+_C_[a]" label="A_[c]_ppi+_C_[a]"><att name="rxnconID" value="A_[c]_ppi+_C_[a]" type="string"/><att name="type" value="reaction" type="string"/></node>
                        <node id="C_p+_A_[(x)]" label="C_p+_A_[(x)]"><att name="rxnconID" value="C_p+_A_[(x)]" type="string"/><att name="type" value="reaction" type="string"/></node>
                        <node id="C_[(x)]-{p}" label="C_[(x)]-{p}"><att name="rxnconID" value="C_[(x)]-{p}" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[(x)]-{p}" label="A_[(x)]-{p}"><att name="rxnconID" value="A_[(x)]-{p}" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[d]--D_[a]" label="A_[d]--D_[a]"><att name="rxnconID" value="A_[d]--D_[a]" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="D_[a]--0" label="D_[a]--0"><att name="rxnconID" value="D_[a]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[(x)]-{0}" label="A_[(x)]-{0}"><att name="rxnconID" value="A_[(x)]-{0}" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="B_[a]--0" label="B_[a]--0"><att name="rxnconID" value="B_[a]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <edge source="A_[d]--0" target="A_[d]_ppi+_D_[a]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="D_p+_C_[(x)]" target="C_[(x)]-{p}"><att name="interaction" value="produce" type="string"/></edge>
                        <edge source="D_p+_C_[(x)]" target="C_[(x)]-{0}"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="A_[c]--0" target="A_[c]_ppi+_C_[a]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="A_[b]--0" target="A_[b]_ppi+_B_[a]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="A_[b]_ppi+_B_[a]" target="A_[b]--B_[a]"><att name="interaction" value="produce" type="string"/></edge>
                        <edge source="A_[b]_ppi+_B_[a]" target="B_[a]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="A_[b]_ppi+_B_[a]" target="A_[b]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="A_[c]--C_[a]" target="comp"><att name="interaction" value="AND" type="string"/></edge>
                        <edge source="comp" target="A_[b]_ppi+_B_[a]"><att name="interaction" value="!" type="string"/></edge>
                        <edge source="A_[d]_ppi+_D_[a]" target="A_[d]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="A_[d]_ppi+_D_[a]" target="A_[d]--D_[a]"><att name="interaction" value="produce" type="string"/></edge>
                        <edge source="A_[d]_ppi+_D_[a]" target="D_[a]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="C_[a]--0" target="A_[c]_ppi+_C_[a]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="C_[(x)]-{0}" target="D_p+_C_[(x)]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="A_[c]_ppi+_C_[a]" target="A_[c]--C_[a]"><att name="interaction" value="produce" type="string"/></edge>
                        <edge source="A_[c]_ppi+_C_[a]" target="A_[c]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="A_[c]_ppi+_C_[a]" target="C_[a]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="C_p+_A_[(x)]" target="A_[(x)]-{0}"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="C_p+_A_[(x)]" target="A_[(x)]-{p}"><att name="interaction" value="produce" type="string"/></edge>
                        <edge source="C_[(x)]-{p}" target="A_[b]_ppi+_B_[a]"><att name="interaction" value="!" type="string"/></edge>
                        <edge source="A_[(x)]-{p}" target="comp"><att name="interaction" value="AND" type="string"/></edge>
                        <edge source="A_[d]--D_[a]" target="comp"><att name="interaction" value="AND" type="string"/></edge>
                        <edge source="D_[a]--0" target="A_[d]_ppi+_D_[a]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="A_[(x)]-{0}" target="C_p+_A_[(x)]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="B_[a]--0" target="A_[b]_ppi+_B_[a]"><att name="interaction" value="ss" type="string"/></edge>
                        </graph>'''

    assert is_digraph_to_xgmml_correct(digraph, expected_xgmml)


def test_graph_from_complex_boolean() -> None:
    """
    Testing a nested boolean.

     Note:
         A system of 6 reactions, 1 boolean contingency.

    Returns:
        None

    Raises:
        AsssertionError: If the created xgmml differes from the expected.

    """
    digraph = _digraph_from_quick_str("""A_ppi+_B; ! <comp>
                                      <comp>; AND <comp1>; AND <comp2>
                                      <comp1>; OR <comp3>; OR A--C
                                      <comp2>; AND A--D; AND A--E
                                      <comp3>; AND A--F; AND A--G
                                      A_ppi+_C
                                      A_ppi+_D
                                      A_ppi+_E
                                      A_ppi+_F
                                      A_ppi+_G""")

    expected_xgmml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
                        <graph directed="1"  xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns="http://www.cs.rpi.edu/XGMML">
                        <att name="selected" value="1" type="boolean" />
                        <att name="name" value="test_graph" type="string"/>
                        <att name="shared name" value="test_graph" type="string"/>

                        <node id="A_[C]--C_[A]" label="A_[C]--C_[A]"><att name="rxnconID" value="A_[C]--C_[A]" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="C_[A]--0" label="C_[A]--0"><att name="rxnconID" value="C_[A]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[D]--0" label="A_[D]--0"><att name="rxnconID" value="A_[D]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[F]--F_[A]" label="A_[F]--F_[A]"><att name="rxnconID" value="A_[F]--F_[A]" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[E]--E_[A]" label="A_[E]--E_[A]"><att name="rxnconID" value="A_[E]--E_[A]" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="G_[A]--0" label="G_[A]--0"><att name="rxnconID" value="G_[A]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[B]--B_[A]" label="A_[B]--B_[A]"><att name="rxnconID" value="A_[B]--B_[A]" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="comp2" label="comp2"><att name="rxnconID" value="comp2" type="string"/><att name="type" value="boolean_and" type="string"/></node>
                        <node id="A_[G]--G_[A]" label="A_[G]--G_[A]"><att name="rxnconID" value="A_[G]--G_[A]" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[E]--0" label="A_[E]--0"><att name="rxnconID" value="A_[E]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[B]--0" label="A_[B]--0"><att name="rxnconID" value="A_[B]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="comp" label="comp"><att name="rxnconID" value="comp" type="string"/><att name="type" value="boolean_and" type="string"/></node>
                        <node id="A_[B]_ppi+_B_[A]" label="A_[B]_ppi+_B_[A]"><att name="rxnconID" value="A_[B]_ppi+_B_[A]" type="string"/><att name="type" value="reaction" type="string"/></node>
                        <node id="B_[A]--0" label="B_[A]--0"><att name="rxnconID" value="B_[A]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[D]--D_[A]" label="A_[D]--D_[A]"><att name="rxnconID" value="A_[D]--D_[A]" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[C]--0" label="A_[C]--0"><att name="rxnconID" value="A_[C]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="comp1" label="comp1"><att name="rxnconID" value="comp1" type="string"/><att name="type" value="boolean_or" type="string"/></node>
                        <node id="F_[A]--0" label="F_[A]--0"><att name="rxnconID" value="F_[A]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[E]_ppi+_E_[A]" label="A_[E]_ppi+_E_[A]"><att name="rxnconID" value="A_[E]_ppi+_E_[A]" type="string"/><att name="type" value="reaction" type="string"/></node>
                        <node id="A_[G]--0" label="A_[G]--0"><att name="rxnconID" value="A_[G]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[C]_ppi+_C_[A]" label="A_[C]_ppi+_C_[A]"><att name="rxnconID" value="A_[C]_ppi+_C_[A]" type="string"/><att name="type" value="reaction" type="string"/></node>
                        <node id="A_[D]_ppi+_D_[A]" label="A_[D]_ppi+_D_[A]"><att name="rxnconID" value="A_[D]_ppi+_D_[A]" type="string"/><att name="type" value="reaction" type="string"/></node>
                        <node id="comp3" label="comp3"><att name="rxnconID" value="comp3" type="string"/><att name="type" value="boolean_and" type="string"/></node>
                        <node id="A_[F]_ppi+_F_[A]" label="A_[F]_ppi+_F_[A]"><att name="rxnconID" value="A_[F]_ppi+_F_[A]" type="string"/><att name="type" value="reaction" type="string"/></node>
                        <node id="D_[A]--0" label="D_[A]--0"><att name="rxnconID" value="D_[A]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[F]--0" label="A_[F]--0"><att name="rxnconID" value="A_[F]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="E_[A]--0" label="E_[A]--0"><att name="rxnconID" value="E_[A]--0" type="string"/><att name="type" value="state" type="string"/></node>
                        <node id="A_[G]_ppi+_G_[A]" label="A_[G]_ppi+_G_[A]"><att name="rxnconID" value="A_[G]_ppi+_G_[A]" type="string"/><att name="type" value="reaction" type="string"/></node>
                        <edge source="A_[C]--C_[A]" target="comp1"><att name="interaction" value="OR" type="string"/></edge>
                        <edge source="C_[A]--0" target="A_[C]_ppi+_C_[A]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="A_[D]--0" target="A_[D]_ppi+_D_[A]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="A_[F]--F_[A]" target="comp3"><att name="interaction" value="AND" type="string"/></edge>
                        <edge source="A_[E]--E_[A]" target="comp2"><att name="interaction" value="AND" type="string"/></edge>
                        <edge source="G_[A]--0" target="A_[G]_ppi+_G_[A]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="comp2" target="comp"><att name="interaction" value="AND" type="string"/></edge>
                        <edge source="A_[G]--G_[A]" target="comp3"><att name="interaction" value="AND" type="string"/></edge>
                        <edge source="A_[E]--0" target="A_[E]_ppi+_E_[A]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="A_[B]--0" target="A_[B]_ppi+_B_[A]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="comp" target="A_[B]_ppi+_B_[A]"><att name="interaction" value="!" type="string"/></edge>
                        <edge source="A_[B]_ppi+_B_[A]" target="B_[A]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="A_[B]_ppi+_B_[A]" target="A_[B]--B_[A]"><att name="interaction" value="produce" type="string"/></edge>
                        <edge source="A_[B]_ppi+_B_[A]" target="A_[B]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="B_[A]--0" target="A_[B]_ppi+_B_[A]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="A_[D]--D_[A]" target="comp2"><att name="interaction" value="AND" type="string"/></edge>
                        <edge source="A_[C]--0" target="A_[C]_ppi+_C_[A]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="comp1" target="comp"><att name="interaction" value="AND" type="string"/></edge>
                        <edge source="F_[A]--0" target="A_[F]_ppi+_F_[A]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="A_[E]_ppi+_E_[A]" target="A_[E]--E_[A]"><att name="interaction" value="produce" type="string"/></edge>
                        <edge source="A_[E]_ppi+_E_[A]" target="A_[E]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="A_[E]_ppi+_E_[A]" target="E_[A]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="A_[G]--0" target="A_[G]_ppi+_G_[A]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="A_[C]_ppi+_C_[A]" target="A_[C]--C_[A]"><att name="interaction" value="produce" type="string"/></edge>
                        <edge source="A_[C]_ppi+_C_[A]" target="C_[A]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="A_[C]_ppi+_C_[A]" target="A_[C]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="A_[D]_ppi+_D_[A]" target="A_[D]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="A_[D]_ppi+_D_[A]" target="D_[A]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="A_[D]_ppi+_D_[A]" target="A_[D]--D_[A]"><att name="interaction" value="produce" type="string"/></edge>
                        <edge source="comp3" target="comp1"><att name="interaction" value="OR" type="string"/></edge>
                        <edge source="A_[F]_ppi+_F_[A]" target="A_[F]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="A_[F]_ppi+_F_[A]" target="F_[A]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="A_[F]_ppi+_F_[A]" target="A_[F]--F_[A]"><att name="interaction" value="produce" type="string"/></edge>
                        <edge source="D_[A]--0" target="A_[D]_ppi+_D_[A]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="A_[F]--0" target="A_[F]_ppi+_F_[A]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="E_[A]--0" target="A_[E]_ppi+_E_[A]"><att name="interaction" value="ss" type="string"/></edge>
                        <edge source="A_[G]_ppi+_G_[A]" target="A_[G]--G_[A]"><att name="interaction" value="produce" type="string"/></edge>
                        <edge source="A_[G]_ppi+_G_[A]" target="G_[A]--0"><att name="interaction" value="consume" type="string"/></edge>
                        <edge source="A_[G]_ppi+_G_[A]" target="A_[G]--0"><att name="interaction" value="consume" type="string"/></edge>
                        </graph>"""

    assert is_digraph_to_xgmml_correct(digraph, expected_xgmml)


def test_format_attributes() -> None:
    """
    Testing attribute annotation in xgmml.

    Returns:
        None

    Raises:
        AssertionError if an expected attribute type is not within the node attributes.

    """
    graph = DiGraph()
    graph.add_node('A', str_value='A', int_value=1, float_value=1.0)
    format_attribute_system_str = XGMML(graph, "graph").to_string()

    expected_attribute_type = ['integer', 'double', 'string']
    xmldoc_actual = minidom.parseString(format_attribute_system_str)
    actual_node_attribute_types = [_get_node_att_type(xmldoc_actual.getElementsByTagName('node')[0], attribute_name=attribute)
                                   for attribute in ['int_value', 'float_value', 'str_value']]
    while expected_attribute_type:
        attribute_type = expected_attribute_type.pop()
        assert attribute_type in actual_node_attribute_types


def test_layout_remains_correct_when_adding_node() -> None:
    """
    Testing if the new graph gets the correct layout of the old graph if we add an additional node.

    Returns:
        None

    Raises:
        AssertionError: If the node coordinates of the new graph differ from expected coordinates.

    """
    test_case = Quick("""A_[b]_ppi+_B_[a]; ! A_[(c)]-{p}
                         C_p+_A_[(c)]; ! C_[d]--D_[c]
                         C_[d]_ppi+_D_[c]""")

    reg_graph = SpeciesReactionGraph(test_case.rxncon_system)
    gml_system = XGMML(reg_graph.to_graph(), "add_node")
    mapped_layout = map_layout2xgmml(gml_system.to_string(), NODE_LAYOUT_MANIPULATION)
    xmldoc_no_layout = minidom.parseString(mapped_layout)
    xmldoc_no_layout_info = _get_labels_and_coordinates_dict(xmldoc_no_layout)
    xmldoc_expected_layout = minidom.parse(NODE_LAYOUT_MANIPULATION)
    xmldoc_expected_layout_info = _get_labels_and_coordinates_dict(xmldoc_expected_layout)

    assert all(xmldoc_no_layout_info[no_layout_node] == xmldoc_expected_layout_info[no_layout_node]
               for no_layout_node in xmldoc_no_layout_info)
    assert all(xmldoc_no_layout_info[expected_layout_node] == xmldoc_expected_layout_info[expected_layout_node]
               for expected_layout_node in xmldoc_expected_layout_info)


def test_layout_remains_correct_when_removing_node() -> None:
    """
    Testing if the new graph gets tie correct layout of the old graph if we remove a node.

    Returns:
        None

    Raises:
        AssertionError: If the node coordinates of the new graph differ from expected coordinates.

    """
    test_case = Quick("""A_[b]_ppi+_B_[a]""")

    reg_graph = SpeciesReactionGraph(test_case.rxncon_system)
    gml_system = XGMML(reg_graph.to_graph(), "remove_node")
    mapped_layout = map_layout2xgmml(gml_system.to_string(), NODE_LAYOUT_MANIPULATION)
    xmldoc_no_layout = minidom.parseString(mapped_layout)
    xmldoc_no_layout_info = _get_labels_and_coordinates_dict(xmldoc_no_layout)
    xmldoc_expected_layout = minidom.parse(NODE_LAYOUT_MANIPULATION)
    xmldoc_expected_layout_info = _get_labels_and_coordinates_dict(xmldoc_expected_layout)

    assert all(xmldoc_no_layout_info[no_layout_node] == xmldoc_expected_layout_info[no_layout_node]
               for no_layout_node in xmldoc_no_layout_info)


def is_digraph_to_xgmml_correct(digraph: DiGraph, expected_xgmml: str) -> bool:
    """

    Args:
        test_case_quick:
        expected_xgmml: expected xgmml content.

    Returns:
        bool: True if the created xgmml is the same as the expected, False otherwise.

    """
    gml_system = XGMML(digraph, "test_graph")

    expected_xgmml_lines = [line.strip() for line in expected_xgmml.split("\n") if line.strip()]
    produced_xgmml_lines = [line.strip() for line in gml_system.to_string().split("\n") if line.strip()]

    if not sorted(produced_xgmml_lines) == sorted(expected_xgmml_lines):
        for i, line in enumerate(produced_xgmml_lines):
            if line not in expected_xgmml_lines:
                print("generated line: {}".format(line))
                print("expected line:  {}".format(expected_xgmml_lines[i]))

    return sorted(produced_xgmml_lines) == sorted(expected_xgmml_lines)


def _digraph_from_quick_str(quick_str: str) -> DiGraph:
    """
    Creating a regulatory graph.

    Args:
        quick_str: A rxncon system in quick format.

    Returns:
        A regulatory graph.

    """
    actual_system = Quick(quick_str)
    reg_system = SpeciesReactionGraph(actual_system.rxncon_system)
    return reg_system.to_graph()


def _get_node_att_type(element: minidom.Element, attribute_name: str) -> str:
    for child in element.childNodes:
        if child.attributes and child.getAttribute('name') == attribute_name:
            return child.getAttribute('type')
    raise AssertionError
