import pytest
import networkx as nex
from collections import namedtuple
import rxncon.input.quick.quick as qui
import rxncon.simulation.rule_graph.regulatory_graph as reg
import rxncon.simulation.rule_graph.graphML as gml

RuleTestCase = namedtuple('RuleTestCase', ['quick_string', 'reaction_node_strings', 'state_node_strings',
                                           'boolean_state_node_strings','edge_tuples'])

def test_graph_generation(the_cases):
    for test_case in the_cases:
        assert is_graph_test_case_correct(test_case)


@pytest.fixture
def the_cases(case_and_expected_graph):
    return case_and_expected_graph


@pytest.fixture
def case_and_expected_graph():
    return [
        RuleTestCase('''A_ppi_B; ! A-{p}
                        C_p+_A''',
                     ['A_ppi_B', 'C_p+_A'],
                     ['A--B', 'A-{p}'],
                     [],
                     [('A_ppi_B', 'A--B', reg.EdgeInteractionType.produce.value),
                      ('C_p+_A', 'A-{p}', reg.EdgeInteractionType.produce.value),
                      ('A-{p}', 'A_ppi_B', reg.EdgeInteractionType.required.value)]),
        RuleTestCase('''A_ppi_B; ! <comp>; ! C-{p}
                        <comp>; AND A-{p}; AND A--C
                        A_ppi_C
                        C_p+_A
                        D_p+_C''',
                     ['A_ppi_B', 'A_ppi_C', 'C_p+_A', 'D_p+_C'],
                     ['A--B', 'A--C', 'A-{p}', 'C-{p}'],
                     ['comp#AND'],
                     [('A_ppi_B', 'A--B', reg.EdgeInteractionType.produce.value),
                      ('A_ppi_C', 'A--C', reg.EdgeInteractionType.produce.value),
                      ('C_p+_A', 'A-{p}', reg.EdgeInteractionType.produce.value),
                      ('D_p+_C', 'C-{p}', reg.EdgeInteractionType.produce.value),
                      ('C-{p}', 'A_ppi_B', reg.EdgeInteractionType.required.value),
                      ('comp', 'A_ppi_B', reg.EdgeInteractionType.required.value),
                      ('A-{p}', 'comp', reg.EdgeInteractionType.AND.value),
                      ('A--C', 'comp', reg.EdgeInteractionType.AND.value)]),
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
                     ['A_ppi_B', 'A_ppi_C', 'A_ppi_D', 'A_ppi_E', 'A_ppi_F', 'A_ppi_G' ],
                     ['A--B', 'A--C', 'A--D', 'A--E', 'A--F', 'A--G'],
                     ['comp#AND', 'comp1#OR', 'comp2#AND', 'comp3#AND'],
                     [('A_ppi_B', 'A--B', reg.EdgeInteractionType.produce.value),
                      ('A_ppi_C', 'A--C', reg.EdgeInteractionType.produce.value),
                      ('A_ppi_D', 'A--D', reg.EdgeInteractionType.produce.value),
                      ('A_ppi_E', 'A--E', reg.EdgeInteractionType.produce.value),
                      ('A_ppi_F', 'A--F', reg.EdgeInteractionType.produce.value),
                      ('A_ppi_G', 'A--G', reg.EdgeInteractionType.produce.value),
                      ('comp', 'A_ppi_B', reg.EdgeInteractionType.required.value),
                      ('comp1', 'comp', reg.EdgeInteractionType.AND.value),
                      ('comp2', 'comp', reg.EdgeInteractionType.AND.value),
                      ('comp3', 'comp1', reg.EdgeInteractionType.OR.value),
                      ('A--C', 'comp1', reg.EdgeInteractionType.OR.value),
                      ('A--D', 'comp2', reg.EdgeInteractionType.AND.value),
                      ('A--E', 'comp2', reg.EdgeInteractionType.AND.value),
                      ('A--F', 'comp3', reg.EdgeInteractionType.AND.value),
                      ('A--G', 'comp3', reg.EdgeInteractionType.AND.value),]
                     )
    ]


def is_graph_test_case_correct(test_case) -> bool:
    actual_system = qui.Quick(test_case.quick_string)
    reg_system = reg.RegulatoryGraph(actual_system.rxncon_system)
    actual_graph = reg_system.to_graph()
    xgmml_system = gml.XGMML(actual_graph, 'test_graph')
    xgmml_system.xgmmlwriter("/home/thiemese/project/rxncon/graphml/test.xgmml")
    expected_graph = nex.DiGraph()
    [expected_graph.add_node(node, type=reg.NodeType.reaction.value) for node in test_case.reaction_node_strings]
    [expected_graph.add_node(node, type=reg.NodeType.state.value) for node in test_case.state_node_strings]
    [expected_graph.add_node(node.split('#AND')[0], type=reg.NodeType.AND.value) for node in test_case.boolean_state_node_strings if '#AND' in node]
    [expected_graph.add_node(node.split('#OR')[0], type=reg.NodeType.OR.value) for node in test_case.boolean_state_node_strings if '#OR' in node]
    [expected_graph.add_edge(source, target, interaction=interaction) for source, target, interaction in test_case.edge_tuples]

    return expected_graph.node == actual_graph.node and expected_graph.edge == actual_graph.edge
