import pytest
from collections import namedtuple
from networkx import DiGraph
from rxncon.input.quick.quick import Quick
from rxncon.simulation.rule_graph.reaction_graph import ReactionGraph, EdgeWith, EdgeType, NodeType

RuleTestCase = namedtuple('RuleTestCase', ['quick_string', 'node_tuples', 'edge_tuples'])

def _get_state_nodes(node_tuples, expected_graph):
    """
    Adding nodes to the expected graph.

    Args:
        test_case: Holding information of the current test.
        expected_graph: The graph defined in the test_case.

    Returns:
        expected graph.

    """
    for node_id, node_label, node_type in node_tuples:
        expected_graph.add_node(node_id, type=node_type.value, label=node_label)
    return expected_graph

def _is_graph_test_case_correct(actual_graph, test_case) -> bool:
    """
    Checking if the created and expected graph are equal.

    Args:
        actual_graph: The graph which should be tested.
        test_case: Holding information of the current test e.g. the expected graph..

    Returns:
        bool: True if all the nodes and all the edges as expected, otherwise False.

    """

    expected_graph = DiGraph()

    expected_graph = _get_state_nodes(test_case.node_tuples, expected_graph)

    [expected_graph.add_edge(source, target, interaction=interaction.value, width=width.value)
     for source, target, width, interaction in test_case.edge_tuples]

    for edge in expected_graph.edge:
        assert expected_graph.edge[edge] == actual_graph.edge[edge]
    return expected_graph.node == actual_graph.node and expected_graph.edge == actual_graph.edge

def _create_reaction_graph(quick_string):
    """
    Creating a regulatory graph.

    Args:
        quick_string: A rxncon system in quick format.

    Returns:
        A regulatory graph.

    """
    actual_system = Quick(quick_string)
    reaction_system = ReactionGraph(actual_system.rxncon_system)
    return reaction_system.to_graph()

def test_simple_system():

    test_case = RuleTestCase('''A_[b]_ppi+_B_[a]''',
                             [('A', 'A', NodeType.component), ('A_[b]', 'b', NodeType.domain),
                              ('B', 'B', NodeType.component), ('B_[a]', 'a', NodeType.domain)],
                             [('A', 'A_[b]', EdgeWith.internal, EdgeType.interaction),
                              ('B', 'B_[a]', EdgeWith.internal, EdgeType.interaction),
                              ('A_[b]', 'B_[a]', EdgeWith.external, EdgeType.interaction)])

    assert _is_graph_test_case_correct(_create_reaction_graph(test_case.quick_string), test_case)

