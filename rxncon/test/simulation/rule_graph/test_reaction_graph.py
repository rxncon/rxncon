import pytest
from typing import Tuple, List
from collections import namedtuple
from networkx import DiGraph
from rxncon.input.quick.quick import Quick
from rxncon.simulation.rule_graph.reaction_graph import ReactionGraph, EdgeWith, EdgeType, NodeType, rxngraph_from_rxncon_system

RuleTestCase = namedtuple('RuleTestCase', ['quick_string', 'node_tuples', 'edge_tuples'])


def _get_state_nodes(node_tuples: List[Tuple[str, str, NodeType]], expected_graph: DiGraph) -> DiGraph:
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


def _is_graph_test_case_correct(actual_graph: DiGraph, test_case: RuleTestCase) -> bool:
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


def _create_reaction_graph(quick_string: str) -> DiGraph:
    """
    Creating a regulatory graph.

    Args:
        quick_string: A rxncon system in quick format.

    Returns:
        A regulatory graph.

    """
    actual_system = Quick(quick_string)
    return rxngraph_from_rxncon_system(actual_system.rxncon_system).reaction_graph


def test_ppi() -> None:
    """
    Testing protein-protein interaction ppi.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RuleTestCase('''A_[b]_ppi+_B_[a]''',
                             [('A', 'A', NodeType.component), ('A_[b]', 'b', NodeType.domain),
                              ('B', 'B', NodeType.component), ('B_[a]', 'a', NodeType.domain)],
                             [('A', 'A_[b]', EdgeWith.internal, EdgeType.interaction),
                              ('B', 'B_[a]', EdgeWith.internal, EdgeType.interaction),
                              ('A_[b]', 'B_[a]', EdgeWith.external, EdgeType.interaction)])

    assert _is_graph_test_case_correct(_create_reaction_graph(test_case.quick_string), test_case)


def test_ipi() -> None:
    """
    Testing intra-protein interaction (ipi.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RuleTestCase('''A_[a1]_ipi+_A_[a2]''',
                             [('A', 'A', NodeType.component), ('A_[a1]', 'a1', NodeType.domain),
                              ('A_[a2]', 'a2', NodeType.domain)],
                             [('A', 'A_[a1]', EdgeWith.internal, EdgeType.interaction),
                              ('A', 'A_[a2]', EdgeWith.internal, EdgeType.interaction),
                              ('A_[a1]', 'A_[a2]', EdgeWith.external, EdgeType.interaction)])

    assert _is_graph_test_case_correct(_create_reaction_graph(test_case.quick_string), test_case)


def test_trans_modification() -> None:
    """
    Testing trans-modification.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RuleTestCase('''A_p+_B_[(a)]''',
                             [('A', 'A', NodeType.component), ('B_[(a)]', 'a', NodeType.residue),
                              ('B', 'B', NodeType.component)],
                             [('B', 'B_[(a)]', EdgeWith.internal, EdgeType.interaction),
                              ('A', 'B_[(a)]', EdgeWith.external, EdgeType.modification)])

    assert _is_graph_test_case_correct(_create_reaction_graph(test_case.quick_string), test_case)


def test_cis_modification() -> None:
    """
    Testing cis-modification.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RuleTestCase('''A_ap+_A_[(a)]''',
                             [('A', 'A', NodeType.component), ('A_[(a)]', 'a', NodeType.residue)],
                             [('A', 'A_[(a)]', EdgeWith.internal, EdgeType.interaction),
                              ('A', 'A_[(a)]', EdgeWith.external, EdgeType.modification)])

    assert _is_graph_test_case_correct(_create_reaction_graph(test_case.quick_string), test_case)


def test_bidirectional_modification() -> None:
    """
    Testing bi-directional modification.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RuleTestCase('''A_[(a)]_pt_B_[(b)]''',
                             [('A', 'A', NodeType.component), ('A_[(a)]', 'a', NodeType.residue),
                              ('B', 'B', NodeType.component), ('B_[(b)]', 'b', NodeType.residue)],
                             [('A', 'A_[(a)]', EdgeWith.internal, EdgeType.interaction),
                              ('B', 'B_[(b)]', EdgeWith.internal, EdgeType.interaction),
                              ('A_[(a)]', 'B_[(b)]', EdgeWith.external, EdgeType.bimodification)])

    assert _is_graph_test_case_correct(_create_reaction_graph(test_case.quick_string), test_case)


def test_degradation() -> None:
    """
    Testing degradation reaction.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RuleTestCase('''A_deg_B''',
                             [('A', 'A', NodeType.component), ('B', 'B', NodeType.component)],
                             [('A', 'B', EdgeWith.external, EdgeType.degradation)])

    assert _is_graph_test_case_correct(_create_reaction_graph(test_case.quick_string), test_case)


def test_synthesis() -> None:
    """
    Testing general synthesis reaction.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RuleTestCase('''A_syn_B''',
                             [('A', 'A', NodeType.component), ('B', 'B', NodeType.component)],
                             [('A', 'B', EdgeWith.external, EdgeType.synthesis)])

    assert _is_graph_test_case_correct(_create_reaction_graph(test_case.quick_string), test_case)


def test_trsc() -> None:
    """
    Testing general synthesis reaction.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RuleTestCase('''A_trsc_B''',
                             [('A', 'A', NodeType.component), ('BGene', 'BGene', NodeType.component)],
                             [('A', 'BGene', EdgeWith.external, EdgeType.synthesis)])

    assert _is_graph_test_case_correct(_create_reaction_graph(test_case.quick_string), test_case)


def test_trsl() -> None:
    """
    Testing general synthesis reaction.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RuleTestCase('''A_trsl_B''',
                             [('A', 'A', NodeType.component), ('BmRNA', 'BmRNA', NodeType.component)],
                             [('A', 'BmRNA', EdgeWith.external, EdgeType.synthesis)])

    assert _is_graph_test_case_correct(_create_reaction_graph(test_case.quick_string), test_case)


def test_multiple_reactions() -> None:
    """
    Testing multiple reactions in combination.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RuleTestCase('''A_[b]_ppi+_B_[a]
                                A_[a1]_ipi+_A_[a2]
                                A_p+_B_[(a)]
                                A_ap+_A_[(a)]
                                A_[(a)]_pt_B_[(b)]''',
                             [('A', 'A', NodeType.component), ('A_[b]', 'b', NodeType.domain),
                              ('A_[a1]', 'a1', NodeType.domain), ('A_[a2]', 'a2', NodeType.domain),
                              ('A_[(a)]', 'a', NodeType.residue),
                              ('B', 'B', NodeType.component), ('B_[a]', 'a', NodeType.domain),
                              ('B_[(a)]', 'a', NodeType.residue), ('B_[(b)]', 'b', NodeType.residue)],
                             [('A', 'A_[b]', EdgeWith.internal, EdgeType.interaction),
                              ('A', 'A_[(a)]', EdgeWith.internal, EdgeType.interaction),
                              ('A', 'A_[a1]', EdgeWith.internal, EdgeType.interaction),
                              ('A', 'A_[a2]', EdgeWith.internal, EdgeType.interaction),
                              ('B', 'B_[(b)]', EdgeWith.internal, EdgeType.interaction),
                              ('B', 'B_[a]', EdgeWith.internal, EdgeType.interaction),
                              ('B', 'B_[(a)]', EdgeWith.internal, EdgeType.interaction),
                              ('A', 'A_[(a)]', EdgeWith.external, EdgeType.modification),
                              ('A', 'B_[(a)]', EdgeWith.external, EdgeType.modification),
                              ('A_[b]', 'B_[a]', EdgeWith.external, EdgeType.interaction),
                              ('A_[a1]', 'A_[a2]', EdgeWith.external, EdgeType.interaction),
                              ('A_[(a)]', 'B_[(b)]', EdgeWith.external, EdgeType.bimodification)])

    assert _is_graph_test_case_correct(_create_reaction_graph(test_case.quick_string), test_case)