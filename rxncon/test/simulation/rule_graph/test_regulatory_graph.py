
from networkx import DiGraph
from collections import namedtuple
import rxncon.input.quick.quick as qui
from rxncon.simulation.rule_graph.regulatory_graph import RegulatoryGraph, NodeType, EdgeInteractionType

from rxncon.simulation.rule_graph.graphML import map_layout2xgmml, XGMML

RuleTestCase = namedtuple('RuleTestCase', ['quick_string', 'reaction_node_strings', 'state_node_strings',
                                           'boolean_state_node_tuple', 'edge_tuples'])


def _get_state_nodes(test_case, expected_graph):
    """
    Adding state nodes to the expected graph.

    Args:
        test_case: Holding information of the current test.
        expected_graph: The graph defined in the test_case.

    Returns:
        expected graph.

    """
    for node in test_case.state_node_strings:
        if '#in' in node:
            expected_graph.add_node(node.split('#in')[0], type=NodeType.input.value, label=node.split('#in')[0])
        else:
            expected_graph.add_node(node, type=NodeType.state.value, label=node)
    return expected_graph


def _get_reaction_nodes(test_case, expected_graph):
    for reaction_id in test_case.reaction_node_strings:

        if '#out' in reaction_id:
            expected_graph.add_node(reaction_id.split('#out')[0], type=NodeType.output.value)
        else:
            reaction_label = reaction_id.split('#')[0]
            expected_graph.add_node(reaction_id, type=NodeType.reaction.value, label=reaction_label)
    return expected_graph

def _get_boolean_complex_state_nodes(test_case, expected_graph):
    """
    Adding boolean nodes to the expected graph.

    Args:
        test_case: Holding information of the current test.
        expected_graph: The graph defined in the test_case.

    Returns:
        Expected graph.
    """
    if not test_case.boolean_state_node_tuple:
        return expected_graph

    for id, label, type in test_case.boolean_state_node_tuple:

        if 'AND' == type:
            expected_graph.add_node(id, type=NodeType.AND.value, label=label)
        elif 'OR' == type:
            expected_graph.add_node(id, type=NodeType.OR.value, label=label)
        elif 'NOT' == type:
            expected_graph.add_node(id, type=NodeType.NOT.value, label=label)
        else:
            raise AssertionError
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

    expected_graph = _get_reaction_nodes(test_case, expected_graph)


    expected_graph = _get_state_nodes(test_case, expected_graph)
    expected_graph = _get_boolean_complex_state_nodes(test_case, expected_graph)

    [expected_graph.add_edge(source, target, interaction=interaction.value) for source, target, interaction in test_case.edge_tuples]

    return expected_graph.node == actual_graph.node and expected_graph.edge == actual_graph.edge


def _create_regulatory_graph(quick_string):
    """
    Creating a regulatory graph.

    Args:
        quick_string: A rxncon system in quick format.

    Returns:
        A regulatory graph.

    """
    actual_system = qui.Quick(quick_string)
    reg_system = RegulatoryGraph(actual_system.rxncon_system)
    return reg_system.to_graph()


def test_regulatory_graph_for_two_reactions_one_contingency():
    """
    Testing 2 reactions and 1 contingency.

    Returns:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RuleTestCase('''A_[b]_ppi+_B_[a]; ! A-{p}
                                C_p+_A_[(c)]''',
                             ['A_[b]_ppi+_B_[a]', 'C_p+_A_[(c)]'],
                             ['A_[b]--B_[a]', 'A_[(c)]-{p}', 'B_[a]--0', 'A_[b]--0', 'A_[(c)]-{0}'],
                             [],
                             [('A_[b]_ppi+_B_[a]', 'A_[b]--B_[a]', EdgeInteractionType.produce),
                              ('A_[b]_ppi+_B_[a]', 'B_[a]--0', EdgeInteractionType.consume),
                              ('A_[b]_ppi+_B_[a]', 'A_[b]--0', EdgeInteractionType.consume),
                              ('B_[a]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[b]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{0}', EdgeInteractionType.consume),
                              ('A_[(c)]-{0}', 'C_p+_A_[(c)]', EdgeInteractionType.source_state),
                              ('A_[(c)]-{p}', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.required)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_regulatory_graph_for_synthesis_reaction():
    """
    Testing regulatory graph for synthesis reaction.

    Note:
        3 reaction, 0 Contingencies.

    Returns:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RuleTestCase('''C_syn_A
                                A_[b]_ppi+_B_[a]
                                C_p+_A_[(c)]''',
                             ['A_[b]_ppi+_B_[a]', 'C_p+_A_[(c)]', 'C_syn_A'],
                             ['A_[b]--B_[a]', 'A_[(c)]-{p}', 'B_[a]--0', 'A_[b]--0', 'A_[(c)]-{0}'],
                             [],
                             [('A_[b]_ppi+_B_[a]', 'A_[b]--B_[a]', EdgeInteractionType.produce),
                              ('A_[b]_ppi+_B_[a]', 'B_[a]--0', EdgeInteractionType.consume),
                              ('A_[b]_ppi+_B_[a]', 'A_[b]--0', EdgeInteractionType.consume),
                              ('B_[a]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[b]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('C_syn_A', 'A_[b]--0', EdgeInteractionType.produce),
                              ('C_syn_A', 'A_[(c)]-{0}', EdgeInteractionType.produce),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{0}', EdgeInteractionType.consume),
                              ('A_[(c)]-{0}', 'C_p+_A_[(c)]', EdgeInteractionType.source_state)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_regulatory_graph_for_boolean_AND():
    """
    Testing a regulatory graph with boolean AND.

    Note:
        A system of 4 reactions, 1 boolean contingency, 1 state contingency.

    Returns:
        AssertionError: If generated graph differs from expected graph.

    """

    test_case = RuleTestCase('''A_[b]_ppi+_B_[a]; ! <comp>; ! C-{p}
                                <comp>; AND A_[(c)]-{p}; AND A_[c]--C_[a]
                                A_[c]_ppi+_C_[a]
                                C_p+_A_[(c)]
                                D_p+_C_[(d)]''',
                             ['A_[b]_ppi+_B_[a]', 'A_[c]_ppi+_C_[a]', 'C_p+_A_[(c)]', 'D_p+_C_[(d)]'],
                             ['A_[b]--B_[a]', 'A_[b]--0', 'B_[a]--0', 'A_[c]--C_[a]', 'A_[c]--0', 'C_[a]--0', 'A_[(c)]-{p}',
                              'A_[(c)]-{0}', 'C_[(d)]-{p}', 'C_[(d)]-{0}'],
                             [('comp', 'comp', 'AND')],
                             [('A_[b]_ppi+_B_[a]', 'A_[b]--B_[a]', EdgeInteractionType.produce),
                              ('A_[b]_ppi+_B_[a]', 'A_[b]--0', EdgeInteractionType.consume),
                              ('A_[b]_ppi+_B_[a]', 'B_[a]--0', EdgeInteractionType.consume),
                              ('A_[b]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('B_[a]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--C_[a]', EdgeInteractionType.produce),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--0', EdgeInteractionType.consume),
                              ('A_[c]_ppi+_C_[a]', 'C_[a]--0', EdgeInteractionType.consume),
                              ('A_[c]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_[a]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{0}', EdgeInteractionType.consume),
                              ('A_[(c)]-{0}', 'C_p+_A_[(c)]', EdgeInteractionType.source_state),
                              ('D_p+_C_[(d)]', 'C_[(d)]-{p}', EdgeInteractionType.produce),
                              ('D_p+_C_[(d)]', 'C_[(d)]-{0}', EdgeInteractionType.consume),
                              ('C_[(d)]-{0}', 'D_p+_C_[(d)]', EdgeInteractionType.source_state),
                              ('C_[(d)]-{p}', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.required),
                              ('comp', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.required),
                              ('A_[(c)]-{p}', 'comp', EdgeInteractionType.AND),
                              ('A_[c]--C_[a]', 'comp', EdgeInteractionType.AND)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_regulatory_graph_for_boolean_AND_OR():
    """
    Testing a regulatory graph with boolean AND OR combination.

    Note:
        A system of 6 reactions, 1 boolean contingency.

    Returns:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RuleTestCase('''A_[b]_ppi+_B_[a]; ! <comp>
                                <comp>; AND <comp1>; AND <comp2>
                                <comp1>; OR <comp3>; OR A--C
                                <comp2>; AND A_[d]--D_[a]; AND A--E
                                <comp3>; AND A_[f]--F_[a]; AND A--G
                                A_[c]_ppi+_C_[a]
                                A_[d]_ppi+_D_[a]
                                A_[e]_ppi+_E_[a]
                                A_[f]_ppi+_F_[a]
                                A_[g]_ppi+_G_[a]''',
                             ['A_[b]_ppi+_B_[a]', 'A_[c]_ppi+_C_[a]', 'A_[d]_ppi+_D_[a]', 'A_[e]_ppi+_E_[a]', 'A_[f]_ppi+_F_[a]',
                              'A_[g]_ppi+_G_[a]'],
                             ['A_[b]--B_[a]', 'A_[b]--0', 'B_[a]--0', 'A_[c]--C_[a]', 'A_[c]--0', 'C_[a]--0', 'A_[d]--D_[a]',
                              'A_[d]--0', 'D_[a]--0', 'A_[e]--E_[a]', 'A_[e]--0', 'E_[a]--0', 'A_[f]--F_[a]', 'A_[f]--0',
                              'F_[a]--0',
                              'A_[g]--G_[a]', 'A_[g]--0', 'G_[a]--0'],
                             [('comp', 'comp', 'AND'), ('comp1', 'comp1', 'OR'), ('comp2','comp2', 'AND'), ('comp3','comp3', 'AND')],
                             [('A_[b]_ppi+_B_[a]', 'A_[b]--B_[a]', EdgeInteractionType.produce),
                              ('A_[b]_ppi+_B_[a]', 'A_[b]--0', EdgeInteractionType.consume),
                              ('A_[b]_ppi+_B_[a]', 'B_[a]--0', EdgeInteractionType.consume),
                              ('A_[b]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('B_[a]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--C_[a]', EdgeInteractionType.produce),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--0', EdgeInteractionType.consume),
                              ('A_[c]_ppi+_C_[a]', 'C_[a]--0', EdgeInteractionType.consume),
                              ('A_[c]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_[a]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('A_[d]_ppi+_D_[a]', 'A_[d]--D_[a]', EdgeInteractionType.produce),
                              ('A_[d]_ppi+_D_[a]', 'A_[d]--0', EdgeInteractionType.consume),
                              ('A_[d]_ppi+_D_[a]', 'D_[a]--0', EdgeInteractionType.consume),
                              ('A_[d]--0', 'A_[d]_ppi+_D_[a]', EdgeInteractionType.source_state),
                              ('D_[a]--0', 'A_[d]_ppi+_D_[a]', EdgeInteractionType.source_state),
                              ('A_[e]_ppi+_E_[a]', 'A_[e]--E_[a]', EdgeInteractionType.produce),
                              ('A_[e]_ppi+_E_[a]', 'A_[e]--0', EdgeInteractionType.consume),
                              ('A_[e]_ppi+_E_[a]', 'E_[a]--0', EdgeInteractionType.consume),
                              ('A_[e]--0', 'A_[e]_ppi+_E_[a]', EdgeInteractionType.source_state),
                              ('E_[a]--0', 'A_[e]_ppi+_E_[a]', EdgeInteractionType.source_state),
                              ('A_[f]_ppi+_F_[a]', 'A_[f]--F_[a]', EdgeInteractionType.produce),
                              ('A_[f]_ppi+_F_[a]', 'A_[f]--0', EdgeInteractionType.consume),
                              ('A_[f]_ppi+_F_[a]', 'F_[a]--0', EdgeInteractionType.consume),
                              ('A_[f]--0', 'A_[f]_ppi+_F_[a]', EdgeInteractionType.source_state),
                              ('F_[a]--0', 'A_[f]_ppi+_F_[a]', EdgeInteractionType.source_state),
                              ('A_[g]_ppi+_G_[a]', 'A_[g]--G_[a]', EdgeInteractionType.produce),
                              ('A_[g]_ppi+_G_[a]', 'A_[g]--0', EdgeInteractionType.consume),
                              ('A_[g]_ppi+_G_[a]', 'G_[a]--0', EdgeInteractionType.consume),
                              ('A_[g]--0', 'A_[g]_ppi+_G_[a]', EdgeInteractionType.source_state),
                              ('G_[a]--0', 'A_[g]_ppi+_G_[a]', EdgeInteractionType.source_state),
                              ('comp', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.required),
                              ('comp1', 'comp', EdgeInteractionType.AND),
                              ('comp2', 'comp', EdgeInteractionType.AND),
                              ('comp3', 'comp1', EdgeInteractionType.OR),
                              ('A_[c]--C_[a]', 'comp1', EdgeInteractionType.OR),
                              ('A_[d]--D_[a]', 'comp2', EdgeInteractionType.AND),
                              ('A_[e]--E_[a]', 'comp2', EdgeInteractionType.AND),
                              ('A_[f]--F_[a]', 'comp3', EdgeInteractionType.AND),
                              ('A_[g]--G_[a]', 'comp3', EdgeInteractionType.AND), ])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_regulatory_graph_for_boolean_AND_OR_NOT():
    """
    Testing regulatory graph for boolean AND OR NOT combination.

    Note:
        A system of 5 reactions, 1 boolean contingency.

    Returns:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RuleTestCase('''A_[b]_ppi+_B_[a]; ! <comp>
                                <comp>; AND <comp1>; AND <Notcomp2>
                                <comp1>; OR <comp3>; OR A--C
                                <Notcomp2>; NOT A--D
                                <comp3>; AND A--F; AND A--G
                                A_[c]_ppi+_C_[a]
                                A_[d]_ppi+_D_[a]
                                A_[f]_ppi+_F_[a]
                                A_[g]_ppi+_G_[a]''',
                             ['A_[b]_ppi+_B_[a]', 'A_[c]_ppi+_C_[a]', 'A_[d]_ppi+_D_[a]', 'A_[f]_ppi+_F_[a]', 'A_[g]_ppi+_G_[a]'],
                             ['A_[b]--B_[a]', 'A_[b]--0', 'B_[a]--0', 'A_[c]--C_[a]', 'A_[c]--0', 'C_[a]--0', 'A_[d]--D_[a]',
                              'A_[d]--0', 'D_[a]--0', 'A_[f]--F_[a]', 'A_[f]--0', 'F_[a]--0', 'A_[g]--G_[a]', 'A_[g]--0',
                              'G_[a]--0'],
                             [('comp', 'comp', 'AND'), ('comp1', 'comp1', 'OR'), ('Notcomp2', 'Notcomp2', 'NOT'), ('comp3', 'comp3', 'AND')],
                             [('A_[b]_ppi+_B_[a]', 'A_[b]--B_[a]', EdgeInteractionType.produce),
                              ('A_[b]_ppi+_B_[a]', 'A_[b]--0', EdgeInteractionType.consume),
                              ('A_[b]_ppi+_B_[a]', 'B_[a]--0', EdgeInteractionType.consume),
                              ('A_[b]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('B_[a]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--C_[a]', EdgeInteractionType.produce),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--0', EdgeInteractionType.consume),
                              ('A_[c]_ppi+_C_[a]', 'C_[a]--0', EdgeInteractionType.consume),
                              ('A_[c]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_[a]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('A_[d]_ppi+_D_[a]', 'A_[d]--D_[a]', EdgeInteractionType.produce),
                              ('A_[d]_ppi+_D_[a]', 'A_[d]--0', EdgeInteractionType.consume),
                              ('A_[d]_ppi+_D_[a]', 'D_[a]--0', EdgeInteractionType.consume),
                              ('A_[d]--0', 'A_[d]_ppi+_D_[a]', EdgeInteractionType.source_state),
                              ('D_[a]--0', 'A_[d]_ppi+_D_[a]', EdgeInteractionType.source_state),
                              ('A_[f]_ppi+_F_[a]', 'A_[f]--F_[a]', EdgeInteractionType.produce),
                              ('A_[f]_ppi+_F_[a]', 'A_[f]--0', EdgeInteractionType.consume),
                              ('A_[f]_ppi+_F_[a]', 'F_[a]--0', EdgeInteractionType.consume),
                              ('A_[f]--0', 'A_[f]_ppi+_F_[a]', EdgeInteractionType.source_state),
                              ('F_[a]--0', 'A_[f]_ppi+_F_[a]', EdgeInteractionType.source_state),
                              ('A_[g]_ppi+_G_[a]', 'A_[g]--G_[a]', EdgeInteractionType.produce),
                              ('A_[g]_ppi+_G_[a]', 'A_[g]--0', EdgeInteractionType.consume),
                              ('A_[g]_ppi+_G_[a]', 'G_[a]--0', EdgeInteractionType.consume),
                              ('A_[g]--0', 'A_[g]_ppi+_G_[a]', EdgeInteractionType.source_state),
                              ('G_[a]--0', 'A_[g]_ppi+_G_[a]', EdgeInteractionType.source_state),
                              ('comp', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.required),
                              ('comp1', 'comp', EdgeInteractionType.AND),
                              ('Notcomp2', 'comp', EdgeInteractionType.AND),
                              ('comp3', 'comp1', EdgeInteractionType.OR),
                              ('A_[c]--C_[a]', 'comp1', EdgeInteractionType.OR),
                              ('A_[d]--D_[a]', 'Notcomp2', EdgeInteractionType.NOT),
                              ('A_[f]--F_[a]', 'comp3', EdgeInteractionType.AND),
                              ('A_[g]--G_[a]', 'comp3', EdgeInteractionType.AND), ])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_regulatory_graph_for_production_consumption():
    """
    Testing a regulatory graph for production and consumption reactions.

    Note:
        A system of 2 reactions, 0 contingencies.

    Returns:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RuleTestCase('''A_p+_B_[(a)]
                                C_p-_B_[(a)]''',
                             ['A_p+_B_[(a)]', 'C_p-_B_[(a)]'],
                             ['B_[(a)]-{p}', 'B_[(a)]-{0}'],
                             [],
                             [('A_p+_B_[(a)]', 'B_[(a)]-{p}', EdgeInteractionType.produce),
                              ('A_p+_B_[(a)]', 'B_[(a)]-{0}', EdgeInteractionType.consume),
                              ('B_[(a)]-{0}', 'A_p+_B_[(a)]', EdgeInteractionType.source_state),
                              ('C_p-_B_[(a)]', 'B_[(a)]-{p}', EdgeInteractionType.consume),
                              ('B_[(a)]-{p}', 'C_p-_B_[(a)]', EdgeInteractionType.source_state),
                              ('C_p-_B_[(a)]', 'B_[(a)]-{0}', EdgeInteractionType.produce)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_regulatory_graph_for_INPUT():
    """
    Testing a regulatory graph for INPUT contingency.

    Note:
        A system of 1 reaction, 1 input contingency.

    Returns:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RuleTestCase('''A_[B]_ppi+_B_[A]; ! [Input]''',
                             ['A_[B]_ppi+_B_[A]'],
                             ['A_[B]--B_[A]', 'A_[B]--0', 'B_[A]--0', '[Input]#in'],
                             [],
                             [('A_[B]_ppi+_B_[A]', 'A_[B]--B_[A]', EdgeInteractionType.produce),
                              ('A_[B]_ppi+_B_[A]', 'A_[B]--0', EdgeInteractionType.consume),
                              ('A_[B]_ppi+_B_[A]', 'B_[A]--0', EdgeInteractionType.consume),
                              ('A_[B]--0', 'A_[B]_ppi+_B_[A]', EdgeInteractionType.source_state),
                              ('B_[A]--0', 'A_[B]_ppi+_B_[A]', EdgeInteractionType.source_state),
                              ('[Input]', 'A_[B]_ppi+_B_[A]', EdgeInteractionType.required)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_regulatory_graph_for_boolean_with_INPUT():
    """
    Testing a regulatory graph for a boolean contingency containing an input.

    Note:
        A system of 3 reactions, 1 boolean contingency.

    Returns:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RuleTestCase('''A_ppi+_B; ! <comp>
                                <comp>; AND <comp1>; AND [Input]
                                <comp1>; OR A--D; OR A--C
                                A_ppi+_D
                                A_ppi+_C''',
                             ['A_[B]_ppi+_B_[A]', 'A_[C]_ppi+_C_[A]', 'A_[D]_ppi+_D_[A]'],
                             ['A_[B]--B_[A]', 'A_[B]--0', 'B_[A]--0', 'A_[C]--C_[A]', 'A_[C]--0', 'C_[A]--0', 'A_[D]--D_[A]',
                              'A_[D]--0', 'D_[A]--0', '[Input]#in'],
                             [('comp', 'comp', 'AND'), ('comp1', 'comp1', 'OR')],
                             [('A_[B]_ppi+_B_[A]', 'A_[B]--B_[A]', EdgeInteractionType.produce),
                              ('A_[B]_ppi+_B_[A]', 'A_[B]--0', EdgeInteractionType.consume),
                              ('A_[B]_ppi+_B_[A]', 'B_[A]--0', EdgeInteractionType.consume),
                              ('A_[B]--0', 'A_[B]_ppi+_B_[A]', EdgeInteractionType.source_state),
                              ('B_[A]--0', 'A_[B]_ppi+_B_[A]', EdgeInteractionType.source_state),
                              ('A_[C]_ppi+_C_[A]', 'A_[C]--C_[A]', EdgeInteractionType.produce),
                              ('A_[C]_ppi+_C_[A]', 'A_[C]--0', EdgeInteractionType.consume),
                              ('A_[C]_ppi+_C_[A]', 'C_[A]--0', EdgeInteractionType.consume),
                              ('A_[C]--0', 'A_[C]_ppi+_C_[A]', EdgeInteractionType.source_state),
                              ('C_[A]--0', 'A_[C]_ppi+_C_[A]', EdgeInteractionType.source_state),
                              ('A_[D]_ppi+_D_[A]', 'A_[D]--D_[A]', EdgeInteractionType.produce),
                              ('A_[D]_ppi+_D_[A]', 'A_[D]--0', EdgeInteractionType.consume),
                              ('A_[D]_ppi+_D_[A]', 'D_[A]--0', EdgeInteractionType.consume),
                              ('A_[D]--0', 'A_[D]_ppi+_D_[A]', EdgeInteractionType.source_state),
                              ('D_[A]--0', 'A_[D]_ppi+_D_[A]', EdgeInteractionType.source_state),
                              ('comp', 'A_[B]_ppi+_B_[A]', EdgeInteractionType.required),
                              ('[Input]', 'comp', EdgeInteractionType.AND),
                              ('comp1', 'comp', EdgeInteractionType.AND),
                              ('A_[D]--D_[A]', 'comp1', EdgeInteractionType.OR),
                              ('A_[C]--C_[A]', 'comp1', EdgeInteractionType.OR)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_regulatory_graph_for_structured_boolean():
    """
    Testing a regulatory graph for structured boolean contingency.

    Note:
        A system of 2 reactions, 1 boolean contingency.

    Returns:
        AssertionError: If generated graph differs from expected graph.

    """

    test_case = RuleTestCase("""IR_[lig]_i+_insulin_[IR]; ! <IR-empty>
                                <IR-empty>; AND IR@0--IR@2; AND IR@0_[lig]--0; AND IR@2_[lig]--0
                                IR_[IR]_ppi+_IR_[IR]""",
                             ['IR_[lig]_i+_insulin_[IR]', 'IR_[IR]_ppi+_IR_[IR]'],
                             ['IR_[IR]--IR_[IR]', 'IR_[IR]--0', 'IR_[lig]--insulin_[IR]', 'IR_[lig]--0',
                              'insulin_[IR]--0'],
                             [('IR-empty', 'IR-empty', 'AND')],
                             [('IR_[IR]_ppi+_IR_[IR]', 'IR_[IR]--IR_[IR]', EdgeInteractionType.produce),
                              ('IR_[IR]_ppi+_IR_[IR]', 'IR_[IR]--0', EdgeInteractionType.consume),
                             ('IR_[IR]--0', 'IR_[IR]_ppi+_IR_[IR]', EdgeInteractionType.source_state),
                             ('IR_[lig]_i+_insulin_[IR]', 'IR_[lig]--insulin_[IR]', EdgeInteractionType.produce),
                             ('IR_[lig]_i+_insulin_[IR]', 'IR_[lig]--0', EdgeInteractionType.consume),
                             ('IR_[lig]_i+_insulin_[IR]', 'insulin_[IR]--0', EdgeInteractionType.consume),
                             ('IR_[lig]--0', 'IR_[lig]_i+_insulin_[IR]', EdgeInteractionType.source_state),
                             ('insulin_[IR]--0', 'IR_[lig]_i+_insulin_[IR]', EdgeInteractionType.source_state),
                             ('IR-empty', 'IR_[lig]_i+_insulin_[IR]', EdgeInteractionType.required),
                             ('IR_[IR]--IR_[IR]', 'IR-empty', EdgeInteractionType.AND),
                             ('IR_[lig]--0', 'IR-empty', EdgeInteractionType.AND),
                             ('IR_[lig]--0', 'IR-empty', EdgeInteractionType.AND)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)

def test_regulatory_graph_for_degradation_no_contingency():
    test_case = RuleTestCase('''A_[b]_ppi+_B_[a]
                                C_p+_A_[(c)]
                                D_deg_A''',
<<<<<<< ffbf8228ef51283ed384c7258021bad0a79472b2

=======
>>>>>>> added first degradation reaction and added test.
                             ['A_[b]_ppi+_B_[a]', 'C_p+_A_[(c)]', 'D_deg_A'],
                             ['A_[b]--B_[a]', 'A_[(c)]-{p}', 'B_[a]--0', 'A_[b]--0', 'A_[(c)]-{0}'],
                             [('D_deg_A_AND_A_[b]--B_[a]', " ", 'AND')],
                             [('A_[b]_ppi+_B_[a]', 'A_[b]--B_[a]', EdgeInteractionType.produce),
                              ('A_[b]_ppi+_B_[a]', 'B_[a]--0', EdgeInteractionType.consume),
                              ('A_[b]_ppi+_B_[a]', 'A_[b]--0', EdgeInteractionType.consume),
                              ('B_[a]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[b]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{0}', EdgeInteractionType.consume),
                              ('A_[(c)]-{0}', 'C_p+_A_[(c)]', EdgeInteractionType.source_state),

                              ('D_deg_A', 'A_[b]--B_[a]', EdgeInteractionType.consume),
                              ('D_deg_A', 'A_[b]--0', EdgeInteractionType.consume),
                              ('D_deg_A', 'A_[(c)]-{p}', EdgeInteractionType.consume),
                              ('D_deg_A', 'A_[(c)]-{0}', EdgeInteractionType.consume),
                              ('D_deg_A', 'D_deg_A_AND_A_[b]--B_[a]', EdgeInteractionType.AND),
                              ('A_[b]--B_[a]', 'D_deg_A_AND_A_[b]--B_[a]', EdgeInteractionType.AND),
                              ('D_deg_A_AND_A_[b]--B_[a]', 'B_[a]--0', EdgeInteractionType.produce)])


    reg_graph = _create_regulatory_graph(test_case.quick_string)
    gml_system = XGMML(reg_graph, "reactions_only")
#    gml_system.to_file("test_deg.xgmml")
    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)

def test_degradation_with_contingency():
    test_case = RuleTestCase('''A_[b]_ppi+_B_[a]
                                A_[c]_ppi+_C_[a]
                                C_p+_A_[(c)]
                                C_p+_A_[(d)]
                                D_deg_A; ! <comp>
                                <comp>; AND A_[(c)]-{p}; AND A_[b]--B_[a]''',
                             ['A_[b]_ppi+_B_[a]', 'A_[c]_ppi+_C_[a]', 'C_p+_A_[(c)]', 'C_p+_A_[(d)]', 'D_deg_A#0'],
                             ['A_[b]--B_[a]', 'A_[c]--C_[a]', 'A_[(c)]-{p}', 'A_[(d)]-{p}', 'B_[a]--0', 'C_[a]--0',
                              'A_[b]--0', 'A_[c]--0', 'A_[(c)]-{0}', 'A_[(d)]-{0}',],
                             [('D_deg_A#0_AND_A_[b]--B_[a]', " ", 'AND'), ('comp', 'comp', 'AND')],
                             [('A_[b]_ppi+_B_[a]', 'A_[b]--B_[a]', EdgeInteractionType.produce),
                              ('A_[b]_ppi+_B_[a]', 'B_[a]--0', EdgeInteractionType.consume),
                              ('A_[b]_ppi+_B_[a]', 'A_[b]--0', EdgeInteractionType.consume),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--C_[a]', EdgeInteractionType.produce),
                              ('A_[c]_ppi+_C_[a]', 'C_[a]--0', EdgeInteractionType.consume),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--0', EdgeInteractionType.consume),
                              ('B_[a]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[b]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[c]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_[a]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{0}', EdgeInteractionType.consume),
                              ('A_[(c)]-{0}', 'C_p+_A_[(c)]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(d)]', 'A_[(d)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(d)]', 'A_[(d)]-{0}', EdgeInteractionType.consume),
                              ('A_[(d)]-{0}', 'C_p+_A_[(d)]', EdgeInteractionType.source_state),
                              ('D_deg_A#0', 'A_[b]--B_[a]', EdgeInteractionType.consume),
                              ('D_deg_A#0', 'A_[(c)]-{p}', EdgeInteractionType.consume),
                              ('D_deg_A#0', 'D_deg_A#0_AND_A_[b]--B_[a]', EdgeInteractionType.AND),
                              ('A_[b]--B_[a]', 'D_deg_A#0_AND_A_[b]--B_[a]', EdgeInteractionType.AND),
                              ('D_deg_A#0_AND_A_[b]--B_[a]', 'B_[a]--0', EdgeInteractionType.produce),
                              ('A_[(c)]-{p}', 'comp', EdgeInteractionType.AND),
                              ('A_[b]--B_[a]', 'comp', EdgeInteractionType.AND),
                              ('comp', 'D_deg_A#0', EdgeInteractionType.required)])

    reg_graph = _create_regulatory_graph(test_case.quick_string)
    gml_system = XGMML(reg_graph, "reactions_only")
    gml_system.to_file("test_deg_bool_cont_AND.xgmml")
    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_degradation_with_boolean_contingency_OR():
    test_case = RuleTestCase('''A_[b]_ppi+_B_[a]
                                A_[b]_ppi+_C_[a]
                                C_p+_A_[(c)]
                                C_p+_A_[(d)]
                                D_deg_A; ! <comp>
                                <comp>; OR A_[(c)]-{p}; OR A_[b]--B_[a]''',
                             [],
                             [],
                             [],
                             [])

    reg_graph = _create_regulatory_graph(test_case.quick_string)
    gml_system = XGMML(reg_graph, "reactions_only")
    #gml_system.to_file("test_deg_bool_cont_OR.xgmml")
    #assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)



def test_degradation_with_contingency2():
    test_case = RuleTestCase('''A_[b]_ppi+_B_[a]
                                A_[b]_ppi+_C_[a]
                                C_p+_A_[(c)]
                                C_p+_A_[(d)]
                                D_deg_A; ! <comp>
                                <comp>; AND <comp1>; AND A_[(c)]-{p}
                                <comp1>; OR <comp2>; OR  A_[b]--B_[a]
                                <comp2>; AND A_[(d)]-{p}; AND AND A_[b]--C_[a]''',
                             [],
                             [],
                             [],
                             [])

    _create_regulatory_graph(test_case.quick_string)
    reg_graph = _create_regulatory_graph(test_case.quick_string)
    gml_system = XGMML(reg_graph, "reactions_only")
    #gml_system.to_file("test_deg_bool_cont_AND_OR.xgmml")
    #assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)
