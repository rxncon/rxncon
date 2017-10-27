from collections import namedtuple

from networkx import DiGraph

import rxncon.input.quick.quick as qui
import rxncon.input.excel_book.excel_book as excel

import rxncon.visualization.graphML as graphML

from rxncon.visualization.regulatory_graph import SpeciesReactionGraph, NodeType, EdgeInteractionType, RegulatoryGraph


RegulatoryGraphTestCase = namedtuple('RegulatoryGraphTestCase',
                                     ['quick_string', 'reaction_node_strings', 'state_node_strings',
                                      'boolean_state_node_tuple', 'edge_tuples'])


def test_simple_system() -> None:
    """
    Testing 2 reactions and 1 contingency.

    Returns:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_[b]_ppi+_B_[a]; ! A-{p}
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


def test_production_consumption() -> None:
    """
    Testing a regulatory graph for production and consumption reactions.

    Note:
        A system of 2 reactions, 0 contingencies.

    Returns:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_p+_B_[(a)]
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


def test_input_state() -> None:
    """
    Testing a regulatory graph for INPUT contingency.

    Note:
        A system of 1 reaction, 1 input contingency.

    Returns:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_[B]_ppi+_B_[A]; ! [Input]''',
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


def test_boolean_AND() -> None:
    """
    Testing a regulatory graph with boolean AND.

    Note:
        A system of 4 reactions, 1 boolean contingency, 1 state contingency.

    Returns:
        AssertionError: If generated graph differs from expected graph.

    """

    test_case = RegulatoryGraphTestCase('''A_[b]_ppi+_B_[a]; ! <comp>; ! C-{p}
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


def test_boolean_AND_OR() -> None:
    """
    Testing a regulatory graph with boolean AND OR combination.

    Note:
        A system of 6 reactions, 1 boolean contingency.

    Returns:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_[b]_ppi+_B_[a]; ! <comp>
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


def test_boolean_AND_OR_NOT() -> None:
    """
    Testing regulatory graph for boolean AND OR NOT combination.

    Note:
        A system of 5 reactions, 1 boolean contingency.

    Returns:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_[b]_ppi+_B_[a]; ! <comp>
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


def test_boolean_with_input_state() -> None:
    """
    Testing a regulatory graph for a boolean contingency containing an input.

    Note:
        A system of 3 reactions, 1 boolean contingency.

    Returns:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_ppi+_B; ! <comp>
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


def test_degradation_no_contingency() -> None:
    """
    Testing a regulatory graph for degradation reaction without contingencies.

    Note:
        A system of 3 reactions.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_[b]_ppi+_B_[a]
                                C_p+_A_[(c)]
                                D_deg_A''',
                                        ['A_[b]_ppi+_B_[a]', 'C_p+_A_[(c)]', 'D_deg_A'],
                                        ['A_[b]--B_[a]', 'A_[(c)]-{p}', 'B_[a]--0', 'A_[b]--0', 'A_[(c)]-{0}'],
                                        [('D_deg_A_ON_A_[b]--B_[a]', " ", 'AND')],

                                        [('A_[b]_ppi+_B_[a]', 'A_[b]--B_[a]', EdgeInteractionType.produce),
                              ('A_[b]_ppi+_B_[a]', 'B_[a]--0', EdgeInteractionType.consume),
                              ('A_[b]_ppi+_B_[a]', 'A_[b]--0', EdgeInteractionType.consume),
                              ('B_[a]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[b]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{0}', EdgeInteractionType.consume),
                              ('A_[(c)]-{0}', 'C_p+_A_[(c)]', EdgeInteractionType.source_state),


                              ('D_deg_A', 'A_[b]--B_[a]', EdgeInteractionType.degrade),
                              ('D_deg_A', 'A_[b]--0', EdgeInteractionType.degrade),
                              ('D_deg_A', 'A_[(c)]-{p}', EdgeInteractionType.degrade),
                              ('D_deg_A', 'A_[(c)]-{0}', EdgeInteractionType.degrade),

                              ('D_deg_A', 'D_deg_A_ON_A_[b]--B_[a]', EdgeInteractionType.AND),
                              ('A_[b]--B_[a]', 'D_deg_A_ON_A_[b]--B_[a]', EdgeInteractionType.AND),
                              ('D_deg_A_ON_A_[b]--B_[a]', 'B_[a]--0', EdgeInteractionType.produce)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_degradation_inhibited_by_interaction() -> None:
    """
    Testing regulatory graph with degradation and contingencies

    Note:
        A system of 4 reactions and one inhibition contingency.
        In this scenario we also create edges to optional/possible degradation targets.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_[b]_ppi+_B_[a]
                                A_[c]_ppi+_C_[a]
                                C_p+_A_[(c)]
                                D_deg_A; x A_[b]--B_[a]''',
                                        ['A_[b]_ppi+_B_[a]', 'A_[c]_ppi+_C_[a]', 'C_p+_A_[(c)]', 'D_deg_A'],
                                        ['A_[b]--B_[a]', 'A_[c]--C_[a]', 'A_[(c)]-{p}', 'B_[a]--0', 'C_[a]--0', 'A_[b]--0', 'A_[c]--0', 'A_[(c)]-{0}'],
                                        [('D_deg_A_ON_A_[c]--C_[a]', " ", 'AND')],

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
                              ('A_[b]--B_[a]', 'D_deg_A', EdgeInteractionType.inhibited),
                              ('D_deg_A', 'A_[b]--0', EdgeInteractionType.degrade),
                              ('D_deg_A', 'A_[(c)]-{p}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[(c)]-{0}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[c]--C_[a]', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[c]--0', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'D_deg_A_ON_A_[c]--C_[a]', EdgeInteractionType.AND),
                              ('A_[c]--C_[a]', 'D_deg_A_ON_A_[c]--C_[a]', EdgeInteractionType.AND),
                              ('D_deg_A_ON_A_[c]--C_[a]', 'C_[a]--0', EdgeInteractionType.produce)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_degradation_inhibited_and_mutually_exclusivity() -> None:
    """
    Testing regulatory graph with degradation and contingencies.

    Note:
        A system of 4 reactions and 1 contingency.
        The complement of A_[(c)]-{P} should be maybe_degraded but not A_[(c)]-{P}.
        Mutually exclusive states should not be degraded.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_[c]_ppi+_C_[a]
                                C_p+_A_[(c)]
                                C_ub+_A_[(c)]
                                D_deg_A; x A_[(c)]-{P}''',
                                        ['A_[c]_ppi+_C_[a]', 'C_p+_A_[(c)]', 'D_deg_A', 'C_ub+_A_[(c)]'],
                                        ['A_[c]--C_[a]', 'A_[(c)]-{p}', 'A_[(c)]-{ub}', 'C_[a]--0', 'A_[c]--0', 'A_[(c)]-{0}'],
                                        [('D_deg_A_ON_A_[c]--C_[a]', " ", 'AND')],

                                        [('A_[c]_ppi+_C_[a]', 'A_[c]--C_[a]', EdgeInteractionType.produce),
                              ('A_[c]_ppi+_C_[a]', 'C_[a]--0', EdgeInteractionType.consume),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--0', EdgeInteractionType.consume),
                              ('A_[c]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_[a]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{0}', EdgeInteractionType.consume),
                              ('C_ub+_A_[(c)]', 'A_[(c)]-{ub}', EdgeInteractionType.produce),
                              ('C_ub+_A_[(c)]', 'A_[(c)]-{0}', EdgeInteractionType.consume),
                              ('A_[(c)]-{0}', 'C_p+_A_[(c)]', EdgeInteractionType.source_state),
                              ('A_[(c)]-{0}', 'C_ub+_A_[(c)]', EdgeInteractionType.source_state),
                              ('A_[(c)]-{p}', 'D_deg_A', EdgeInteractionType.inhibited),
                              ('D_deg_A', 'A_[(c)]-{0}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[(c)]-{ub}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[c]--C_[a]', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[c]--0', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'D_deg_A_ON_A_[c]--C_[a]', EdgeInteractionType.AND),
                              ('A_[c]--C_[a]', 'D_deg_A_ON_A_[c]--C_[a]', EdgeInteractionType.AND),
                              ('D_deg_A_ON_A_[c]--C_[a]', 'C_[a]--0', EdgeInteractionType.produce)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_degradation_with_contingency() -> None:
    """
    Testing regulatory graph with degradation and contingencies.

    Note:
        The system contains 5 reactions and 1 boolean contingency.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_[b]_ppi+_B_[a]
                                A_[c]_ppi+_C_[a]
                                C_p+_A_[(c)]
                                C_p+_A_[(d)]
                                D_deg_A; ! <comp>
                                <comp>; AND A_[(c)]-{p}; AND A_[b]--B_[a]''',
                                        ['A_[b]_ppi+_B_[a]', 'A_[c]_ppi+_C_[a]', 'C_p+_A_[(c)]', 'C_p+_A_[(d)]', 'D_deg_A'],
                                        ['A_[b]--B_[a]', 'A_[c]--C_[a]', 'A_[(c)]-{p}', 'A_[(d)]-{p}', 'B_[a]--0', 'C_[a]--0',
                              'A_[b]--0', 'A_[c]--0', 'A_[(c)]-{0}', 'A_[(d)]-{0}',],
                                        [('D_deg_A_ON_A_[b]--B_[a]', " ", 'AND'), ('D_deg_A_ON_A_[c]--C_[a]', " ", 'AND'),
                              ('comp', 'comp', 'AND')],
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
                              ('D_deg_A', 'A_[b]--B_[a]', EdgeInteractionType.degrade),
                              ('D_deg_A', 'A_[(c)]-{p}', EdgeInteractionType.degrade),
                              ('D_deg_A', 'A_[(d)]-{p}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[(d)]-{0}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[c]--0', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[c]--C_[a]', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'D_deg_A_ON_A_[b]--B_[a]', EdgeInteractionType.AND),
                              ('A_[b]--B_[a]', 'D_deg_A_ON_A_[b]--B_[a]', EdgeInteractionType.AND),
                              ('D_deg_A_ON_A_[b]--B_[a]', 'B_[a]--0', EdgeInteractionType.produce),
                              ('D_deg_A', 'D_deg_A_ON_A_[c]--C_[a]', EdgeInteractionType.AND),
                              ('A_[c]--C_[a]', 'D_deg_A_ON_A_[c]--C_[a]', EdgeInteractionType.AND),
                              ('D_deg_A_ON_A_[c]--C_[a]', 'C_[a]--0', EdgeInteractionType.produce),
                              ('A_[(c)]-{p}', 'comp', EdgeInteractionType.AND),
                              ('A_[b]--B_[a]', 'comp', EdgeInteractionType.AND),
                              ('comp', 'D_deg_A', EdgeInteractionType.required),])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_degradation_boolean_AND_NOT() -> None:
    """
    Testing the degradation of a boolean AND NOT combination.

    Note:

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_[x]_ppi+_B_[a]
                                A_[c]_ppi+_C_[a]
                                C_p+_A_[(c)]
                                C_p+_A_[(d)]
                                D_deg_A; ! <comp>
                                <comp>; AND A_[(c)]-{p}; AND <NOT>
                                <NOT>; NOT A_[x]--B_[a]''',
                                        ['A_[x]_ppi+_B_[a]', 'A_[c]_ppi+_C_[a]', 'C_p+_A_[(c)]', 'C_p+_A_[(d)]', 'D_deg_A'],
                                        ['A_[x]--B_[a]', 'A_[c]--C_[a]', 'A_[(c)]-{p}', 'A_[(d)]-{p}', 'B_[a]--0', 'C_[a]--0',
                              'A_[x]--0', 'A_[c]--0', 'A_[(c)]-{0}', 'A_[(d)]-{0}'],
                                        [('comp', 'comp', 'AND'), ('NOT', 'NOT', 'NOT'), ('D_deg_A_ON_A_[c]--C_[a]', ' ', 'AND')],
                                        [('A_[x]_ppi+_B_[a]', 'A_[x]--B_[a]', EdgeInteractionType.produce),
                              ('A_[x]_ppi+_B_[a]', 'B_[a]--0', EdgeInteractionType.consume),
                              ('A_[x]_ppi+_B_[a]', 'A_[x]--0', EdgeInteractionType.consume),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--C_[a]', EdgeInteractionType.produce),
                              ('A_[c]_ppi+_C_[a]', 'C_[a]--0', EdgeInteractionType.consume),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--0', EdgeInteractionType.consume),
                              ('B_[a]--0', 'A_[x]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[x]--0', 'A_[x]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[c]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_[a]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{0}', EdgeInteractionType.consume),
                              ('A_[(c)]-{0}', 'C_p+_A_[(c)]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(d)]', 'A_[(d)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(d)]', 'A_[(d)]-{0}', EdgeInteractionType.consume),
                              ('A_[(d)]-{0}', 'C_p+_A_[(d)]', EdgeInteractionType.source_state),
                              ('D_deg_A', 'A_[(c)]-{p}', EdgeInteractionType.degrade),
                              ('D_deg_A', 'A_[x]--0', EdgeInteractionType.degrade),
                              ('D_deg_A', 'A_[c]--0', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[(d)]-{p}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[(d)]-{0}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[c]--C_[a]', EdgeInteractionType.maybe_degraded),

                              ('D_deg_A', 'D_deg_A_ON_A_[c]--C_[a]', EdgeInteractionType.AND),
                              ('A_[c]--C_[a]', 'D_deg_A_ON_A_[c]--C_[a]', EdgeInteractionType.AND),
                              ('D_deg_A_ON_A_[c]--C_[a]', 'C_[a]--0', EdgeInteractionType.produce),
                              ('A_[x]--B_[a]', 'NOT', EdgeInteractionType.NOT),
                              ('NOT', 'comp', EdgeInteractionType.AND),
                              ('A_[(c)]-{p}', 'comp', EdgeInteractionType.AND),
                              ('comp', 'D_deg_A', EdgeInteractionType.required)])


    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_degradation_boolean_multiple_complementary_states() -> None:
    """
    Testing the degradation of a boolean AND NOT combination.

    Note:
        A special case is that the complement of A_[x]--B_[a] is A_[x]--0 and A_[x]--C_[a]. This leads to a
        maybe_degraded edge between D_deg_A and A_[x]--0 and D_deg_A and A_[x]--C_[a].

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_[x]_ppi+_B_[a]
                                A_[x]_ppi+_C_[a]
                                C_p+_A_[(c)]
                                C_p+_A_[(d)]
                                D_deg_A; ! <comp>
                                <comp>; AND A_[(c)]-{p}; AND <NOT>
                                <NOT>; NOT A_[x]--B_[a]''',
                                        ['A_[x]_ppi+_B_[a]', 'A_[x]_ppi+_C_[a]', 'C_p+_A_[(c)]', 'C_p+_A_[(d)]', 'D_deg_A'],
                                        ['A_[x]--B_[a]', 'A_[x]--C_[a]', 'A_[(c)]-{p}', 'A_[(d)]-{p}', 'B_[a]--0', 'C_[a]--0',
                              'A_[x]--0', 'A_[(c)]-{0}', 'A_[(d)]-{0}'],
                                        [('comp', 'comp', 'AND'), ('NOT', 'NOT', 'NOT'), ('D_deg_A_ON_A_[x]--C_[a]', ' ', 'AND')],
                                        [('A_[x]_ppi+_B_[a]', 'A_[x]--B_[a]', EdgeInteractionType.produce),
                              ('A_[x]_ppi+_B_[a]', 'B_[a]--0', EdgeInteractionType.consume),
                              ('A_[x]_ppi+_B_[a]', 'A_[x]--0', EdgeInteractionType.consume),
                              ('A_[x]_ppi+_C_[a]', 'A_[x]--C_[a]', EdgeInteractionType.produce),
                              ('A_[x]_ppi+_C_[a]', 'C_[a]--0', EdgeInteractionType.consume),
                              ('A_[x]_ppi+_C_[a]', 'A_[x]--0', EdgeInteractionType.consume),
                              ('B_[a]--0', 'A_[x]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[x]--0', 'A_[x]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[x]--0', 'A_[x]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_[a]--0', 'A_[x]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{0}', EdgeInteractionType.consume),
                              ('A_[(c)]-{0}', 'C_p+_A_[(c)]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(d)]', 'A_[(d)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(d)]', 'A_[(d)]-{0}', EdgeInteractionType.consume),
                              ('A_[(d)]-{0}', 'C_p+_A_[(d)]', EdgeInteractionType.source_state),
                              ('D_deg_A', 'A_[(c)]-{p}', EdgeInteractionType.degrade),
                              ('D_deg_A', 'A_[x]--0', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[(d)]-{p}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[(d)]-{0}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[x]--C_[a]', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'D_deg_A_ON_A_[x]--C_[a]', EdgeInteractionType.AND),
                              ('A_[x]--C_[a]', 'D_deg_A_ON_A_[x]--C_[a]', EdgeInteractionType.AND),
                              ('D_deg_A_ON_A_[x]--C_[a]', 'C_[a]--0', EdgeInteractionType.produce),
                              ('A_[x]--B_[a]', 'NOT', EdgeInteractionType.NOT),
                              ('NOT', 'comp', EdgeInteractionType.AND),
                              ('A_[(c)]-{p}', 'comp', EdgeInteractionType.AND),
                              ('comp', 'D_deg_A', EdgeInteractionType.required)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_degradation_boolean_NOT_AND() -> None:
    """
    Testing degradation with AND combination of NOT boolean.

    Note:
        We have the special case that we have a NOT of ANDs. This will result in OR of NOTs. Meaning that in this case we
        have NOT A_[x]--B_[a] OR NOT A_[(d)]-{p}.

        In principle we are splitting the degradation during the interpretation step if we have to consider ORs.

        Hence, for the case NOT A_[x]--B_[a].
        A_[x]--B_[a] is protected. We get a degraded edge to A_[x]--0 and maybe_degraded edges to
        A_[(d)]-{0}, A_[(d)]-{p}.

        For the case NOT A_[(d)]-{p}. A_[(d)]-{p} is protected. We we get a degraded edge to A_[(d)]-{0} and
        maybe_degraded edges to A_[x]--B_[a], A_[x]--0.

        Since we are not splitting the degradation reaction in the regulatory graph, the degraded edges will be
        transformed into maybe_degraded edges, because it can be that the state is degraded but not in all cases.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_[x]_ppi+_B_[a]
                                A_[c]_ppi+_C_[a]
                                C_p+_A_[(c)]
                                C_p+_A_[(d)]
                                D_deg_A; ! <comp>
                                <comp>; AND A_[(c)]-{p}; AND <NOT>
                                <NOT>; NOT <comp1>
                                <comp1>; AND A_[x]--B_[a]; AND A_[(d)]-{p}''',
                                        ['A_[x]_ppi+_B_[a]', 'A_[c]_ppi+_C_[a]', 'C_p+_A_[(c)]', 'C_p+_A_[(d)]', 'D_deg_A'],
                                        ['A_[x]--B_[a]', 'A_[c]--C_[a]', 'A_[(c)]-{p}', 'A_[(d)]-{p}', 'B_[a]--0', 'C_[a]--0',
                              'A_[x]--0', 'A_[c]--0', 'A_[(c)]-{0}', 'A_[(d)]-{0}'],
                                        [('comp', 'comp', 'AND'), ('comp1', 'comp1', 'AND'), ('NOT', 'NOT', 'NOT'),
                              ('D_deg_A_ON_A_[c]--C_[a]', ' ', 'AND'), ('D_deg_A_ON_A_[x]--B_[a]', ' ', 'AND')],
                                        [('A_[x]_ppi+_B_[a]', 'A_[x]--B_[a]', EdgeInteractionType.produce),
                              ('A_[x]_ppi+_B_[a]', 'B_[a]--0', EdgeInteractionType.consume),
                              ('A_[x]_ppi+_B_[a]', 'A_[x]--0', EdgeInteractionType.consume),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--C_[a]', EdgeInteractionType.produce),
                              ('A_[c]_ppi+_C_[a]', 'C_[a]--0', EdgeInteractionType.consume),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--0', EdgeInteractionType.consume),
                              ('B_[a]--0', 'A_[x]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[x]--0', 'A_[x]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[c]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_[a]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{0}', EdgeInteractionType.consume),
                              ('A_[(c)]-{0}', 'C_p+_A_[(c)]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(d)]', 'A_[(d)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(d)]', 'A_[(d)]-{0}', EdgeInteractionType.consume),
                              ('A_[(d)]-{0}', 'C_p+_A_[(d)]', EdgeInteractionType.source_state),

                              ('D_deg_A', 'A_[(c)]-{p}', EdgeInteractionType.degrade),
                              ('D_deg_A', 'A_[(d)]-{0}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[(d)]-{p}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[x]--0', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[c]--0', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[c]--C_[a]', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'D_deg_A_ON_A_[c]--C_[a]', EdgeInteractionType.AND),
                              ('A_[c]--C_[a]', 'D_deg_A_ON_A_[c]--C_[a]', EdgeInteractionType.AND),
                              ('D_deg_A_ON_A_[c]--C_[a]', 'C_[a]--0', EdgeInteractionType.produce),
                              ('D_deg_A', 'A_[x]--B_[a]', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'D_deg_A_ON_A_[x]--B_[a]', EdgeInteractionType.AND),
                              ('A_[x]--B_[a]', 'D_deg_A_ON_A_[x]--B_[a]', EdgeInteractionType.AND),
                              ('D_deg_A_ON_A_[x]--B_[a]', 'B_[a]--0', EdgeInteractionType.produce),
                              ('A_[x]--B_[a]', 'comp1', EdgeInteractionType.AND),
                              ('A_[(d)]-{p}', 'comp1', EdgeInteractionType.AND),
                              ('comp1', 'NOT', EdgeInteractionType.NOT),
                              ('NOT', 'comp', EdgeInteractionType.AND),
                              ('A_[(c)]-{p}', 'comp', EdgeInteractionType.AND),
                              ('comp', 'D_deg_A', EdgeInteractionType.required)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_degradation_boolean_OR() -> None:
    """
    Testing degradation with OR statement

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_[x]_ppi+_B_[a]
                                A_[c]_ppi+_C_[a]
                                C_p+_A_[(c)]
                                C_p+_A_[(d)]
                                D_deg_A; ! <comp>
                                <comp>; OR A_[(c)]-{p}; OR A_[x]--B_[a]''',
                                        ['A_[x]_ppi+_B_[a]', 'A_[c]_ppi+_C_[a]', 'C_p+_A_[(c)]', 'C_p+_A_[(d)]', 'D_deg_A'],
                                        ['A_[x]--B_[a]', 'A_[c]--C_[a]', 'A_[(c)]-{p}', 'A_[(d)]-{p}', 'B_[a]--0', 'C_[a]--0',
                              'A_[x]--0', 'A_[c]--0', 'A_[(c)]-{0}', 'A_[(d)]-{0}'],
                                        [('comp', 'comp', 'OR'), ('D_deg_A_ON_A_[c]--C_[a]', ' ', 'AND'),
                              ('D_deg_A_ON_A_[x]--B_[a]', ' ', 'AND')],
                                        [('A_[x]_ppi+_B_[a]', 'A_[x]--B_[a]', EdgeInteractionType.produce),
                              ('A_[x]_ppi+_B_[a]', 'B_[a]--0', EdgeInteractionType.consume),
                              ('A_[x]_ppi+_B_[a]', 'A_[x]--0', EdgeInteractionType.consume),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--C_[a]', EdgeInteractionType.produce),
                              ('A_[c]_ppi+_C_[a]', 'C_[a]--0', EdgeInteractionType.consume),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--0', EdgeInteractionType.consume),
                              ('B_[a]--0', 'A_[x]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[x]--0', 'A_[x]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[c]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_[a]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{0}', EdgeInteractionType.consume),
                              ('A_[(c)]-{0}', 'C_p+_A_[(c)]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(d)]', 'A_[(d)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(d)]', 'A_[(d)]-{0}', EdgeInteractionType.consume),
                              ('A_[(d)]-{0}', 'C_p+_A_[(d)]', EdgeInteractionType.source_state),

                              ('D_deg_A', 'A_[(c)]-{p}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[(c)]-{0}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[x]--B_[a]', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[x]--0', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'D_deg_A_ON_A_[x]--B_[a]', EdgeInteractionType.AND),
                              ('A_[x]--B_[a]', 'D_deg_A_ON_A_[x]--B_[a]', EdgeInteractionType.AND),
                              ('D_deg_A_ON_A_[x]--B_[a]', 'B_[a]--0', EdgeInteractionType.produce),
                              ('D_deg_A', 'A_[(d)]-{0}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[(d)]-{p}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[c]--0', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[c]--C_[a]', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'D_deg_A_ON_A_[c]--C_[a]', EdgeInteractionType.AND),
                              ('A_[c]--C_[a]', 'D_deg_A_ON_A_[c]--C_[a]', EdgeInteractionType.AND),
                              ('D_deg_A_ON_A_[c]--C_[a]', 'C_[a]--0', EdgeInteractionType.produce),
                              ('A_[x]--B_[a]', 'comp', EdgeInteractionType.OR),
                              ('A_[(c)]-{p}', 'comp', EdgeInteractionType.OR),
                              ('comp', 'D_deg_A', EdgeInteractionType.required)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_degradation_boolean_NOT_OR() -> None:
    """
    Testing degradation with NOT of ORs

    Note:
        The NOT of ORS will be translated into a requirement of ANDs
    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_[x]_ppi+_B_[a]
                                A_[c]_ppi+_C_[a]
                                C_p+_A_[(c)]
                                C_p+_A_[(d)]
                                D_deg_A; x <comp>
                                <comp>; OR A_[(c)]-{p}; OR A_[x]--B_[a]''',
                                        ['A_[x]_ppi+_B_[a]', 'A_[c]_ppi+_C_[a]', 'C_p+_A_[(c)]', 'C_p+_A_[(d)]', 'D_deg_A'],
                                        ['A_[x]--B_[a]', 'A_[c]--C_[a]', 'A_[(c)]-{p}', 'A_[(d)]-{p}', 'B_[a]--0', 'C_[a]--0',
                              'A_[x]--0', 'A_[c]--0', 'A_[(c)]-{0}', 'A_[(d)]-{0}'],
                                        [('comp', 'comp', 'OR'), ('D_deg_A_ON_A_[c]--C_[a]', ' ', 'AND')],
                                        [('A_[x]_ppi+_B_[a]', 'A_[x]--B_[a]', EdgeInteractionType.produce),
                              ('A_[x]_ppi+_B_[a]', 'B_[a]--0', EdgeInteractionType.consume),
                              ('A_[x]_ppi+_B_[a]', 'A_[x]--0', EdgeInteractionType.consume),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--C_[a]', EdgeInteractionType.produce),
                              ('A_[c]_ppi+_C_[a]', 'C_[a]--0', EdgeInteractionType.consume),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--0', EdgeInteractionType.consume),
                              ('B_[a]--0', 'A_[x]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[x]--0', 'A_[x]_ppi+_B_[a]', EdgeInteractionType.source_state),
                              ('A_[c]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_[a]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{0}', EdgeInteractionType.consume),
                              ('A_[(c)]-{0}', 'C_p+_A_[(c)]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(d)]', 'A_[(d)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(d)]', 'A_[(d)]-{0}', EdgeInteractionType.consume),
                              ('A_[(d)]-{0}', 'C_p+_A_[(d)]', EdgeInteractionType.source_state),

                              ('D_deg_A', 'A_[(c)]-{0}', EdgeInteractionType.degrade),
                              ('D_deg_A', 'A_[x]--0', EdgeInteractionType.degrade),

                              ('D_deg_A', 'A_[(d)]-{0}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[(d)]-{p}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[c]--0', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[c]--C_[a]', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'D_deg_A_ON_A_[c]--C_[a]', EdgeInteractionType.AND),
                              ('A_[c]--C_[a]', 'D_deg_A_ON_A_[c]--C_[a]', EdgeInteractionType.AND),
                              ('D_deg_A_ON_A_[c]--C_[a]', 'C_[a]--0', EdgeInteractionType.produce),
                              ('A_[x]--B_[a]', 'comp', EdgeInteractionType.OR),
                              ('A_[(c)]-{p}', 'comp', EdgeInteractionType.OR),
                              ('comp', 'D_deg_A', EdgeInteractionType.inhibited)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_degradation_boolean_double_negation() -> None:
    test_case = RegulatoryGraphTestCase('''A_[c]_ppi+_C_[a]
                                C_p+_A_[(c)]
                                C_p+_A_[(d)]
                                D_deg_A; x <comp>
                                <comp>; OR A_[c]--C_[a]; OR <NOT>
                                <NOT>; NOT A_[(c)]-{p}''',
                                        ['A_[c]_ppi+_C_[a]', 'C_p+_A_[(c)]', 'C_p+_A_[(d)]', 'D_deg_A'],
                                        ['A_[c]--C_[a]', 'A_[(c)]-{p}', 'C_[a]--0', 'A_[(d)]-{0}', 'A_[(d)]-{p}', 'A_[c]--0',
                              'A_[(c)]-{0}'],
                                        [('comp', 'comp', 'OR'), ('NOT', 'NOT', 'NOT')],
                                        [
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--C_[a]', EdgeInteractionType.produce),
                              ('A_[c]_ppi+_C_[a]', 'C_[a]--0', EdgeInteractionType.consume),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--0', EdgeInteractionType.consume),
                              ('A_[c]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_[a]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{0}', EdgeInteractionType.consume),
                              ('A_[(c)]-{0}', 'C_p+_A_[(c)]', EdgeInteractionType.source_state),
                              ('C_p+_A_[(d)]', 'A_[(d)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(d)]', 'A_[(d)]-{0}', EdgeInteractionType.consume),
                              ('A_[(d)]-{0}', 'C_p+_A_[(d)]', EdgeInteractionType.source_state),

                              ('D_deg_A', 'A_[(c)]-{p}', EdgeInteractionType.degrade),

                              ('D_deg_A', 'A_[c]--0', EdgeInteractionType.degrade),

                              ('D_deg_A', 'A_[(d)]-{0}', EdgeInteractionType.maybe_degraded),
                              ('D_deg_A', 'A_[(d)]-{p}', EdgeInteractionType.maybe_degraded),
                              ('A_[c]--C_[a]', 'comp', EdgeInteractionType.OR),
                              ('NOT', 'comp', EdgeInteractionType.OR),
                              ('A_[(c)]-{p}', 'NOT', EdgeInteractionType.NOT),
                              ('comp', 'D_deg_A', EdgeInteractionType.inhibited)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_degradation_homodimer() -> None:
    """
    Testing degradation of a homodimer

    Note:
        The homodimer should be degraded completely without any production. Only one Boolean complex, producing B_[a]--0,
        should be created.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_[a1]_ppi+_A_[a2]
                                           A_[b]_ppi+_B_[a]
                                           C_deg_A''',
                                        ['A_[a1]_ppi+_A_[a2]', 'A_[b]_ppi+_B_[a]', 'C_deg_A'],
                                        ['A_[a1]--A_[a2]', 'A_[a1]--0', 'A_[a2]--0', 'A_[b]--B_[a]', 'B_[a]--0', 'A_[b]--0'],
                                        [('C_deg_A_ON_A_[b]--B_[a]', " ", 'AND')],
                                        [
                                            ('A_[a1]_ppi+_A_[a2]', 'A_[a1]--A_[a2]', EdgeInteractionType.produce),
                                            ('A_[a1]_ppi+_A_[a2]', 'A_[a1]--0', EdgeInteractionType.consume),
                                            ('A_[a1]_ppi+_A_[a2]', 'A_[a2]--0', EdgeInteractionType.consume),
                                            ('A_[a1]--0', 'A_[a1]_ppi+_A_[a2]', EdgeInteractionType.source_state),
                                            ('A_[a2]--0', 'A_[a1]_ppi+_A_[a2]', EdgeInteractionType.source_state),
                                            ('A_[b]_ppi+_B_[a]', 'A_[b]--B_[a]', EdgeInteractionType.produce),
                                            ('A_[b]_ppi+_B_[a]', 'B_[a]--0', EdgeInteractionType.consume),
                                            ('A_[b]_ppi+_B_[a]', 'A_[b]--0', EdgeInteractionType.consume),
                                            ('B_[a]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                                            ('A_[b]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),

                                            ('C_deg_A', 'A_[b]--B_[a]', EdgeInteractionType.degrade),
                                            ('C_deg_A', 'A_[b]--0', EdgeInteractionType.degrade),
                                            ('C_deg_A', 'A_[a1]--0', EdgeInteractionType.degrade),
                                            ('C_deg_A', 'A_[a2]--0', EdgeInteractionType.degrade),
                                            ('C_deg_A', 'A_[a1]--A_[a2]', EdgeInteractionType.degrade),

                                            ('C_deg_A', 'C_deg_A_ON_A_[b]--B_[a]', EdgeInteractionType.AND),
                                            ('A_[b]--B_[a]', 'C_deg_A_ON_A_[b]--B_[a]', EdgeInteractionType.AND),
                                            ('C_deg_A_ON_A_[b]--B_[a]', 'B_[a]--0', EdgeInteractionType.produce)])
    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_degradation_required_homodimer() -> None:
    """
    Testing degradation of homodimer

    Note:
        The homodimer should be degraded completely without any production. Only B_[a]--0 should be produced by
        the degradation.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_[a1]_ppi+_A_[a2]
                                           A_[b]_ppi+_B_[a]
                                           C_deg_A; ! A@0_[a1]--A@1_[a2]''',
                                        ['A_[a1]_ppi+_A_[a2]', 'A_[b]_ppi+_B_[a]', 'C_deg_A'],
                                        ['A_[a1]--A_[a2]', 'A_[a1]--0', 'A_[a2]--0', 'A_[b]--B_[a]', 'B_[a]--0', 'A_[b]--0'],
                                        [('C_deg_A_ON_A_[b]--B_[a]', " ", 'AND')],
                                        [
                                            ('A_[a1]_ppi+_A_[a2]', 'A_[a1]--A_[a2]', EdgeInteractionType.produce),
                                            ('A_[a1]--A_[a2]', 'C_deg_A', EdgeInteractionType.required),
                                            ('A_[a1]_ppi+_A_[a2]', 'A_[a1]--0', EdgeInteractionType.consume),
                                            ('A_[a1]_ppi+_A_[a2]', 'A_[a2]--0', EdgeInteractionType.consume),
                                            ('A_[a1]--0', 'A_[a1]_ppi+_A_[a2]', EdgeInteractionType.source_state),
                                            ('A_[a2]--0', 'A_[a1]_ppi+_A_[a2]', EdgeInteractionType.source_state),
                                            ('A_[b]_ppi+_B_[a]', 'A_[b]--B_[a]', EdgeInteractionType.produce),
                                            ('A_[b]_ppi+_B_[a]', 'B_[a]--0', EdgeInteractionType.consume),
                                            ('A_[b]_ppi+_B_[a]', 'A_[b]--0', EdgeInteractionType.consume),
                                            ('B_[a]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                                            ('A_[b]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),

                                            ('C_deg_A', 'A_[b]--B_[a]', EdgeInteractionType.maybe_degraded),
                                            ('C_deg_A', 'A_[b]--0', EdgeInteractionType.maybe_degraded),
                                            ('C_deg_A', 'A_[a1]--A_[a2]', EdgeInteractionType.degrade),

                                            ('C_deg_A', 'C_deg_A_ON_A_[b]--B_[a]', EdgeInteractionType.AND),
                                            ('A_[b]--B_[a]', 'C_deg_A_ON_A_[b]--B_[a]', EdgeInteractionType.AND),
                                            ('C_deg_A_ON_A_[b]--B_[a]', 'B_[a]--0', EdgeInteractionType.produce)])
    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_degradation_maybe_homodimer() -> None:
    """
    Testing degradation of homodimer

    Note:
        The homodimer should maybe degraded completely without any production. B_[a]--0 should be produced by the
        degradation.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_[a1]_ppi+_A_[a2]
                                           A_[b]_ppi+_B_[a]
                                           C_deg_A; ! A_[b]--B_[a]''',
                                            ['A_[a1]_ppi+_A_[a2]', 'A_[b]_ppi+_B_[a]', 'C_deg_A'],
                                            ['A_[a1]--A_[a2]', 'A_[a1]--0', 'A_[a2]--0', 'A_[b]--B_[a]', 'B_[a]--0', 'A_[b]--0'],
                                            [('C_deg_A_ON_A_[b]--B_[a]', " ", 'AND')],
                                            [
                                                ('A_[a1]_ppi+_A_[a2]', 'A_[a1]--A_[a2]', EdgeInteractionType.produce),

                                                ('A_[a1]_ppi+_A_[a2]', 'A_[a1]--0', EdgeInteractionType.consume),
                                                ('A_[a1]_ppi+_A_[a2]', 'A_[a2]--0', EdgeInteractionType.consume),
                                                ('A_[a1]--0', 'A_[a1]_ppi+_A_[a2]', EdgeInteractionType.source_state),
                                                ('A_[a2]--0', 'A_[a1]_ppi+_A_[a2]', EdgeInteractionType.source_state),
                                                ('A_[b]_ppi+_B_[a]', 'A_[b]--B_[a]', EdgeInteractionType.produce),
                                                ('A_[b]_ppi+_B_[a]', 'B_[a]--0', EdgeInteractionType.consume),
                                                ('A_[b]_ppi+_B_[a]', 'A_[b]--0', EdgeInteractionType.consume),
                                                ('B_[a]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                                                ('A_[b]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),

                                                ('C_deg_A', 'A_[b]--B_[a]', EdgeInteractionType.degrade),
                                                ('C_deg_A', 'A_[a1]--0', EdgeInteractionType.maybe_degraded),
                                                ('C_deg_A', 'A_[a2]--0', EdgeInteractionType.maybe_degraded),
                                                ('C_deg_A', 'A_[a1]--A_[a2]', EdgeInteractionType.maybe_degraded),

                                                ('A_[b]--B_[a]', 'C_deg_A', EdgeInteractionType.required),

                                                ('C_deg_A', 'C_deg_A_ON_A_[b]--B_[a]', EdgeInteractionType.AND),
                                                ('A_[b]--B_[a]', 'C_deg_A_ON_A_[b]--B_[a]', EdgeInteractionType.AND),
                                                ('C_deg_A_ON_A_[b]--B_[a]', 'B_[a]--0', EdgeInteractionType.produce)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_degradation_homodimer_inhibited_heterodimer() -> None:
    """
    Testing possible degradation of homodimer if the degradation is inhibited by a heterodimer.

    Note:
        The homodimer should maybe degraded completely without any production.

    Returns:
        None

    Raises:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''A_[a1]_ppi+_A_[a2]
                                           A_[b]_ppi+_B_[a]
                                           C_deg_A; x A_[b]--B_[a]''',
                                            ['A_[a1]_ppi+_A_[a2]', 'A_[b]_ppi+_B_[a]', 'C_deg_A'],
                                            ['A_[a1]--A_[a2]', 'A_[a1]--0', 'A_[a2]--0', 'A_[b]--B_[a]', 'B_[a]--0', 'A_[b]--0'],
                                            [],
                                            [
                                                ('A_[a1]_ppi+_A_[a2]', 'A_[a1]--A_[a2]', EdgeInteractionType.produce),

                                                ('A_[a1]_ppi+_A_[a2]', 'A_[a1]--0', EdgeInteractionType.consume),
                                                ('A_[a1]_ppi+_A_[a2]', 'A_[a2]--0', EdgeInteractionType.consume),
                                                ('A_[a1]--0', 'A_[a1]_ppi+_A_[a2]', EdgeInteractionType.source_state),
                                                ('A_[a2]--0', 'A_[a1]_ppi+_A_[a2]', EdgeInteractionType.source_state),
                                                ('A_[b]_ppi+_B_[a]', 'A_[b]--B_[a]', EdgeInteractionType.produce),
                                                ('A_[b]_ppi+_B_[a]', 'B_[a]--0', EdgeInteractionType.consume),
                                                ('A_[b]_ppi+_B_[a]', 'A_[b]--0', EdgeInteractionType.consume),
                                                ('B_[a]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),
                                                ('A_[b]--0', 'A_[b]_ppi+_B_[a]', EdgeInteractionType.source_state),

                                                ('C_deg_A', 'A_[b]--0', EdgeInteractionType.degrade),
                                                ('C_deg_A', 'A_[a1]--0', EdgeInteractionType.maybe_degraded),
                                                ('C_deg_A', 'A_[a2]--0', EdgeInteractionType.maybe_degraded),
                                                ('C_deg_A', 'A_[a1]--A_[a2]', EdgeInteractionType.maybe_degraded),

                                                ('A_[b]--B_[a]', 'C_deg_A', EdgeInteractionType.inhibited)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_synthesis_of_neutral_state() -> None:
    """
    Testing regulatory graph for synthesis reaction.

    Note:
        3 reaction, 0 Contingencies.

    Returns:
        AssertionError: If generated graph differs from expected graph.

    """
    test_case = RegulatoryGraphTestCase('''C_syn_A
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
                              ('C_syn_A', 'A_[b]--0', EdgeInteractionType.synthesis),
                              ('C_syn_A', 'A_[(c)]-{0}', EdgeInteractionType.synthesis),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{0}', EdgeInteractionType.consume),
                              ('A_[(c)]-{0}', 'C_p+_A_[(c)]', EdgeInteractionType.source_state)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_synthesis_of_component_without_states() -> None:
    test_case = RegulatoryGraphTestCase('''D_trsc_A''',
                                        ['D_trsc_AGene'],
                                        ['AmRNA#component'],
                                        [],
                                        [('D_trsc_AGene', 'AmRNA', EdgeInteractionType.synthesis)])
    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_trsc_trsl_of_neutral_state() -> None:
    test_case = RegulatoryGraphTestCase('''D_trsc_A
                                E_trsl_A
                                C_p+_A_[(c)]
                                A_[c]_ppi+_C_[a]''',
                                        ['D_trsc_AGene', 'E_trsl_AmRNA', 'C_p+_A_[(c)]', 'A_[c]_ppi+_C_[a]'],
                                        ['AmRNA#component', 'A_[(c)]-{0}', 'A_[(c)]-{p}', 'A_[c]--C_[a]', 'A_[c]--0', 'C_[a]--0'],
                                        [],
                                        [('D_trsc_AGene', 'AmRNA', EdgeInteractionType.synthesis),
                              ('E_trsl_AmRNA', 'A_[(c)]-{0}', EdgeInteractionType.synthesis),
                              ('E_trsl_AmRNA', 'A_[c]--0', EdgeInteractionType.synthesis),
                              ('AmRNA', 'E_trsl_AmRNA', EdgeInteractionType.input_state),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{p}', EdgeInteractionType.produce),
                              ('C_p+_A_[(c)]', 'A_[(c)]-{0}', EdgeInteractionType.consume),
                              ('A_[(c)]-{0}', 'C_p+_A_[(c)]', EdgeInteractionType.source_state),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--C_[a]', EdgeInteractionType.produce),
                              ('A_[c]_ppi+_C_[a]', 'C_[a]--0', EdgeInteractionType.consume),
                              ('A_[c]_ppi+_C_[a]', 'A_[c]--0', EdgeInteractionType.consume),
                              ('A_[c]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state),
                              ('C_[a]--0', 'A_[c]_ppi+_C_[a]', EdgeInteractionType.source_state)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_trsc_trsl_of_component_without_states() -> None:
    test_case = RegulatoryGraphTestCase('''D_trsc_A
                                E_trsl_A''',
                                        ['D_trsc_AGene', 'E_trsl_AmRNA'],
                                        ['AmRNA#component', 'A#comp'],
                                        [],
                                        [('D_trsc_AGene', 'AmRNA', EdgeInteractionType.synthesis),
                              ('E_trsl_AmRNA', 'A', EdgeInteractionType.synthesis),
                              ('AmRNA', 'E_trsl_AmRNA', EdgeInteractionType.input_state)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def test_synthesis_degradation_of_component_without_states() -> None:
    test_case = RegulatoryGraphTestCase('''D_syn_A
                                E_deg_A''',
                                        ['D_syn_A', 'E_deg_A'],
                                        ['A#component'],
                                        [],
                                        [('D_syn_A', 'A', EdgeInteractionType.synthesis),
                              ('E_deg_A', 'A', EdgeInteractionType.degrade),
                              ('A', 'E_deg_A', EdgeInteractionType.source_state)])

    assert _is_graph_test_case_correct(_create_regulatory_graph(test_case.quick_string), test_case)


def _get_state_nodes(test_case: RegulatoryGraphTestCase, expected_graph: DiGraph) -> DiGraph:
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
        elif '#comp' in node:
            expected_graph.add_node(node.split('#comp')[0], type=NodeType.component.value, label=node.split('#comp')[0])
        else:
            expected_graph.add_node(node, type=NodeType.state.value, label=node)
    return expected_graph


def _get_reaction_nodes(test_case: RegulatoryGraphTestCase, expected_graph: DiGraph) -> DiGraph:
    for reaction_id in test_case.reaction_node_strings:

        if '#out' in reaction_id:
            expected_graph.add_node(reaction_id.split('#out')[0], type=NodeType.output.value)
        else:
            reaction_label = reaction_id.split('#')[0]
            expected_graph.add_node(reaction_id, type=NodeType.reaction.value, label=reaction_label)
    return expected_graph


def _get_boolean_complex_state_nodes(test_case: RegulatoryGraphTestCase, expected_graph: DiGraph) -> DiGraph:
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
        elif 'boolean' == type:
            expected_graph.add_node(id, type=NodeType.boolean.value, label=label)
        else:
            raise AssertionError
    return expected_graph


def _is_graph_test_case_correct(actual_graph: DiGraph, test_case: RegulatoryGraphTestCase) -> bool:
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
    for edge in expected_graph.edges:
        assert expected_graph.edges[edge] == actual_graph.edges[edge]
    return expected_graph.node == actual_graph.node


def _create_regulatory_graph(quick_string: str) -> DiGraph:
    """
    Creating a regulatory graph.

    Args:
        quick_string: A rxncon system in quick format.

    Returns:
        A regulatory graph.

    """
    actual_system = qui.Quick(quick_string)
    reg_system = SpeciesReactionGraph(actual_system.rxncon_system)
    return reg_system.to_graph()


