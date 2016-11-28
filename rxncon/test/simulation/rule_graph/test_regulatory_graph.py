import pytest
import networkx as nex
from collections import namedtuple
import rxncon.input.quick.quick as qui
import rxncon.simulation.rule_graph.regulatory_graph as reg
from rxncon.simulation.rule_graph.graphML import XGMML

RuleTestCase = namedtuple('RuleTestCase', ['quick_string', 'reaction_node_strings', 'state_node_strings',
                                           'boolean_state_node_strings', 'edge_tuples'])


def test_regulatory_graph_generation(the_cases_regulatory_graph):
    for test_case in the_cases_regulatory_graph:
        actual_system = qui.Quick(test_case.quick_string)
        reg_system = reg.RegulatoryGraph(actual_system.rxncon_system)
        actual_graph = reg_system.to_graph()
        assert is_graph_test_case_correct(actual_graph, test_case)


@pytest.fixture
def the_cases_regulatory_graph(case_and_expected_regulatory_graph):
    return case_and_expected_regulatory_graph


@pytest.fixture
def case_and_expected_regulatory_graph():
    return [
        RuleTestCase('''A_[b]_ppi+_B_[a]; ! A-{p}
                        C_p+_A_[(c)]''',
                     ['A_[b]_ppi+_B_[a]', 'C_p+_A_[(c)]'],
                     ['A_[b]--B_[a]', 'A_[(c)]-{p}', 'B_[a]--0', 'A_[b]--0', 'A_[(c)]-{0}'],
                     [],
                     [('A_[b]_ppi+_B_[a]', 'A_[b]--B_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[b]_ppi+_B_[a]', 'B_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[b]_ppi+_B_[a]', 'A_[b]--0', reg.EdgeInteractionType.consume.value),
                      ('B_[a]--0', 'A_[b]_ppi+_B_[a]', reg.EdgeInteractionType.source_state.value),
                      ('A_[b]--0', 'A_[b]_ppi+_B_[a]', reg.EdgeInteractionType.source_state.value),
                      ('C_p+_A_[(c)]', 'A_[(c)]-{p}', reg.EdgeInteractionType.produce.value),
                      ('C_p+_A_[(c)]', 'A_[(c)]-{0}', reg.EdgeInteractionType.consume.value),
                      ('A_[(c)]-{0}', 'C_p+_A_[(c)]', reg.EdgeInteractionType.source_state.value),
                      ('A_[(c)]-{p}', 'A_[b]_ppi+_B_[a]', reg.EdgeInteractionType.required.value)]),

        RuleTestCase('''C_syn_A
                        A_[b]_ppi+_B_[a]
                        C_p+_A_[(c)]''',
                     ['A_[b]_ppi+_B_[a]', 'C_p+_A_[(c)]', 'C_syn_A'],
                     ['A_[b]--B_[a]', 'A_[(c)]-{p}', 'B_[a]--0', 'A_[b]--0', 'A_[(c)]-{0}'],
                     [],
                     [('A_[b]_ppi+_B_[a]', 'A_[b]--B_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[b]_ppi+_B_[a]', 'B_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[b]_ppi+_B_[a]', 'A_[b]--0', reg.EdgeInteractionType.consume.value),
                      ('B_[a]--0', 'A_[b]_ppi+_B_[a]', reg.EdgeInteractionType.source_state.value),
                      ('A_[b]--0', 'A_[b]_ppi+_B_[a]', reg.EdgeInteractionType.source_state.value),
                      ('C_syn_A', 'A_[b]--0', reg.EdgeInteractionType.produce.value),
                      ('C_syn_A', 'A_[(c)]-{0}', reg.EdgeInteractionType.produce.value),
                      ('C_p+_A_[(c)]', 'A_[(c)]-{p}', reg.EdgeInteractionType.produce.value),
                      ('C_p+_A_[(c)]', 'A_[(c)]-{0}', reg.EdgeInteractionType.consume.value),
                      ('A_[(c)]-{0}', 'C_p+_A_[(c)]', reg.EdgeInteractionType.source_state.value)]),

        RuleTestCase('''A_[b]_ppi+_B_[a]; ! <comp>; ! C-{p}
                        <comp>; AND A-{p}; AND A--C
                        A_[c]_ppi+_C_[a]
                        C_p+_A_[(c)]
                        D_p+_C_[(d)]''',
                     ['A_[b]_ppi+_B_[a]', 'A_[c]_ppi+_C_[a]', 'C_p+_A_[(c)]', 'D_p+_C_[(d)]'],
                     ['A_[b]--B_[a]', 'A_[b]--0', 'B_[a]--0', 'A_[c]--C_[a]','A_[c]--0', 'C_[a]--0', 'A_[(c)]-{p}', 'A_[(c)]-{0}' , 'C_[(d)]-{p}', 'C_[(d)]-{0}'],
                     ['comp#AND'],
                     [('A_[b]_ppi+_B_[a]', 'A_[b]--B_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[b]_ppi+_B_[a]', 'A_[b]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[b]_ppi+_B_[a]', 'B_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[b]--0', 'A_[b]_ppi+_B_[a]', reg.EdgeInteractionType.source_state.value),
                      ('B_[a]--0', 'A_[b]_ppi+_B_[a]', reg.EdgeInteractionType.source_state.value),
                      ('A_[c]_ppi+_C_[a]', 'A_[c]--C_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[c]_ppi+_C_[a]', 'A_[c]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[c]_ppi+_C_[a]', 'C_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[c]--0', 'A_[c]_ppi+_C_[a]', reg.EdgeInteractionType.source_state.value),
                      ('C_[a]--0', 'A_[c]_ppi+_C_[a]', reg.EdgeInteractionType.source_state.value),
                      ('C_p+_A_[(c)]', 'A_[(c)]-{p}', reg.EdgeInteractionType.produce.value),
                      ('C_p+_A_[(c)]', 'A_[(c)]-{0}', reg.EdgeInteractionType.consume.value),
                      ('A_[(c)]-{0}', 'C_p+_A_[(c)]', reg.EdgeInteractionType.source_state.value),
                      ('D_p+_C_[(d)]', 'C_[(d)]-{p}', reg.EdgeInteractionType.produce.value),
                      ('D_p+_C_[(d)]', 'C_[(d)]-{0}', reg.EdgeInteractionType.consume.value),
                      ('C_[(d)]-{0}', 'D_p+_C_[(d)]', reg.EdgeInteractionType.source_state.value),
                      ('C_[(d)]-{p}', 'A_[b]_ppi+_B_[a]', reg.EdgeInteractionType.required.value),
                      ('comp', 'A_[b]_ppi+_B_[a]', reg.EdgeInteractionType.required.value),
                      ('A_[(c)]-{p}', 'comp', reg.EdgeInteractionType.AND.value),
                      ('A_[c]--C_[a]', 'comp', reg.EdgeInteractionType.AND.value)]),
        # # #
        RuleTestCase('''A_[b]_ppi+_B_[a]; ! <comp>
                        <comp>; AND <comp1>; AND <comp2>
                        <comp1>; OR <comp3>; OR A--C
                        <comp2>; AND A_[d]--D_[a]; AND A--E
                        <comp3>; AND A_[f]--F_[a]; AND A--G
                        A_[c]_ppi+_C_[a]
                        A_[d]_ppi+_D_[a]
                        A_[e]_ppi+_E_[a]
                        A_[f]_ppi+_F_[a]
                        A_[g]_ppi+_G_[a]''',
                     ['A_[b]_ppi+_B_[a]', 'A_[c]_ppi+_C_[a]', 'A_[d]_ppi+_D_[a]', 'A_[e]_ppi+_E_[a]', 'A_[f]_ppi+_F_[a]', 'A_[g]_ppi+_G_[a]' ],
                     ['A_[b]--B_[a]', 'A_[b]--0', 'B_[a]--0', 'A_[c]--C_[a]', 'A_[c]--0', 'C_[a]--0',  'A_[d]--D_[a]', 'A_[d]--0', 'D_[a]--0', 'A_[e]--E_[a]', 'A_[e]--0', 'E_[a]--0', 'A_[f]--F_[a]', 'A_[f]--0', 'F_[a]--0',
                      'A_[g]--G_[a]', 'A_[g]--0', 'G_[a]--0'],
                     ['comp#AND', 'comp1#OR', 'comp2#AND', 'comp3#AND'],
                     [('A_[b]_ppi+_B_[a]', 'A_[b]--B_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[b]_ppi+_B_[a]', 'A_[b]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[b]_ppi+_B_[a]', 'B_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[b]--0', 'A_[b]_ppi+_B_[a]', reg.EdgeInteractionType.source_state.value),
                      ('B_[a]--0', 'A_[b]_ppi+_B_[a]', reg.EdgeInteractionType.source_state.value),
                      ('A_[c]_ppi+_C_[a]', 'A_[c]--C_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[c]_ppi+_C_[a]', 'A_[c]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[c]_ppi+_C_[a]', 'C_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[c]--0', 'A_[c]_ppi+_C_[a]', reg.EdgeInteractionType.source_state.value),
                      ('C_[a]--0', 'A_[c]_ppi+_C_[a]', reg.EdgeInteractionType.source_state.value),
                      ('A_[d]_ppi+_D_[a]', 'A_[d]--D_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[d]_ppi+_D_[a]', 'A_[d]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[d]_ppi+_D_[a]', 'D_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[d]--0', 'A_[d]_ppi+_D_[a]', reg.EdgeInteractionType.source_state.value),
                      ('D_[a]--0', 'A_[d]_ppi+_D_[a]', reg.EdgeInteractionType.source_state.value),
                      ('A_[e]_ppi+_E_[a]', 'A_[e]--E_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[e]_ppi+_E_[a]', 'A_[e]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[e]_ppi+_E_[a]', 'E_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[e]--0', 'A_[e]_ppi+_E_[a]', reg.EdgeInteractionType.source_state.value),
                      ('E_[a]--0', 'A_[e]_ppi+_E_[a]', reg.EdgeInteractionType.source_state.value),
                      ('A_[f]_ppi+_F_[a]', 'A_[f]--F_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[f]_ppi+_F_[a]', 'A_[f]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[f]_ppi+_F_[a]', 'F_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[f]--0', 'A_[f]_ppi+_F_[a]', reg.EdgeInteractionType.source_state.value),
                      ('F_[a]--0', 'A_[f]_ppi+_F_[a]', reg.EdgeInteractionType.source_state.value),
                      ('A_[g]_ppi+_G_[a]', 'A_[g]--G_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[g]_ppi+_G_[a]', 'A_[g]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[g]_ppi+_G_[a]', 'G_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[g]--0', 'A_[g]_ppi+_G_[a]', reg.EdgeInteractionType.source_state.value),
                      ('G_[a]--0', 'A_[g]_ppi+_G_[a]', reg.EdgeInteractionType.source_state.value),
                      ('comp', 'A_[b]_ppi+_B_[a]', reg.EdgeInteractionType.required.value),
                      ('comp1', 'comp', reg.EdgeInteractionType.AND.value),
                      ('comp2', 'comp', reg.EdgeInteractionType.AND.value),
                      ('comp3', 'comp1', reg.EdgeInteractionType.OR.value),
                      ('A_[c]--C_[a]', 'comp1', reg.EdgeInteractionType.OR.value),
                      ('A_[d]--D_[a]', 'comp2', reg.EdgeInteractionType.AND.value),
                      ('A_[e]--E_[a]', 'comp2', reg.EdgeInteractionType.AND.value),
                      ('A_[f]--F_[a]', 'comp3', reg.EdgeInteractionType.AND.value),
                      ('A_[g]--G_[a]', 'comp3', reg.EdgeInteractionType.AND.value),]
                     ),

        RuleTestCase('''A_[b]_ppi+_B_[a]; ! <comp>
                        <comp>; AND <comp1>; AND <Notcomp2>
                        <comp1>; OR <comp3>; OR A--C
                        <Notcomp2>; NOT A--D
                        <comp3>; AND A--F; AND A--G
                        A_[c]_ppi+_C_[a]
                        A_[d]_ppi+_D_[a]
                        A_[f]_ppi+_F_[a]
                        A_[g]_ppi+_G_[a]''',
                     ['A_[b]_ppi+_B_[a]', 'A_[c]_ppi+_C_[a]', 'A_[d]_ppi+_D_[a]', 'A_[f]_ppi+_F_[a]', 'A_[g]_ppi+_G_[a]' ],
                     ['A_[b]--B_[a]', 'A_[b]--0', 'B_[a]--0', 'A_[c]--C_[a]', 'A_[c]--0', 'C_[a]--0', 'A_[d]--D_[a]',
                      'A_[d]--0', 'D_[a]--0', 'A_[f]--F_[a]', 'A_[f]--0', 'F_[a]--0', 'A_[g]--G_[a]', 'A_[g]--0',
                      'G_[a]--0'],
                     ['comp#AND', 'comp1#OR', 'Notcomp2#NOT', 'comp3#AND'],
                     [('A_[b]_ppi+_B_[a]', 'A_[b]--B_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[b]_ppi+_B_[a]', 'A_[b]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[b]_ppi+_B_[a]', 'B_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[b]--0', 'A_[b]_ppi+_B_[a]', reg.EdgeInteractionType.source_state.value),
                      ('B_[a]--0', 'A_[b]_ppi+_B_[a]', reg.EdgeInteractionType.source_state.value),
                      ('A_[c]_ppi+_C_[a]', 'A_[c]--C_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[c]_ppi+_C_[a]', 'A_[c]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[c]_ppi+_C_[a]', 'C_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[c]--0', 'A_[c]_ppi+_C_[a]', reg.EdgeInteractionType.source_state.value),
                      ('C_[a]--0', 'A_[c]_ppi+_C_[a]', reg.EdgeInteractionType.source_state.value),
                      ('A_[d]_ppi+_D_[a]', 'A_[d]--D_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[d]_ppi+_D_[a]', 'A_[d]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[d]_ppi+_D_[a]', 'D_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[d]--0', 'A_[d]_ppi+_D_[a]', reg.EdgeInteractionType.source_state.value),
                      ('D_[a]--0', 'A_[d]_ppi+_D_[a]', reg.EdgeInteractionType.source_state.value),
                      ('A_[f]_ppi+_F_[a]', 'A_[f]--F_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[f]_ppi+_F_[a]', 'A_[f]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[f]_ppi+_F_[a]', 'F_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[f]--0', 'A_[f]_ppi+_F_[a]', reg.EdgeInteractionType.source_state.value),
                      ('F_[a]--0', 'A_[f]_ppi+_F_[a]', reg.EdgeInteractionType.source_state.value),
                      ('A_[g]_ppi+_G_[a]', 'A_[g]--G_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[g]_ppi+_G_[a]', 'A_[g]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[g]_ppi+_G_[a]', 'G_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[g]--0', 'A_[g]_ppi+_G_[a]', reg.EdgeInteractionType.source_state.value),
                      ('G_[a]--0', 'A_[g]_ppi+_G_[a]', reg.EdgeInteractionType.source_state.value),
                      ('comp', 'A_[b]_ppi+_B_[a]', reg.EdgeInteractionType.required.value),
                      ('comp1', 'comp', reg.EdgeInteractionType.AND.value),
                      ('Notcomp2', 'comp', reg.EdgeInteractionType.AND.value),
                      ('comp3', 'comp1', reg.EdgeInteractionType.OR.value),
                      ('A_[c]--C_[a]', 'comp1', reg.EdgeInteractionType.OR.value),
                      ('A_[d]--D_[a]', 'Notcomp2', reg.EdgeInteractionType.NOT.value),
                      ('A_[f]--F_[a]', 'comp3', reg.EdgeInteractionType.AND.value),
                      ('A_[g]--G_[a]', 'comp3', reg.EdgeInteractionType.AND.value),]
                     ),
        #
        RuleTestCase('''A_p+_B_[(a)]
                        C_p-_B_[(a)]''',
                     ['A_p+_B_[(a)]', 'C_p-_B_[(a)]'],
                     ['B_[(a)]-{p}', 'B_[(a)]-{0}'],
                     [],
                     [('A_p+_B_[(a)]', 'B_[(a)]-{p}', reg.EdgeInteractionType.produce.value),
                      ('A_p+_B_[(a)]', 'B_[(a)]-{0}', reg.EdgeInteractionType.consume.value),
                      ('B_[(a)]-{0}', 'A_p+_B_[(a)]', reg.EdgeInteractionType.source_state.value),
                      ('C_p-_B_[(a)]', 'B_[(a)]-{p}', reg.EdgeInteractionType.consume.value),
                      ('B_[(a)]-{p}', 'C_p-_B_[(a)]', reg.EdgeInteractionType.source_state.value),
                      ('C_p-_B_[(a)]', 'B_[(a)]-{0}', reg.EdgeInteractionType.produce.value)]),

        #
        RuleTestCase('''A_[B]_ppi+_B_[A]; ! [Input]''',
                     ['A_[B]_ppi+_B_[A]'],
                     ['A_[B]--B_[A]', 'A_[B]--0', 'B_[A]--0', '[Input]#in'],
                     [],
                     [('A_[B]_ppi+_B_[A]', 'A_[B]--B_[A]', reg.EdgeInteractionType.produce.value),
                      ('A_[B]_ppi+_B_[A]', 'A_[B]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[B]_ppi+_B_[A]', 'B_[A]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[B]--0', 'A_[B]_ppi+_B_[A]', reg.EdgeInteractionType.source_state.value),
                      ('B_[A]--0', 'A_[B]_ppi+_B_[A]', reg.EdgeInteractionType.source_state.value),
                      ('[Input]', 'A_[B]_ppi+_B_[A]', reg.EdgeInteractionType.required.value)]),

        RuleTestCase('''A_ppi+_B; ! <comp>
                        <comp>; AND <comp1>; AND [Input]
                        <comp1>; OR A--D; OR A--C
                        A_ppi+_D
                        A_ppi+_C''',
                     ['A_[B]_ppi+_B_[A]', 'A_[C]_ppi+_C_[A]', 'A_[D]_ppi+_D_[A]'],
                     ['A_[B]--B_[A]','A_[B]--0', 'B_[A]--0', 'A_[C]--C_[A]', 'A_[C]--0', 'C_[A]--0', 'A_[D]--D_[A]',
                      'A_[D]--0', 'D_[A]--0', '[Input]#in'],
                     ['comp#AND', 'comp1#OR'],
                     [('A_[B]_ppi+_B_[A]', 'A_[B]--B_[A]', reg.EdgeInteractionType.produce.value),
                      ('A_[B]_ppi+_B_[A]', 'A_[B]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[B]_ppi+_B_[A]', 'B_[A]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[B]--0', 'A_[B]_ppi+_B_[A]', reg.EdgeInteractionType.source_state.value),
                      ('B_[A]--0', 'A_[B]_ppi+_B_[A]', reg.EdgeInteractionType.source_state.value),
                      ('A_[C]_ppi+_C_[A]', 'A_[C]--C_[A]', reg.EdgeInteractionType.produce.value),
                      ('A_[C]_ppi+_C_[A]', 'A_[C]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[C]_ppi+_C_[A]', 'C_[A]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[C]--0', 'A_[C]_ppi+_C_[A]', reg.EdgeInteractionType.source_state.value),
                      ('C_[A]--0', 'A_[C]_ppi+_C_[A]', reg.EdgeInteractionType.source_state.value),
                      ('A_[D]_ppi+_D_[A]', 'A_[D]--D_[A]', reg.EdgeInteractionType.produce.value),
                      ('A_[D]_ppi+_D_[A]', 'A_[D]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[D]_ppi+_D_[A]', 'D_[A]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[D]--0', 'A_[D]_ppi+_D_[A]', reg.EdgeInteractionType.source_state.value),
                      ('D_[A]--0', 'A_[D]_ppi+_D_[A]', reg.EdgeInteractionType.source_state.value),
                      ('comp', 'A_[B]_ppi+_B_[A]', reg.EdgeInteractionType.required.value),
                      ('[Input]', 'comp', reg.EdgeInteractionType.AND.value),
                      ('comp1', 'comp', reg.EdgeInteractionType.AND.value),
                      ('A_[D]--D_[A]', 'comp1', reg.EdgeInteractionType.OR.value),
                      ('A_[C]--C_[A]', 'comp1', reg.EdgeInteractionType.OR.value)]),

    ]

def get_state_nodes(test_case, expected_graph):
    for node in test_case.state_node_strings:
        if '#comp' in node:
            expected_graph.add_node(node.split('#comp')[0], type=reg.NodeType.componentstate.value)
        elif '#in' in node:
            expected_graph.add_node(node.split('#in')[0], type=reg.NodeType.input.value)
        else:
            expected_graph.add_node(node, type=reg.NodeType.state.value)
    return expected_graph


def get_boolean_complex_state_nodes(test_case, expected_graph):
    for node in test_case.boolean_state_node_strings:
        if '#AND' in node:
            expected_graph.add_node(node.split('#AND')[0], type=reg.NodeType.AND.value)
        elif '#OR' in node:
            expected_graph.add_node(node.split('#OR')[0], type=reg.NodeType.OR.value)
        elif '#NOT' in node:
            expected_graph.add_node(node.split('#NOT')[0], type=reg.NodeType.NOT.value)
        else:
            raise AssertionError
    return expected_graph


def is_graph_test_case_correct(actual_graph, test_case) -> bool:

    #gml_system = XGMML(actual_graph, "test_graph_NOT")
    #gml_system.to_file('/home/thiemese/data/projects/graph/test_boolean_NOT1.xgmml')
    expected_graph = nex.DiGraph()
    [expected_graph.add_node(node, type=reg.NodeType.reaction.value) if '#out' not in node
     else expected_graph.add_node(node.split('#out')[0], type=reg.NodeType.output.value) for node in test_case.reaction_node_strings]

    expected_graph = get_state_nodes(test_case, expected_graph)
    expected_graph = get_boolean_complex_state_nodes(test_case, expected_graph)

    [expected_graph.add_edge(source, target, interaction=interaction) for source, target, interaction in test_case.edge_tuples]

    return expected_graph.node == actual_graph.node and expected_graph.edge == actual_graph.edge
