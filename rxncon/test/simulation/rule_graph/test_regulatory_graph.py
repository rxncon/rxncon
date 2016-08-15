import pytest
import networkx as nex
from collections import namedtuple
import rxncon.input.quick.quick as qui
import rxncon.simulation.rule_graph.regulatory_graph as reg
import rxncon.core.effector as eff
import rxncon.core.state as sta
import rxncon.core.reaction as rxn

RuleTestCase = namedtuple('RuleTestCase', ['quick_string', 'reaction_node_strings', 'state_node_strings',
                                           'boolean_state_node_strings', 'edge_tuples'])


# def test_get_subspecifications_of_state():
#     quick_system = qui.Quick("""D_p+_A
#                                 D_p+_A_[d]
#                                 D_p+_A_[(r)]""")
#
#     state_a_p = sta.state_from_string('A-{P}')
#     reaction_D_p_A = rxn.reaction_from_string('D_p+_A')
#     reaction_D_p_A_d = rxn.reaction_from_string('D_p+_A_[d]')
#     reaction_D_p_A_r = rxn.reaction_from_string('D_p+_A_[(r)]')
#     expected_subspecifications = {eff.StateEffector(reaction_D_p_A.product), eff.StateEffector(reaction_D_p_A_d.product),
#                                   eff.StateEffector(reaction_D_p_A_r.product)}
#     reg_graph = reg.RegulatoryGraph(quick_system.rxncon_system)
#     actual_subspecifications = reg_graph.get_subset_of_state(state_a_p)
#
#     assert actual_subspecifications == expected_subspecifications

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
        RuleTestCase('''A_[b]_ppi_B_[a]; ! A-{p}
                        C_p+_A_[(c)]''',
                     ['A_[b]_ppi_B_[a]', 'C_p+_A_[(c)]'],
                     ['A_[b]--B_[a]', 'A_[(c)]-{p}', 'B_[a]--0', 'A_[b]--0', 'A_[(c)]-{0}'],
                     [],
                     [('A_[b]_ppi_B_[a]', 'A_[b]--B_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[b]_ppi_B_[a]', 'B_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[b]_ppi_B_[a]', 'A_[b]--0', reg.EdgeInteractionType.consume.value),
                      ('C_p+_A_[(c)]', 'A_[(c)]-{p}', reg.EdgeInteractionType.produce.value),
                      ('C_p+_A_[(c)]', 'A_[(c)]-{0}', reg.EdgeInteractionType.consume.value),
                      ('A_[(c)]-{p}', 'A_[b]_ppi_B_[a]', reg.EdgeInteractionType.required.value)]),

        RuleTestCase('''A_[b]_ppi_B_[a]; ! <comp>; ! C-{p}
                        <comp>; AND A-{p}; AND A--C
                        A_[c]_ppi_C_[a]
                        C_p+_A_[(c)]
                        D_p+_C_[(d)]''',
                     ['A_[b]_ppi_B_[a]', 'A_[c]_ppi_C_[a]', 'C_p+_A_[(c)]', 'D_p+_C_[(d)]'],
                     ['A_[b]--B_[a]', 'A_[b]--0', 'B_[a]--0', 'A_[c]--C_[a]','A_[c]--0', 'C_[a]--0', 'A_[(c)]-{p}', 'A_[(c)]-{0}' , 'C_[(d)]-{p}', 'C_[(d)]-{0}'],
                     ['comp#AND'],
                     [('A_[b]_ppi_B_[a]', 'A_[b]--B_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[b]_ppi_B_[a]', 'A_[b]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[b]_ppi_B_[a]', 'B_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[c]_ppi_C_[a]', 'A_[c]--C_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[c]_ppi_C_[a]', 'A_[c]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[c]_ppi_C_[a]', 'C_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('C_p+_A_[(c)]', 'A_[(c)]-{p}', reg.EdgeInteractionType.produce.value),
                      ('C_p+_A_[(c)]', 'A_[(c)]-{0}', reg.EdgeInteractionType.consume.value),
                      ('D_p+_C_[(d)]', 'C_[(d)]-{p}', reg.EdgeInteractionType.produce.value),
                      ('D_p+_C_[(d)]', 'C_[(d)]-{0}', reg.EdgeInteractionType.consume.value),
                      ('C_[(d)]-{p}', 'A_[b]_ppi_B_[a]', reg.EdgeInteractionType.required.value),
                      ('comp', 'A_[b]_ppi_B_[a]', reg.EdgeInteractionType.required.value),
                      ('A_[(c)]-{p}', 'comp', reg.EdgeInteractionType.AND.value),
                      ('A_[c]--C_[a]', 'comp', reg.EdgeInteractionType.AND.value)]),
        #
        RuleTestCase('''A_[b]_ppi_B_[a]; ! <comp>
                        <comp>; AND <comp1>; AND <comp2>
                        <comp1>; OR <comp3>; OR A--C
                        <comp2>; AND A_[d]--D_[a]; AND A--E
                        <comp3>; AND A_[f]--F_[a]; AND A--G
                        A_[c]_ppi_C_[a]
                        A_[d]_ppi_D_[a]
                        A_[e]_ppi_E_[a]
                        A_[f]_ppi_F_[a]
                        A_[g]_ppi_G_[a]''',
                     ['A_[b]_ppi_B_[a]', 'A_[c]_ppi_C_[a]', 'A_[d]_ppi_D_[a]', 'A_[e]_ppi_E_[a]', 'A_[f]_ppi_F_[a]', 'A_[g]_ppi_G_[a]' ],
                     ['A_[b]--B_[a]', 'A_[b]--0', 'B_[a]--0', 'A_[c]--C_[a]', 'A_[c]--0', 'C_[a]--0',  'A_[d]--D_[a]', 'A_[d]--0', 'D_[a]--0', 'A_[e]--E_[a]', 'A_[e]--0', 'E_[a]--0', 'A_[f]--F_[a]', 'A_[f]--0', 'F_[a]--0',
                      'A_[g]--G_[a]', 'A_[g]--0', 'G_[a]--0'],
                     ['comp#AND', 'comp1#OR', 'comp2#AND', 'comp3#AND'],
                     [('A_[b]_ppi_B_[a]', 'A_[b]--B_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[b]_ppi_B_[a]', 'A_[b]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[b]_ppi_B_[a]', 'B_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[c]_ppi_C_[a]', 'A_[c]--C_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[c]_ppi_C_[a]', 'A_[c]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[c]_ppi_C_[a]', 'C_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[d]_ppi_D_[a]', 'A_[d]--D_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[d]_ppi_D_[a]', 'A_[d]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[d]_ppi_D_[a]', 'D_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[e]_ppi_E_[a]', 'A_[e]--E_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[e]_ppi_E_[a]', 'A_[e]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[e]_ppi_E_[a]', 'E_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[f]_ppi_F_[a]', 'A_[f]--F_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[f]_ppi_F_[a]', 'A_[f]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[f]_ppi_F_[a]', 'F_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[g]_ppi_G_[a]', 'A_[g]--G_[a]', reg.EdgeInteractionType.produce.value),
                      ('A_[g]_ppi_G_[a]', 'A_[g]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[g]_ppi_G_[a]', 'G_[a]--0', reg.EdgeInteractionType.consume.value),
                      ('comp', 'A_[b]_ppi_B_[a]', reg.EdgeInteractionType.required.value),
                      ('comp1', 'comp', reg.EdgeInteractionType.AND.value),
                      ('comp2', 'comp', reg.EdgeInteractionType.AND.value),
                      ('comp3', 'comp1', reg.EdgeInteractionType.OR.value),
                      ('A_[c]--C_[a]', 'comp1', reg.EdgeInteractionType.OR.value),
                      ('A_[d]--D_[a]', 'comp2', reg.EdgeInteractionType.AND.value),
                      ('A_[e]--E_[a]', 'comp2', reg.EdgeInteractionType.AND.value),
                      ('A_[f]--F_[a]', 'comp3', reg.EdgeInteractionType.AND.value),
                      ('A_[g]--G_[a]', 'comp3', reg.EdgeInteractionType.AND.value),]
                     ),

        RuleTestCase('''A_[B]_ppi_B_[A]; ! <comp>
                        <comp>; AND <comp1>; AND <Notcomp2>
                        <comp1>; OR <comp3>; OR A--C
                        <Notcomp2>; NOT A--D
                        <comp3>; AND A--F; AND A--G
                        A_[C]_ppi_C_[A]
                        A_[D]_ppi_D_[A]
                        A_[F]_ppi_F_[A]
                        A_[G]_ppi_G_[A]''',
                     ['A_[B]_ppi_B_[A]', 'A_[C]_ppi_C_[A]', 'A_[D]_ppi_D_[A]', 'A_[F]_ppi_F_[A]', 'A_[G]_ppi_G_[A]' ],
                     ['A_[B]--B_[A]', 'A_[B]--0', 'B_[A]--0', 'A_[C]--C_[A]', 'A_[C]--0', 'C_[A]--0', 'A_[D]--D_[A]',
                      'A_[D]--0', 'D_[A]--0', 'A_[F]--F_[A]', 'A_[F]--0', 'F_[A]--0', 'A_[G]--G_[A]', 'A_[G]--0',
                      'G_[A]--0'],
                     ['comp#AND', 'comp1#OR', 'Notcomp2#NOT', 'comp3#AND'],
                     [('A_[B]_ppi_B_[A]', 'A_[B]--B_[A]', reg.EdgeInteractionType.produce.value),
                      ('A_[B]_ppi_B_[A]', 'A_[B]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[B]_ppi_B_[A]', 'B_[A]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[C]_ppi_C_[A]', 'A_[C]--C_[A]', reg.EdgeInteractionType.produce.value),
                      ('A_[C]_ppi_C_[A]', 'A_[C]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[C]_ppi_C_[A]', 'C_[A]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[D]_ppi_D_[A]', 'A_[D]--D_[A]', reg.EdgeInteractionType.produce.value),
                      ('A_[D]_ppi_D_[A]', 'A_[D]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[D]_ppi_D_[A]', 'D_[A]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[F]_ppi_F_[A]', 'A_[F]--F_[A]', reg.EdgeInteractionType.produce.value),
                      ('A_[F]_ppi_F_[A]', 'A_[F]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[F]_ppi_F_[A]', 'F_[A]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[G]_ppi_G_[A]', 'A_[G]--G_[A]', reg.EdgeInteractionType.produce.value),
                      ('A_[G]_ppi_G_[A]', 'A_[G]--0', reg.EdgeInteractionType.consume.value),
                      ('A_[G]_ppi_G_[A]', 'G_[A]--0', reg.EdgeInteractionType.consume.value),
                      ('comp', 'A_[B]_ppi_B_[A]', reg.EdgeInteractionType.required.value),
                      ('comp1', 'comp', reg.EdgeInteractionType.AND.value),
                      ('Notcomp2', 'comp', reg.EdgeInteractionType.AND.value),
                      ('comp3', 'comp1', reg.EdgeInteractionType.OR.value),
                      ('A_[C]--C_[A]', 'comp1', reg.EdgeInteractionType.OR.value),
                      ('A_[D]--D_[A]', 'Notcomp2', reg.EdgeInteractionType.NOT.value),
                      ('A_[F]--F_[A]', 'comp3', reg.EdgeInteractionType.AND.value),
                      ('A_[G]--G_[A]', 'comp3', reg.EdgeInteractionType.AND.value),]
                     ),

        RuleTestCase('''A_p+_B_[(a)]
                        C_p-_B_[(a)]''',
                     ['A_p+_B_[(a)]', 'C_p-_B_[(a)]'],
                     ['B_[(a)]-{p}', 'B_[(a)]-{0}'],
                     [],
                     [('A_p+_B_[(a)]', 'B_[(a)]-{p}', reg.EdgeInteractionType.produce.value),
                      ('A_p+_B_[(a)]', 'B_[(a)]-{0}', reg.EdgeInteractionType.consume.value),
                      ('C_p-_B_[(a)]', 'B_[(a)]-{p}', reg.EdgeInteractionType.consume.value),
                      ('C_p-_B_[(a)]', 'B_[(a)]-{0}', reg.EdgeInteractionType.produce.value)]),


        # RuleTestCase('''A_[B]_ppi_B_[A]; ! [Input]
        #                 [Output]; x A--B''',
        #              ['A_[B]_ppi_B_[A]', '[Output]#out'],
        #              ['A_[B]--B_[A]', '[Input]#in'],
        #              [],
        #              [('A_[B]_ppi_B_[A]', 'A_[B]--B_[A]', reg.EdgeInteractionType.produce.value),
        #               ('[Input]', 'A_[B]_ppi_B_[A]', reg.EdgeInteractionType.required.value),
        #               ('A_[B]--B_[A]', '[Output]', reg.EdgeInteractionType.inhibition.value)]),

        # RuleTestCase('''A_ppi_B; ! <comp>
        #                 <comp>; AND <comp1>; AND [Input]
        #                 <comp1>; OR A--D; OR A--C
        #                 [Output]; x <comp>
        #                 A_ppi_D
        #                 A_ppi_C''',
        #              ['A_ppi_B', 'A_ppi_C', 'A_ppi_D', '[Output]#out'],
        #              ['A--B', 'A--C', 'A--D', '[Input]#in'],
        #              ['comp#AND', 'comp1#OR'],
        #              [('A_ppi_B', 'A--B', reg.EdgeInteractionType.produce.value),
        #               ('A_ppi_C', 'A--C', reg.EdgeInteractionType.produce.value),
        #               ('A_ppi_D', 'A--D', reg.EdgeInteractionType.produce.value),
        #               ('comp', 'A_ppi_B', reg.EdgeInteractionType.required.value),
        #               ('comp', '[Output]', reg.EdgeInteractionType.inhibition.value),
        #               ('[Input]', 'comp', reg.EdgeInteractionType.AND.value),
        #               ('comp1', 'comp', reg.EdgeInteractionType.AND.value),
        #               ('A--D', 'comp1', reg.EdgeInteractionType.OR.value),
        #               ('A--C', 'comp1', reg.EdgeInteractionType.OR.value)]),

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

    #gml_system = gml.XGMML(actual_graph, "test_graph")
    #gml_system.to_file('/home/thiemese/project/rxncon/graphml/test_boolean.xgmml')
    expected_graph = nex.DiGraph()
    [expected_graph.add_node(node, type=reg.NodeType.reaction.value) if '#out' not in node
     else expected_graph.add_node(node.split('#out')[0], type=reg.NodeType.output.value) for node in test_case.reaction_node_strings]

    expected_graph = get_state_nodes(test_case, expected_graph)
    expected_graph = get_boolean_complex_state_nodes(test_case, expected_graph)

    [expected_graph.add_edge(source, target, interaction=interaction) for source, target, interaction in test_case.edge_tuples]

    return expected_graph.node == actual_graph.node and expected_graph.edge == actual_graph.edge
