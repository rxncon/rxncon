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

        # RuleTestCase('''A_[b]_ppi_B_[a]; ! <comp>; ! C-{p}
        #                 <comp>; AND A-{p}; AND A_[c]--C_[a]
        #                 A_[b]_ppi_C_[a]
        #                 C_p+_A_[(c)]
        #                 D_p+_C_[(d)]''',
        #              ['A_[b]_ppi_B_[a]', 'A_[c]_ppi_C_[a]', 'C_p+_A_[(c)]', 'D_p+_C_[(d)]'],
        #              ['A_[b]--B_[a]', 'A_[c]--C_[a]', 'A_[(c)]-{p}', 'C_[(d)]-{p}'],
        #              ['comp#AND'],
        #              [('A_[b]_ppi_B_[a]', 'A_[b]--B_[a]', reg.EdgeInteractionType.produce.value),
        #               ('A_[c]_ppi_C_[a]', 'A_[c]--C_[a]', reg.EdgeInteractionType.produce.value),
        #               ('C_p+_A_[(c)]', 'A-{p}', reg.EdgeInteractionType.produce.value),
        #               ('D_p+_C_[(d)]', 'C_[(d)]-{p}', reg.EdgeInteractionType.produce.value),
        #               ('C_[(d)]-{p}', 'A_[b]_ppi_B_[a]', reg.EdgeInteractionType.required.value),
        #               ('comp', 'A_[b]_ppi_B_[a]', reg.EdgeInteractionType.required.value),
        #               ('A_[(c)]-{p}', 'comp', reg.EdgeInteractionType.AND.value),
        #               ('A_[c]--C_[a]', 'comp', reg.EdgeInteractionType.AND.value)]),
        #
        # RuleTestCase('''A_[b]_ppi_B_[a]; ! <comp>
        #                 <comp>; AND <comp1>; AND <comp2>
        #                 <comp1>; OR <comp3>; OR A_[c]--C_[a]
        #                 <comp2>; AND A_[d]--D_[a]; AND A_[e]--E_[a]
        #                 <comp3>; AND A_[f]--F_[a]; AND A_[g]--G_[a]
        #                 A_[c]_ppi_C_[a]
        #                 A_[d]_ppi_D_[a]
        #                 A_[e]_ppi_E_[a]
        #                 A_[f]_ppi_F_[a]
        #                 A_[g]_ppi_G_[a]''',
        #              ['A_[b]_ppi_B_[a]', 'A_[c]_ppi_C_[a]', 'A_[d]_ppi_D_[a]', 'A_[e]_ppi_E_[a]', 'A_[f]_ppi_F_[a]', 'A_[g]_ppi_G_[a]' ],
        #              ['A_[b]--B_[a]', 'A_[c]--C_[a]', 'A_[d]--D_[a]', 'A_[e]--E_[a]', 'A_[f]--F_[a]', 'A_[g]--G_[a]'],
        #              ['comp#AND', 'comp1#OR', 'comp2#AND', 'comp3#AND'],
        #              [('A_[b]_ppi_B_[a]', 'A_[b]--B_[a]', reg.EdgeInteractionType.produce.value),
        #               ('A_[c]_ppi_C_[a]', 'A_[c]--C_[a]', reg.EdgeInteractionType.produce.value),
        #               ('A_[d]_ppi_D_[a]', 'A_[d]--D_[a]', reg.EdgeInteractionType.produce.value),
        #               ('A_[e]_ppi_E_[a]', 'A_[e]--E_[a]', reg.EdgeInteractionType.produce.value),
        #               ('A_[f]_ppi_F_[a]', 'A_[f]--F_[a]', reg.EdgeInteractionType.produce.value),
        #               ('A_[g]_ppi_G_[a]', 'A_[g]--G_[a]', reg.EdgeInteractionType.produce.value),
        #               ('comp', 'A_[b]_ppi_B_[a]', reg.EdgeInteractionType.required.value),
        #               ('comp1', 'comp', reg.EdgeInteractionType.AND.value),
        #               ('comp2', 'comp', reg.EdgeInteractionType.AND.value),
        #               ('comp3', 'comp1', reg.EdgeInteractionType.OR.value),
        #               ('A_[c]--C_[a]', 'comp1', reg.EdgeInteractionType.OR.value),
        #               ('A_[d]--D_[a]', 'comp2', reg.EdgeInteractionType.AND.value),
        #               ('A_[e]--E_[a]', 'comp2', reg.EdgeInteractionType.AND.value),
        #               ('A_[f]--F_[a]', 'comp3', reg.EdgeInteractionType.AND.value),
        #               ('A_[g]--G_[a]', 'comp3', reg.EdgeInteractionType.AND.value),]
        #              ),

        # RuleTestCase('''A_[b]_ppi_B; ! <comp>
        #                 <comp>; AND <comp1>; AND <Notcomp2>
        #                 <comp1>; OR <comp3>; OR A--C
        #                 <Notcomp2>; NOT A--D
        #                 <comp3>; AND A--F; AND A--G
        #                 A_ppi_C
        #                 A_ppi_D
        #                 A_ppi_F
        #                 A_ppi_G''',
        #              ['A_[b]_ppi_B', 'A_ppi_C', 'A_ppi_D', 'A_ppi_F', 'A_ppi_G' ],
        #              ['A_[b]--B', 'A--C', 'A--D', 'A--F', 'A--G'],
        #              ['comp#AND', 'comp1#OR', 'Notcomp2#NOT', 'comp3#AND'],
        #              [('A_[b]_ppi_B', 'A_[b]--B', reg.EdgeInteractionType.produce.value),
        #               ('A_ppi_C', 'A--C', reg.EdgeInteractionType.produce.value),
        #               ('A_ppi_D', 'A--D', reg.EdgeInteractionType.produce.value),
        #               ('A_ppi_F', 'A--F', reg.EdgeInteractionType.produce.value),
        #               ('A_ppi_G', 'A--G', reg.EdgeInteractionType.produce.value),
        #               ('comp', 'A_[b]_ppi_B', reg.EdgeInteractionType.required.value),
        #               ('comp1', 'comp', reg.EdgeInteractionType.AND.value),
        #               ('Notcomp2', 'comp', reg.EdgeInteractionType.AND.value),
        #               ('comp3', 'comp1', reg.EdgeInteractionType.OR.value),
        #               ('A--C', 'comp1', reg.EdgeInteractionType.OR.value),
        #               ('A--D', 'Notcomp2', reg.EdgeInteractionType.NOT.value),
        #               ('A--F', 'comp3', reg.EdgeInteractionType.AND.value),
        #               ('A--G', 'comp3', reg.EdgeInteractionType.AND.value),]
        #              ),
        # RuleTestCase('''A_ppi_B; ! [Input]
        #                 [Output]; x A--B''',
        #              ['A_ppi_B', '[Output]#out'],
        #              ['A--B', '[Input]#in'],
        #              [],
        #              [('A_ppi_B', 'A--B', reg.EdgeInteractionType.produce.value),
        #               ('[Input]', 'A_ppi_B', reg.EdgeInteractionType.required.value),
        #               ('A--B', '[Output]', reg.EdgeInteractionType.inhibition.value)]),

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

        # RuleTestCase('''A_p+_B
        #                 C_p-_B''',
        #              ['A_p+_B', 'C_p-_B'],
        #              ['B-{p}'],
        #              [],
        #              [('A_p+_B', 'B-{p}', reg.EdgeInteractionType.produce.value),
        #               ('C_p-_B', 'B-{p}', reg.EdgeInteractionType.consume.value)]),#
        # RuleTestCase('''A_p+_B
        #                 A_p+_B_[d]
        #                 E_p+_B_[(r)]
        #                 C_ppi_B; ! B-{P}
        #                 C_p-_B''',
        #              ['A_p+_B', 'A_p+_B_[d]', 'E_p+_B_[(r)]', 'C_ppi_B', 'C_p-_B'],
        #              ['B-{p}', 'B_[d]-{p}', 'B_[(r)]-{p}', 'C--B'],
        #              [],
        #              [('A_p+_B', 'B-{p}', reg.EdgeInteractionType.produce.value),
        #               ('A_p+_B_[d]', 'B_[d]-{p}', reg.EdgeInteractionType.produce.value),
        #               ('E_p+_B_[(r)]', 'B_[(r)]-{p}', reg.EdgeInteractionType.produce.value),
        #               ('C_p-_B', 'B-{p}', reg.EdgeInteractionType.consume.value),
        #               ('C_p-_B', 'B_[d]-{p}', reg.EdgeInteractionType.consume.value),
        #               ('C_p-_B', 'B_[(r)]-{p}', reg.EdgeInteractionType.consume.value),
        #               ('C_ppi_B', 'C--B', reg.EdgeInteractionType.produce.value),
        #               ('B-{p}', 'C_ppi_B', reg.EdgeInteractionType.required.value),
        #               ('B_[d]-{p}', 'C_ppi_B', reg.EdgeInteractionType.required.value),
        #               ('B_[(r)]-{p}', 'C_ppi_B', reg.EdgeInteractionType.required.value)
        #               ]
        #              )

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
