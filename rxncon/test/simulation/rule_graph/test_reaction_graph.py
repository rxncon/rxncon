import pytest
import networkx as nex
from collections import namedtuple
import rxncon.input.quick.quick as qui
import rxncon.simulation.rule_graph.reaction_graph as rxng
import rxncon.simulation.rule_graph.graphML as gml

# TestCase = namedtuple('TestCase', ['quick_string', 'component_node_strings', 'domain_node_strings',
#                                            'residue_node_strings', 'edge_tuples'])
#
#
#
# def test_reaction_graph.py(case_and_expected_reaction_graph):
#
#     for the_case in case_and_expected_reaction_graph:
#         actual_system = qui.Quick(the_case.quick_string)
#         reg_system = rxng.ReactionGraph(actual_system.rxncon_system)
#         actual_graph = reg_system.to_graph()
#         return is_reaction_graph_correct(the_case)
#
# @pytest.fixture
# def case_and_expected_reaction_graph():
#     return [
#             TestCase('A_[b]_ppi_B_[a]',
#                      ['A', 'B'],
#                      ['b', 'a'],
#                      [],
#                      ('A', 'b', rxng.EdgeType.internal, rxng.EdgeLength.short),
#                      ('B', 'a', rxng.EdgeType.internal, rxng.EdgeLength.short),
#                      ('b', 'a', rxng.EdgeType.interaction, rxng.EdgeLength.long))
#     ]
#
#
# def is_reaction_graph_correct(the_case):
#     # gml_system = gml.XGMML(actual_graph, "test_graph")
#     # gml_system.to_file('/home/thiemese/project/rxncon/graphml/test_boolean.xgmml')
#     expected_graph = nex.DiGraph()
#     #[expected_graph.add_node(node, dict(size=NodeType.state.value)) for node in the_case.component_node_strings]

def test_a():
    pass

def test_regulatory_graph_generation():
    print('adad')

def test_test():
    graph = nex.DiGraph()
    graph.add_node('A', [dict(name='size', value=40.0), dict(name='name', value='A')])
    graph.add_node('b', [dict(name='size', value=20.0), dict(name='name', value='b')])
    graph.add_node('B', [dict(name='size', value=40.0), dict(name='name', value='B')])
    graph.add_node('a', [dict(name='size', value=20.0), dict(name='name', value='a')])
    graph.add_edge('a', 'b', interaction='i')
    gml_system = gml.XGMML(graph, "test_graph")
    pass


