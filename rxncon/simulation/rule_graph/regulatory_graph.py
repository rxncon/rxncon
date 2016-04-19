import typing as tp
import re
from enum import Enum, unique
import networkx as nex
import rxncon.core.rxncon_system as rxs
import rxncon.core.reaction as rxn
import rxncon.core.effector as eff
import rxncon.core.contingency as con

# for cytoscape export:
# reaction graph of rxncon: visalisation of specifications and there relationships
# for each spec one node of certain size
# for each domain/subdomain/residue one node of certain size with edge length 0 to it's main node,
# the length of a edge should be defined by it's weight since there is no algorithm considering the length, which makes
# the length of an edge dependent on the node position only

#@unique
class NodeType(Enum):
    reaction = 'reaction'
    state = 'state'
    AND = "boolean_and"
    OR = "boolean_or"
    NOT = "boolean_not"
    input = 'input'
    output = 'output'


@unique
class EdgeInteractionType(Enum):
    produce = 'produce'
    consume = 'consume'
    required = '!'
    inhibition = "x"
    positive    = 'k+'
    negative    = 'k-'
    AND = 'AND'
    OR = 'OR'
    NOT = 'NOT'


effector_edge_interaction_type_mapping = {eff.OrEffector: EdgeInteractionType.OR,
                                          eff.AndEffector: EdgeInteractionType.AND,
                                          eff.NotEffector: EdgeInteractionType.NOT}


effocor_node_type_mapping = {eff.OrEffector: NodeType.OR,
                             eff.AndEffector: NodeType.AND,
                             eff.NotEffector: NodeType.NOT}


class RegulatoryGraph():
    def __init__(self,rxncon_system: rxs.RxnConSystem):
        self.rxncon_system = rxncon_system

    def to_graph(self):
        graph = nex.DiGraph()
        for reaction in self.rxncon_system.reactions:
            graph = self.add_reaction_to_graph(reaction, graph)

            contingencies = self.rxncon_system.strict_contingencies_for_reaction(reaction)
            contingencies.extend(self.rxncon_system.quantitative_contingencies_for_reaction(reaction))
            graph = self.add_contingencies_to_graph(contingencies, reaction, graph)
        return graph

    def add_reaction_to_graph(self, reaction: rxn.Reaction, graph: nex.DiGraph):
        graph.add_node(str(reaction), dict(type=NodeType.reaction.value))
        if reaction.product:
            graph.add_node(str(reaction.product), dict(type=NodeType.state.value))
            graph.add_edge(str(reaction), str(reaction.product), interaction=EdgeInteractionType.produce.value)
        elif reaction.source:
            graph.add_edge(str(reaction), str(reaction.source), interaction=EdgeInteractionType.consume.value)
        return graph

    def add_contingencies_to_graph(self, contingencies: tp.List[con.Contingency], reaction: rxn.Reaction, graph: nex.DiGraph):
        for contingency in contingencies:
            if isinstance(contingency.effector, eff.StateEffector):
                graph.add_edge(str(contingency.effector.expr), str(reaction), interaction=contingency.type.value)
            elif isinstance(contingency.effector, eff.AndEffector) or isinstance(contingency.effector, eff.OrEffector):
                graph = self.add_edges_from_contingency(contingency, graph)
        return graph

    def replace_invalid_chars(self, name: str):
        return re.sub('[<>]', '', name)

    def add_edges_from_contingency(self, contingency: con.Contingency, graph: nex.DiGraph):
        if isinstance(contingency.target, rxn.Reaction):
            if isinstance(contingency.effector, eff.AndEffector) or isinstance(contingency.effector, eff.OrEffector):
                graph.add_node(self.replace_invalid_chars(contingency.effector.name),
                               type=effocor_node_type_mapping[type(contingency.effector)].value)
                graph.add_edge(self.replace_invalid_chars(contingency.effector.name), str(contingency.target),
                               interaction=contingency.type.value)
                graph = self.add_edges_from_effector(contingency.effector.left_expr, contingency.effector, graph)
                graph = self.add_edges_from_effector(contingency.effector.right_expr, contingency.effector, graph)
                return graph
            elif isinstance(contingency.effector, eff.StateEffector):
                graph.add_edge(self.replace_invalid_chars(str(contingency.effector.expr)), str(contingency.target),
                               interaction=contingency.type.value)
                return graph
        else:
            raise AssertionError

    def add_edges_from_effector(self, effector: eff.Effector, root, graph: nex.DiGraph):
        if isinstance(effector, eff.StateEffector):
            graph.add_edge(self.replace_invalid_chars(str(effector.expr)), self.replace_invalid_chars(root.name),
                           interaction=effector_edge_interaction_type_mapping[type(root)].value)
            return graph
        elif isinstance(effector, eff.AndEffector) or isinstance(effector, eff.OrEffector):
            graph.add_node(self.replace_invalid_chars(effector.name),
                           type=effocor_node_type_mapping[type(effector)].value)
            graph.add_edge(self.replace_invalid_chars(effector.name), self.replace_invalid_chars(root.name),
                           interaction=effector_edge_interaction_type_mapping[type(root)].value)
            graph = self.add_edges_from_effector(effector.left_expr, effector, graph)
            graph = self.add_edges_from_effector(effector.right_expr, effector, graph)
            return graph
        elif isinstance(effector, eff.NotEffector):
            graph.add_node(self.replace_invalid_chars(effector.name),
                           type=effocor_node_type_mapping[type(effector)].value)
            graph.add_edge(self.replace_invalid_chars(effector.name), self.replace_invalid_chars(root.name),
                           interaction=effector_edge_interaction_type_mapping[type(root)].value)
            graph = self.add_edges_from_effector(effector.expr, effector, graph)
            return graph

        else:
            raise NotImplementedError
