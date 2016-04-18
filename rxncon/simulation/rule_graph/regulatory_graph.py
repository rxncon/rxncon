import typing as tp
from enum import Enum, unique
import networkx as nex
import rxncon.core.rxncon_system as rxs
import rxncon.core.reaction as rxn
import rxncon.venntastic.sets as venn

from rxncon.simulation.bBM.bbm_from_rxncon import _state_set_from_contingencies

@unique
class NodeType(Enum):
    reaction = 'reaction'
    state = 'state'

@unique
class EdgeInteractionType(Enum):
    produce = 'produce'
    consume = 'consume'
    required = '!'
    inhibition = "x"
    positive    = 'k+'
    negative    = 'k-'


class RegulatoryGraph():
    def __init__(self,rxncon_system: rxs.RxnConSystem):
        self.rxncon_system = rxncon_system

    def to_graph(self):
        graph = nex.DiGraph()
        for reaction in self.rxncon_system.reactions:
            graph.add_node(str(reaction), dict(type=NodeType.reaction.value))
            if reaction.product:
                graph.add_node(str(reaction.product), dict(type=NodeType.state.value))
                graph.add_edge(str(reaction), str(reaction.product), interaction=EdgeInteractionType.produce.value)
            elif reaction.source:
                graph.add_edge(str(reaction), str(reaction.source), interaction=EdgeInteractionType.consume.value)

            contingency_state_set = self.contingencies_from_rxncon_system_and_reaction(reaction)
            if contingency_state_set:
                graph = self.add_contingency_edges_to_regulatory_graph(contingency_state_set, reaction, graph)
        return graph

    def add_contingency_edges_to_regulatory_graph(self, vennset: venn.Set, reaction:rxn.Reaction, graph: nex.Graph):

        if isinstance(vennset, venn.EmptySet):
            raise AssertionError
        if isinstance(vennset, venn.PropertySet):
            #todo add node
            graph.add_edge(str(vennset.value), str(reaction), interaction=EdgeInteractionType.required.value)
            return graph
            #return venn.PropertySet(vennset.value)
        if isinstance(vennset, venn.Complement):
            #todo add node complement and edges
            return venn.Complement(self.add_contingency_edges_to_regulatory_graph(vennset.expr, reaction, graph))
        elif isinstance(vennset, venn.Intersection):
            #todo: G.add_node(boolean_and) and edges
            return venn.Intersection(self.add_contingency_edges_to_regulatory_graph(vennset.left_expr, reaction, graph), self.add_contingency_edges_to_regulatory_graph(vennset.right_expr, reaction, graph))
        elif isinstance(vennset, venn.Union):
            #todo: G.add_node(boolean_or) and edges
            return venn.Union(self.add_contingency_edges_to_regulatory_graph(vennset.left_expr, reaction, graph), self.add_contingency_edges_to_regulatory_graph(vennset.right_expr, reaction, graph))
        else:
            raise NotImplementedError

    def contingencies_from_rxncon_system_and_reaction(self, reaction: rxn.Reaction):
        strict_contingency_state_set = _state_set_from_contingencies(self.rxncon_system.strict_contingencies_for_reaction(reaction))
        if isinstance(strict_contingency_state_set.to_full_simplified_form(), venn.EmptySet):
            raise AssertionError("There is no way to fulfill the contingencies: {}".format(strict_contingency_state_set))

        additional_contingency_state_set = _state_set_from_contingencies(self.rxncon_system.quantitative_contingencies_for_reaction(reaction))

        if isinstance(additional_contingency_state_set.to_full_simplified_form(), venn.EmptySet):
            raise AssertionError("There is no way to fulfill the contingencies: {}".format(additional_contingency_state_set))
        vennset = venn.Intersection(strict_contingency_state_set.to_full_simplified_form(), additional_contingency_state_set.to_full_simplified_form())
        return vennset.to_full_simplified_form()
