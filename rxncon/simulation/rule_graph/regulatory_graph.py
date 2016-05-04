import typing as tp
import re
from enum import Enum, unique
import networkx as nex
import rxncon.core.rxncon_system as rxs
import rxncon.core.reaction as rxn
import rxncon.core.effector as eff
import rxncon.core.contingency as con
import rxncon.core.state as sta

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
    componentstate = 'component_state'


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
    component = 'component'
    product_state = 'ps'
    source_state = 'ss'


effector_edge_interaction_type_mapping = {eff.OrEffector: EdgeInteractionType.OR,
                                          eff.AndEffector: EdgeInteractionType.AND,
                                          eff.NotEffector: EdgeInteractionType.NOT,
                                          con.ContingencyType.requirement: EdgeInteractionType.required,
                                          con.ContingencyType.inhibition: EdgeInteractionType.inhibition,
                                          con.ContingencyType.positive: EdgeInteractionType.positive,
                                          con.ContingencyType.negative: EdgeInteractionType.negative
                                          }

# handling Input, Output?
effector_node_type_mapping = { eff.OrEffector: NodeType.OR,
                               eff.AndEffector: NodeType.AND,
                               eff.NotEffector: NodeType.NOT
                               }


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
        if isinstance(reaction, rxn.OutputReaction):
            graph.add_node(str(reaction), dict(type=NodeType.output.value))
            return graph
        elif reaction.product:
            graph.add_node(str(reaction.product), dict(type=NodeType.state.value))
            graph.add_edge(str(reaction), str(reaction.product), interaction=EdgeInteractionType.produce.value)
            return graph
        elif reaction.source:
            for state_effector in self.get_subspecifications_of_state(reaction.source):
                graph.add_edge(str(reaction), str(state_effector.expr), interaction=EdgeInteractionType.consume.value)
            return graph
        else:
            raise AssertionError

    def add_state_effector_edge_or_input_to_graph(self, effector: eff.StateEffector, target_name, edge_type: EdgeInteractionType,
                                                  graph: nex.DiGraph):
        if isinstance(effector.expr, sta.InputState):
            graph.add_node(self.replace_invalid_chars(str(effector.expr)), type=NodeType.input.value)
        graph.add_edge(self.replace_invalid_chars(str(effector.expr)), target_name,
                       interaction=edge_type.value)
        return graph

    def get_subspecifications_of_state(self, state: sta.State) -> tp.List[eff.StateEffector]:
        state_subspecification_effectors = set()
        for react in self.rxncon_system.reactions:
            if react.product is not None and state.is_superspecification_of(react.product):
                state_subspecification_effectors.add(eff.StateEffector(react.product))
            elif react.source is not None and state.is_superspecification_of(react.source):
                state_subspecification_effectors.add(eff.StateEffector(react.source))
        return state_subspecification_effectors

    def add_contingencies_to_graph(self, contingencies: tp.List[con.Contingency], reaction: rxn.Reaction, graph: nex.DiGraph):
        for contingency in contingencies:
            if isinstance(contingency.effector, eff.StateEffector):
                for state_effector in self.get_subspecifications_of_state(contingency.effector.expr):
                    graph = self.add_state_effector_edge_or_input_to_graph(state_effector, str(reaction),
                                                                           effector_edge_interaction_type_mapping[contingency.type], graph)
            elif isinstance(contingency.effector, eff.AndEffector) or isinstance(contingency.effector, eff.OrEffector) \
                    or isinstance(contingency.effector, eff.NotEffector):
                graph = self.add_information_from_effector_reaction_contingency_type_to_graph(contingency.effector, contingency.target,
                                                                                              effector_edge_interaction_type_mapping[contingency.type], graph)

        return graph

    def replace_invalid_chars(self, name: str):
        return re.sub('[<>]', '', name)

    def target_name_from_reaction_or_effector(self, target: tp.Union[eff.Effector, rxn.Reaction]):
        if isinstance(target, rxn.Reaction):
            return self.replace_invalid_chars(str(target))
        elif isinstance(target, eff.Effector):
            return self.replace_invalid_chars(target.name)
        else:
            raise AssertionError

    def add_information_from_effector_reaction_contingency_type_to_graph(self, effector: eff.Effector, target: tp.Union[eff.Effector, rxn.Reaction],
                                                                         edge_type: EdgeInteractionType, graph: nex.DiGraph):

        target_name = self.target_name_from_reaction_or_effector(target)

        if isinstance(effector, eff.AndEffector) or isinstance(effector, eff.OrEffector):
            graph.add_node(self.replace_invalid_chars(effector.name),
                           type=effector_node_type_mapping[type(effector)].value)
            graph.add_edge(self.replace_invalid_chars(effector.name), target_name,
                           interaction=edge_type.value)
            graph = self.add_information_from_effector_reaction_contingency_type_to_graph(effector.left_expr, effector,
                                                                                          effector_edge_interaction_type_mapping[type(effector)], graph)
            graph = self.add_information_from_effector_reaction_contingency_type_to_graph(effector.right_expr, effector,
                                                                                          effector_edge_interaction_type_mapping[type(effector)], graph)
            return graph
        elif isinstance(effector, eff.StateEffector):
            for state_effector in self.get_subspecifications_of_state(effector.expr):
                graph = self.add_state_effector_edge_or_input_to_graph(state_effector, target_name,
                                                                       edge_type, graph)
            return graph
        elif isinstance(effector, eff.NotEffector):
            graph.add_node(self.replace_invalid_chars(effector.name),
                           type=effector_node_type_mapping[type(effector)].value)
            graph.add_edge(self.replace_invalid_chars(effector.name), target_name,
                           interaction=edge_type.value)
            graph = self.add_information_from_effector_reaction_contingency_type_to_graph(effector.expr, effector,
                                                                                          effector_edge_interaction_type_mapping[type(effector)], graph)
            return graph
        else:
            raise AssertionError


class SemanticRegulatoryGraph(RegulatoryGraph):
    def __init__(self, rxncon_system: rxs.RxnConSystem):
        self.rxncon_system = rxncon_system

    def add_reaction_to_graph(self, reaction: rxn.Reaction, graph: nex.DiGraph):
        graph.add_node(str(reaction), dict(type=NodeType.reaction.value))
        if isinstance(reaction, rxn.OutputReaction):
            graph.add_node(str(reaction), dict(type=NodeType.output.value))
            return graph
        elif reaction.product:
            graph.add_node(str(reaction.product), dict(type=NodeType.state.value))
            graph.add_edge(str(reaction), str(reaction.product), interaction=EdgeInteractionType.produce.value)
            graph.add_edge(str(reaction.product), str(reaction), interaction=EdgeInteractionType.product_state.value)
            return graph
        elif reaction.source:
            graph.add_edge(str(reaction), str(reaction.source), interaction=EdgeInteractionType.consume.value)
            graph.add_edge(str(reaction.source), str(reaction), interaction=EdgeInteractionType.source_state.value)
            return graph
        else:
            raise AssertionError

class BooleanReagulatoryGraph(RegulatoryGraph):
    def __init__(self, rxncon_system: rxs.RxnConSystem):
        self.rxncon_system = rxncon_system

    def add_component_states_to_graph(self, reaction: rxn.Reaction, graph: nex.DiGraph):
        #todo: check if a component is produced. As soon as it gets produced change the node type
        #todo: from componentstate to state
        if str(sta.ComponentState(reaction.subject.to_component_specification())) not in graph.nodes():
            graph.add_node(str(sta.ComponentState(reaction.subject.to_component_specification())),
                           dict(type=NodeType.componentstate.value))

        if str(sta.ComponentState(reaction.object.to_component_specification())) not in graph.nodes():
            graph.add_node(str(sta.ComponentState(reaction.object.to_component_specification())),
                           dict(type=NodeType.componentstate.value))
        graph.add_edge(str(sta.ComponentState(reaction.subject.to_component_specification())),
                           str(reaction), interaction=EdgeInteractionType.component.value)
        graph.add_edge(str(sta.ComponentState(reaction.object.to_component_specification())),
                           str(reaction), interaction=EdgeInteractionType.component.value)

        return graph

    def add_reaction_to_graph(self, reaction: rxn.Reaction, graph: nex.DiGraph):
        graph.add_node(str(reaction), dict(type=NodeType.reaction.value))
        if isinstance(reaction, rxn.OutputReaction):
            graph.add_node(str(reaction), dict(type=NodeType.output.value))
            return graph

        graph = self.add_component_states_to_graph(reaction, graph)

        if reaction.product:
            graph.add_node(str(reaction.product), dict(type=NodeType.state.value))
            graph.add_edge(str(reaction), str(reaction.product), interaction=EdgeInteractionType.produce.value)
            graph.add_edge(str(reaction.product), str(reaction), interaction=EdgeInteractionType.product_state.value)
            return graph
        if reaction.source:
            graph.add_node(str(reaction.source), dict(type=NodeType.state.value))
            graph.add_edge(str(reaction), str(reaction.source), interaction=EdgeInteractionType.consume.value)
            graph.add_edge(str(reaction.source), str(reaction), interaction=EdgeInteractionType.source_state.value)
            return graph
        else:
            raise AssertionError