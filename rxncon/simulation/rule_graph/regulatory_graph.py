import typing as tp
import re
import copy
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

#todo: subspecification

class RegulatoryGraph():
    def __init__(self,rxncon_system: rxs.RxnConSystem):
        self.rxncon_system = rxncon_system
        self.graph = nex.DiGraph()

    def to_graph(self):
        for reaction in self.rxncon_system.reactions:
            self.add_reaction_information_to_graph(reaction)
            contingencies = self._get_contingencies(reaction)
            self.add_contingencies_to_graph(contingencies, reaction)
        return self.graph

    def _get_contingencies(self, reaction: rxn.Reaction):
        contingencies = self.rxncon_system.strict_contingencies_for_reaction(reaction)
        contingencies.extend(self.rxncon_system.quantitative_contingencies_for_reaction(reaction))
        return contingencies

    def _add_reaction_reactant_to_graph(self, reaction: rxn.Reaction, reactants: rxn.Reactant, edge_type: EdgeInteractionType):
        for reactant_state in reactants.states:
            self.graph.add_node(str(reactant_state), dict(type=NodeType.state.value))
            self.graph.add_edge(str(reaction), str(reactant_state), interaction=edge_type.value)

    def add_reaction_information_to_graph(self, reaction: rxn.Reaction):
        self.graph.add_node(str(reaction), dict(type=NodeType.reaction.value))

        for reactant_post in reaction.reactants_rhs:
            self._add_reaction_reactant_to_graph(reaction, reactant_post, EdgeInteractionType.produce)

        for reactant_pre in reaction.reactants_lhs:
            self._add_reaction_reactant_to_graph(reaction, reactant_pre, EdgeInteractionType.consume)

    def remove_structure(self, state):
        state = copy.deepcopy(state)
        for var in state.variables:
            state.variables[var].structure_index = None
        return state

    def add_state_effector_edge_or_input_to_graph(self, effector: eff.StateEffector, target_name, edge_type: EdgeInteractionType):

        self.graph.add_edge(self.replace_invalid_chars(str(self.remove_structure(effector.expr))), target_name,
                            interaction=edge_type.value)

    def _check_reactant_subset_of_state(self, state: sta.State, reactant: rxn.Reactant, subset_states: tp.Set):
        for reactant_state in reactant.states:
            if state.is_superset_of(reactant_state):
                subset_states.add(eff.StateEffector(reactant_state))
        return subset_states

    def get_subset_of_state(self, state: sta.State) -> tp.List[eff.StateEffector]:
        subset_states = set()
        for reaction in self.rxncon_system.reactions:
            for reactant_post in reaction.reactants_rhs + reaction.reactants_lhs:
                subset_states = self._check_reactant_subset_of_state(state, reactant_post, subset_states)
        return subset_states

    def add_contingencies_to_graph(self, contingencies: tp.List[con.Contingency], reaction: rxn.Reaction):
        for contingency in contingencies:
            if isinstance(contingency.effector, eff.StateEffector):
                for state_effector in self.get_subset_of_state(contingency.effector.expr):
                    self.add_state_effector_edge_or_input_to_graph(state_effector, str(reaction),
                                                                           effector_edge_interaction_type_mapping[contingency.type])
            elif isinstance(contingency.effector, eff.AndEffector) or isinstance(contingency.effector, eff.OrEffector) \
                    or isinstance(contingency.effector, eff.NotEffector):
                        self.add_information_from_effector_reaction_contingency_type_to_graph(contingency.effector, self.target_name_from_reaction_or_effector(contingency.target),
                                                                                              effector_edge_interaction_type_mapping[contingency.type])

    def replace_invalid_chars(self, name: str):
        return re.sub('[<>]', '', name)

    def target_name_from_reaction_or_effector(self, target: tp.Union[eff.Effector, rxn.Reaction]):
        if isinstance(target, rxn.Reaction):
            return self.replace_invalid_chars(str(target))
        elif isinstance(target, eff.Effector):
            return self.replace_invalid_chars(target.name)
        else:
            raise AssertionError

    def add_information_from_effector_reaction_contingency_type_to_graph(self, effector: eff.Effector, target_name,
                                                                         edge_type: EdgeInteractionType):

        if isinstance(effector, eff.AndEffector) or isinstance(effector, eff.OrEffector):
            if effector.name is not None:
                self.graph.add_node(self.replace_invalid_chars(effector.name),
                                    type=effector_node_type_mapping[type(effector)].value)
                self.graph.add_edge(self.replace_invalid_chars(effector.name), target_name,
                                    interaction=edge_type.value)
                target_name = self.target_name_from_reaction_or_effector(effector)

            self.add_information_from_effector_reaction_contingency_type_to_graph(effector.left_expr, target_name,
                                                                                  effector_edge_interaction_type_mapping[type(effector)])
            self.add_information_from_effector_reaction_contingency_type_to_graph(effector.right_expr, target_name,
                                                                                  effector_edge_interaction_type_mapping[type(effector)])
        elif isinstance(effector, eff.StateEffector):
            for state_effector in self.get_subset_of_state(effector.expr):
                self.add_state_effector_edge_or_input_to_graph(state_effector, target_name, edge_type)
        elif isinstance(effector, eff.NotEffector):
            self.graph.add_node(self.replace_invalid_chars(effector.name),
                                type=effector_node_type_mapping[type(effector)].value)
            self.graph.add_edge(self.replace_invalid_chars(effector.name), target_name,
                                interaction=edge_type.value)
            target_name = self.target_name_from_reaction_or_effector(effector)
            self.add_information_from_effector_reaction_contingency_type_to_graph(effector.expr, target_name,
                                                                                  effector_edge_interaction_type_mapping[type(effector)])
        else:
            raise AssertionError


