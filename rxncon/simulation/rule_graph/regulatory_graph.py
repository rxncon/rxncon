from typing import Union, List
import re
from enum import Enum, unique
from networkx import DiGraph
from rxncon.core.rxncon_system import RxnConSystem, Reaction, ReactionTerm, Contingency
from rxncon.core.effector import StateEffector, AndEffector, NotEffector, OrEffector, Effector


@unique
class NodeType(Enum):
    reaction = 'reaction'
    state    = 'state'
    AND      = "boolean_and"
    OR       = "boolean_or"
    NOT      = "boolean_not"
    input    = 'input'
    output   = 'output'

@unique
class EdgeInteractionType(Enum):
    produce       = 'produce'
    consume       = 'consume'
    required      = '!'
    inhibition    = "x"
    positive      = 'k+'
    negative      = 'k-'
    AND           = 'AND'
    OR            = 'OR'
    NOT           = 'NOT'
    source_state  = 'ss'



class RegulatoryGraph():
    """
    Definition of the regulatory Graph.

    Args:
        rxncon_system: The rxncon system.
    """
    def __init__(self,rxncon_system: RxnConSystem) -> None:
        self.rxncon_system = rxncon_system
        self.graph = DiGraph()

    def to_graph(self) -> DiGraph:
        """
        Creates the regulatory graph.

        Returns:
            DiGraph: The regulatory graph as networkx DiGraph.

        """
        for reaction in self.rxncon_system.reactions:
            self.add_reaction_information_to_graph(reaction)
            self.add_contingency_information_to_graph(self.rxncon_system.contingencies_for_reaction(reaction))
        return self.graph

    def add_reaction_information_to_graph(self, reaction: Reaction) -> None:
        """
        Adds reaction information to the graph.

        Args:
            reaction: A rxncon reaction.

        Mutates:
            self.graph (DiGraph): The regulatory graph.

        Returns:
            None

        """

        def _add_reaction_reactant_to_graph(reaction: Reaction, reactants: ReactionTerm,
                                            edge_type: EdgeInteractionType) -> None:

            """
            Adds the reactants of a reaction to the graph.

            Args:
                reaction: A rxncon reaction.
                reactants: The reactants of the respective reaction.
                edge_type: The type of the edge can be either consuming or producing.

            Mutates:
                self.graph (DiGraph): The regulatory graph.

            Returns:
                None

            """
            for reactant_state in reactants.states:
                self.graph.add_node(str(reactant_state), dict(type=NodeType.state.value))
                self.graph.add_edge(str(reaction), str(reactant_state), interaction=edge_type.value)

        def _add_reaction_source_state_to_graph(reaction: Reaction, reactants: ReactionTerm):
            for reactant_state in reactants.states:
                self.graph.add_edge(str(reactant_state), str(reaction), interaction=EdgeInteractionType.source_state.value)

        self.graph.add_node(str(reaction), dict(type=NodeType.reaction.value))

        for reactant_post in reaction.terms_rhs:
            _add_reaction_reactant_to_graph(reaction, reactant_post, EdgeInteractionType.produce)

        for reactant_pre in reaction.terms_lhs:
            _add_reaction_reactant_to_graph(reaction, reactant_pre, EdgeInteractionType.consume)
            _add_reaction_source_state_to_graph(reaction, reactant_pre)

    def add_contingency_information_to_graph(self, contingencies: List[Contingency]) -> None:
        for contingency in contingencies:
            self._add_information_from_effector_to_graph(contingency.effector, contingency.type,
                                                         self._target_name_from_reaction_or_effector(contingency.target)
                                                         )

    def _replace_invalid_chars(self, name: str) -> str:
        return re.sub('[<>]', '', name)

    def _add_node_and_edge(self, name: str, node_type: NodeType, edge_type: EdgeInteractionType, target_name: str) -> None:
        self.graph.add_node(self._replace_invalid_chars(str(name)),
                            type=node_type.value)

        self.graph.add_edge(self._replace_invalid_chars(str(name)), target_name,
                            interaction=edge_type.value)

    def _target_name_from_reaction_or_effector(self, target: Union[Effector, Reaction]) -> str:
        if isinstance(target, Reaction):
            return self._replace_invalid_chars(str(target))
        elif isinstance(target, Effector):
            return self._replace_invalid_chars(target.name)
        else:
            raise AssertionError

    def _add_information_from_effector_to_graph(self, effector, edge_type, target_name) -> None:
        if isinstance(effector, StateEffector):
            name = str(effector.expr.to_non_structured())
            if re.match("^\[[/w]*\]?", name):
                self._add_node_and_edge(name, NodeType.input, edge_type, target_name)
            else:
                self._add_node_and_edge(name, NodeType.state, edge_type, target_name)

        elif isinstance(effector, NotEffector):
            self._add_node_and_edge(str(effector.name), NodeType.NOT, edge_type, target_name)
            target_name = self._target_name_from_reaction_or_effector(effector)
            self._add_information_from_effector_to_graph(effector.expr, EdgeInteractionType.NOT, target_name)

        elif isinstance(effector, OrEffector):
            if effector.name is not None:
                self._add_node_and_edge(effector.name, NodeType.OR, edge_type, target_name)
                target_name = self._target_name_from_reaction_or_effector(effector)

            self._add_information_from_effector_to_graph(effector.left_expr, EdgeInteractionType.OR, target_name)
            self._add_information_from_effector_to_graph(effector.right_expr, EdgeInteractionType.OR, target_name)

        elif isinstance(effector, AndEffector):
            if effector.name is not None:
                self._add_node_and_edge(effector.name, NodeType.AND, edge_type, target_name)
                target_name = self._target_name_from_reaction_or_effector(effector)

            self._add_information_from_effector_to_graph(effector.left_expr, EdgeInteractionType.AND, target_name)
            self._add_information_from_effector_to_graph(effector.right_expr, EdgeInteractionType.AND, target_name)
        else:
            raise AssertionError
