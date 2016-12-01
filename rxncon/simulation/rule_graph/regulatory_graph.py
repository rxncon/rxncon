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
        def update_degradations_with_contingencies():
            for reaction_target, contingency in reaction_target_to_factor.items():
                # First case: reaction with a non-trivial contingency should degrade only the states appearing
                # in the contingency that are connected to the degraded component.
                if reaction_target.degraded_components and not contingency_factor.is_equivalent_to(UniversalSet()):
                    for degraded_component in reaction_target.degraded_components:
                        if ComponentStateTarget(degraded_component) in component_state_targets:
                            reaction_target.degraded_targets.append(ComponentStateTarget(degraded_component))
                        else:
                            reaction_target.degraded_targets += [state_target for state_target in contingency_factor.values
                                                                 if degraded_component in state_target.components]

                # Second case: reaction with a trivial contingency should degrade all states for the degraded component.
                elif reaction_target.degraded_components and contingency_factor.is_equivalent_to(UniversalSet()):
                    for degraded_component in reaction_target.degraded_components:
                        if ComponentStateTarget(degraded_component) in component_state_targets:
                            reaction_target.degraded_targets.append(ComponentStateTarget(degraded_component))
                        else:
                            reaction_target.degraded_targets += \
                                [StateTarget(x) for x in rxncon_sys.states_for_component(degraded_component)]
        #
        # def update_degradations_for_interaction_states():
        #     new_reactions = {}
        #
        #     for reaction_target, contingency_factor in reaction_target_to_factor.items():
        #         for num, interaction_state in enumerate(state for state in rxncon_sys.states if len(state.components) > 1
        #                                                 and reaction_target.degrades(StateTarget(state))):
        #             neutral_targets = StateTarget(interaction_state).neutral_targets
        #             new_reaction_target = deepcopy(reaction_target)
        #
        #             new_reaction_target.consumed_targets.append(StateTarget(interaction_state))
        #             new_reaction_target.produced_targets += \
        #                 [x for x in neutral_targets if not any(component in new_reaction_target.degraded_components for component in x.components)]
        #
        #             new_reaction_target.interaction_variant_index = num + 1
        #             new_reactions[new_reaction_target] = Intersection(contingency_factor, ValueSet(StateTarget(interaction_state)))

        def _add_reaction_reactant_to_graph(reaction: Reaction, reactants: ReactionTerm,
                                            edge_type: EdgeInteractionType) -> None:

            """
            Adds the reactants of a reaction to the graph.

            Note:
                Adds nodes and edges to the graph.

            Args:
                reaction: A rxncon reaction.
                reactants: The reactants of the respective reaction.
                edge_type: The type of the edge can be either consuming or producing.

            Mutates:
                self.graph (DiGraph): The regulatory graph.

            Returns:
                None

            """
            if not reaction.degraded_components:
                for reactant_state in reactants.states:
                    self.graph.add_node(str(reactant_state), dict(type=NodeType.state.value))
                    self.graph.add_edge(str(reaction), str(reactant_state), interaction=edge_type.value)
            self.degraded_targets    = [x for x in reaction.degraded_states]
        def _add_reaction_source_state_edges_to_graph(reaction: Reaction, reactants: ReactionTerm):
            """
            Adds reaction source states edges to the graph.

            Note:
                Adds only edges.

            Args:
                reaction: A rxncon reaction.
                reactants: The reactants of the respective reaction.

            Mutates:
                self.graph (DiGraph): The regulatory graph.

            Returns:
                None

            """
            for reactant_state in reactants.states:
                self.graph.add_edge(str(reactant_state), str(reaction), interaction=EdgeInteractionType.source_state.value)

        self.graph.add_node(str(reaction), dict(type=NodeType.reaction.value))

        for reactant_post in reaction.terms_rhs:
            _add_reaction_reactant_to_graph(reaction, reactant_post, EdgeInteractionType.produce)

        for reactant_pre in reaction.terms_lhs:
            _add_reaction_reactant_to_graph(reaction, reactant_pre, EdgeInteractionType.consume)
            _add_reaction_source_state_edges_to_graph(reaction, reactant_pre)

    def add_contingency_information_to_graph(self, contingencies: List[Contingency]) -> None:
        """
        Adds contingency information to the regulatory graph.

        Args:
            contingencies: Contextual constraints on reactions.

        Returns:
            None

        """
        for contingency in contingencies:
            self._add_information_from_effector_to_graph(contingency.effector, contingency.type,
                                                         self._target_name_from_reaction_or_effector(contingency.target)
                                                         )

    def _replace_invalid_chars(self, name: str) -> str:
        """
        Replacing invalid characters.

        Args:
            name: String where the invalid characters should be replaced.

        Returns:
            Valid string for xgmml files.

        """
        return re.sub('[<>]', '', name)

    def _add_node_and_edge(self, name: str, node_type: NodeType, edge_type: EdgeInteractionType, target_name: str) -> None:
        """
        Adds a node or edge to the regulatory graph.

        Args:
            name: Name of the node or edge source.

            node_type: Type of the node e.g. reaction, state, input
            edge_type: Type of the edge e.g. contingency type or source state, production, consumption
            target_name: Name of the edge target.

        Mutates:
            self.graph (DiGraph): The regulatory graph.

        Returns:
            None

        """
        self.graph.add_node(self._replace_invalid_chars(str(name)),
                            type=node_type.value)

        self.graph.add_edge(self._replace_invalid_chars(str(name)), target_name,
                            interaction=edge_type.value)

    def _target_name_from_reaction_or_effector(self, target: Union[Effector, Reaction]) -> str:
        """
        Creates a valid target name, which can be used in xgmml files.

        Args:
            target: Either a reaction of a effector.

        Returns:
            A valid name of a reaction or effector.

        Raises:
            AssertionError if the target is neither of type Reaction or Effector.

        """
        if isinstance(target, Reaction):
            return self._replace_invalid_chars(str(target))
        elif isinstance(target, Effector):
            return self._replace_invalid_chars(target.name)
        else:
            raise AssertionError

    def _add_information_from_effector_to_graph(self, effector, edge_type, target_name) -> None:
        """
        Adds the effector information of the contingency to the graph.

        Args:
            effector: Modifier part of the contingency
            edge_type: Type of connecting edge. The contingency type.
            target_name: Either a reaction or a boolean node.

        Mutates:
            self.graph (DiGraph): The regulatory graph.

        Returns:
            None

        Raises:
            AssertionError: If the Effector type is not known, or if the object is not of type effector.

        """
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
