from typing import Union, List
import re
from enum import Enum, unique
from networkx import DiGraph

from rxncon.venntastic.sets import Union as VennUnion, Complement, ValueSet, Set as VennSet
from rxncon.core.contingency import ContingencyType
from rxncon.core.state import State


from rxncon.core.rxncon_system import RxnConSystem, Reaction, ReactionTerm, Contingency
from rxncon.core.effector import StateEffector, AndEffector, NotEffector, OrEffector, Effector
from rxncon.venntastic.sets import Intersection


@unique
class NodeType(Enum):
    reaction = 'reaction'
    state    = 'state'
    AND      = "boolean_and"
    OR       = "boolean_or"
    NOT      = "boolean_not"
    boolean  = "boolean"
    input    = 'input'
    output   = 'output'

@unique
class EdgeInteractionType(Enum):
    produce       = 'produce'
    consume       = 'consume'
    synthesis     = 'synthesis'
    degrade       = 'degrade'
    maybe_degraded = 'maybe_degraded'
    required      = '!'
    inhibited    = "x"
    positive      = 'k+'
    negative      = 'k-'
    AND           = 'AND'
    OR            = 'OR'
    NOT           = 'NOT'
    source_state  = 'ss'

edge_type_mapping = {ContingencyType.requirement: EdgeInteractionType.required,
                     ContingencyType.inhibition: EdgeInteractionType.inhibited,
                     ContingencyType.positive: EdgeInteractionType.positive,
                     ContingencyType.negative: EdgeInteractionType.negative}


class RegulatoryGraph:
    """
    Definition of the regulatory Graph.

    Args:
        rxncon_system: The rxncon system.
    """
    def __init__(self,rxncon_system: RxnConSystem, potential_degradation: bool = False) -> None:
        self.rxncon_system = rxncon_system
        self.potential_degradation = potential_degradation
        self.regulatory_graph = DiGraph()

    def to_graph(self) -> DiGraph:
        """
        Creates the regulatory graph.

        Returns:
            The regulatory graph as networkx DiGraph.

        """
        for reaction in self.rxncon_system.reactions:
            if not reaction.degraded_components:
                self.add_reaction_information_to_graph(reaction)
                self.add_contingency_information_to_graph(self.rxncon_system.contingencies_for_reaction(reaction))
            else:
                self.add_degradation_reaction_information_to_graph(reaction, self.rxncon_system.contingencies_for_reaction(reaction))
        return self.regulatory_graph

    def add_degradation_reaction_information_to_graph(self, reaction, contingencies) -> None:
        """
        Adding degradation information to the graph.

        Args:
            reaction: rxncon reaction.
            contingencies: rxncon contingency.

        Mutates:
            The regulatory Graph.

        Returns:
            None

        """

        def _effector_to_vennset(eff: Effector, con_type: Union[ContingencyType, None]=None) -> VennSet:
            """
            Preprocessing effector. Save information in leafs.

            Note:
                We need the contingency information later on in the process. Since, this information is lost during the
                calculation of the disjunctive normal form we have to add it beforehand. For this we say an inhibiting
                contingency is the complement of everything afterwards.

                The contingency type is only needed at the beginning of the recursion step.

            Args:
                eff: rxncon Effector.
                cont_type: Contingency type.

            Returns:
                Return VennSet

            Raises:
                AssertionError if the eff is not an effector.

            """

            if isinstance(eff, StateEffector):
                if con_type is ContingencyType.inhibition:
                    return Complement(ValueSet(eff.expr.to_non_structured()))
                return ValueSet(eff.expr.to_non_structured())
            elif isinstance(eff, NotEffector):
                if con_type is ContingencyType.inhibition:
                    return _effector_to_vennset(eff.expr)
                return Complement(_effector_to_vennset(eff.expr))
            elif isinstance(eff, OrEffector):
                if con_type is ContingencyType.inhibition:
                    return Complement(VennUnion(*[_effector_to_vennset(expr) for expr in eff.exprs]))
                return VennUnion(*[_effector_to_vennset(expr) for expr in eff.exprs])
            elif isinstance(eff, AndEffector):
                if con_type is ContingencyType.inhibition:
                    return Complement(Intersection(*[_effector_to_vennset(expr) for expr in eff.exprs]))
                return Intersection(*[_effector_to_vennset(expr) for expr in eff.exprs])
            else:
                raise AssertionError

        def _add_interaction_state_for_degradation(state: State, reaction: Reaction, edge_type: EdgeInteractionType = EdgeInteractionType.degrade) ->None:
            """
            Adds interaction state for degradation.

            Note:
                If we degrade an interaction state, we degrade one component of the state and release the other unbound.
                For visualisation purpose we add an additional boolean node reflecting this interpretation.

            Args:
                state: interaction state
                reaction: reaction ID of the degradation reaction.
                edge_type: type of the edge.

            Returns:
                None

            """
            self._add_edge(source=str(reaction), target=str(state), interaction=edge_type)
            boolean_node_id = '{0}_ON_{1}'.format(str(reaction), str(state))
            self._add_node(id=boolean_node_id, label=' ', type=NodeType.boolean)

            self._add_edge(source=str(reaction), target=boolean_node_id, interaction=EdgeInteractionType.AND)
            self._add_edge(source=str(state), target=boolean_node_id, interaction=EdgeInteractionType.AND)

            produced_neutral_states = [neutral_state for neutral_state in state.neutral_states
                                       if not any(component in reaction.degraded_components
                                                  for component in neutral_state.components)]
            for neutral_state in produced_neutral_states:
                self._add_edge(source=boolean_node_id, target=str(neutral_state), interaction=EdgeInteractionType.produce)

        def _update_no_contingency_case() -> None:
            """
            Updating degradation information if no contingencies are given.

            Mutates:
                regulatory graph.

            Returns:
                None

            """
            nonlocal reaction
            degraded_states = [x for degraded_component in reaction.degraded_components
                               for x in self.rxncon_system.states_for_component(degraded_component)]
            for state in degraded_states:
                if len(state.components) > 1:
                    _add_interaction_state_for_degradation(state, reaction)
                else:
                    self._add_edge(source=str(reaction), target=str(state), interaction=EdgeInteractionType.degrade)

        def _add_possible_degraded_states(positive_states: List[State], negative_states: List[State]) -> None:
            """
            Adding optional degradation.

            Note:
                If there is a required contingency on a degradation only the these contingencies will
                have a degradation edge for sure. We are adding optional degradation edges to all other possible
                (non-mutually exclusive) states.

            Args:
                positive_states: List of not negated value sets
                negative_states: List of negated value sets.

            Returns:
                None

            """

            degraded_states = [x for degraded_component in reaction.degraded_components
                               for x in self.rxncon_system.states_for_component(degraded_component)
                               if x not in positive_states and x not in negative_states]

            for state in degraded_states:
                if not any(positive_state for positive_state in positive_states if state.is_mutually_exclusive_with(positive_state)):
                    if len(state.components) > 1:
                        _add_interaction_state_for_degradation(state, reaction, EdgeInteractionType.maybe_degraded)
                    else:
                        self._add_edge(source=str(reaction), target=str(state), interaction=EdgeInteractionType.maybe_degraded)

        def _add_complement_of_state_for_degradation_reaction(state: State) -> None:
            """
            Adding complement of a state to the regulatory graph.

            Args:
                state: rxncon state.

            Mutates:
                regulatory graph.

            Returns:
                None

            """
            for degraded_component in reaction.degraded_components:
                if degraded_component in state.components:
                    complements = self.rxncon_system.complementary_states_for_component(degraded_component, state)
                    # If we have one unique complement of a state we know what gets degraded.
                    if len(complements) == 1:
                        self._add_edge(source=str(reaction), target=str(complements[0]), interaction=EdgeInteractionType.degrade)
                    # If we have more than one complement of a state. We say that its a possible degradation.
                    else:
                        for state in complements:
                            self._add_edge(source=str(reaction), target=str(state), interaction=EdgeInteractionType.maybe_degraded)

        def _update_contingency_information(value_set: ValueSet) -> None:
            """
            Updating the degradation information with its contingencies.

            Args:
                value_set: VennSet object.

            Mutates:
                regulatory_graph.

            Returns:
                None

            """

            state = value_set.value

            for degraded_component in reaction.degraded_components:
                if degraded_component in state.components:
                    if len(state.components) > 1:
                        _add_interaction_state_for_degradation(state, reaction)
                    else:
                        self._add_edge(source=str(reaction), target=str(state), interaction=EdgeInteractionType.degrade)

        def _update_contingency_information_for_complement(complement_value: Complement):
            """
            Updating the contingency information for complements.

            Note:
                If we have a complement of a value set we have to add the complemented states of these value set to
                regulatory graph.

            Args:
                complement_value: Complement of a value set.

            Returns:
                None

            Raises:
                AssertionError if the the complement_value has more than one entry.

            """
            assert len(complement_value.values) == 1
            state = complement_value.values[0]
            _add_complement_of_state_for_degradation_reaction(state)

        def _get_positive_and_negative_states(nested_list: List[Union[ValueSet, Complement]], dnf_of_cont: List[List[Union[ValueSet, Complement]]]):
            """
            Calculating a list of negative states (complements) and positive states (not complements).

            Args:
                nested_list: List if Complements and ValueSets.
                dnf_of_cont: disjunctive normal form of the contingency list belonging to a certain reaction.

            Returns:
                A list of complements (negative_value_set) and positive_value_set.

            """
            negative_value_set = []  # type: List[Complement]
            positive_value_set = []  # type: List[ValueSet]

            for value in nested_list:
                if all(value in nested_list for nested_list in dnf_of_cont):
                    if isinstance(value, Complement):
                        negative_value_set.append(value)
                    else:
                        positive_value_set.append(value)

            return negative_value_set, positive_value_set

        # First case: reaction with a non-trivial contingency should degrade only the states appearing
        # in the contingency that are connected to the degraded component.
        self._add_node(id=str(reaction), label=str(reaction), type=NodeType.reaction)
        if contingencies:

            self.add_contingency_information_to_graph(contingencies)

            cont = Intersection(*(_effector_to_vennset(contingency.effector, contingency.type) for contingency in contingencies)).to_simplified_set()
            dnf_of_cont = cont.to_dnf_nested_list()

            for index, nested_list in enumerate(dnf_of_cont):

                negative_value_set, positive_value_set = _get_positive_and_negative_states(nested_list, dnf_of_cont)

                [_update_contingency_information(value_set) for value_set in positive_value_set]
                [_update_contingency_information_for_complement(value_set) for value_set in negative_value_set]

                positive_states = [value_set.value for value_set in positive_value_set]
                negative_states = [value_set.expr.value for value_set in negative_value_set]
                _add_possible_degraded_states(positive_states, negative_states)

        # Second case: reaction with a trivial contingency should degrade all states for the degraded component.
        else:
            _update_no_contingency_case()

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

            for reactant_state in reactants.states:
                self._add_node(id=str(reactant_state), type=NodeType.state, label=str(reactant_state))
                self._add_edge(source=str(reaction), target=str(reactant_state), interaction=edge_type)

        def _add_reaction_source_state_edges_to_graph(reaction: Reaction, reactants: ReactionTerm) -> None:
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
                self._add_edge(source=str(reactant_state), target=str(reaction), interaction=EdgeInteractionType.source_state)

        self._add_node(id=str(reaction), type=NodeType.reaction, label=str(reaction))

        for reactant_post in reaction.terms_rhs:
            if reaction.synthesised_states:
                _add_reaction_reactant_to_graph(reaction, reactant_post, EdgeInteractionType.synthesis)
            else:
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
        def _target_name_from_reaction_or_effector(target: Union[Effector, Reaction]) -> str:
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

        def _add_information_from_effector_to_graph(effector, edge_type, target_name) -> None:
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
            def add_node_and_edge(name: str, node_type: NodeType, edge_type: EdgeInteractionType, target_name: str) -> None:
                """
                Adds a node or edge to the regulatory graph.

                Args:
                    name: Name of the node or edge source.

                    node_type: Type of the node e.g. state, input

                    edge_type: Type of the edge e.g. contingency type or source state, production, consumption
                    target_name: Name of the edge target.

                Mutates:
                    self.graph (DiGraph): The regulatory graph.


                Returns:
                    None

                """

                self._add_node(str(name), type=node_type, label=str(name))
                self._add_edge(source=str(name), target=target_name, interaction=edge_type)

            if isinstance(effector, StateEffector):
                name = str(effector.expr.to_non_structured())
                if re.match("^\[[/w]*\]?", name):
                    add_node_and_edge(name, NodeType.input, edge_type, target_name)
                else:
                    add_node_and_edge(name, NodeType.state, edge_type, target_name)

            elif isinstance(effector, NotEffector):
                add_node_and_edge(str(effector.name), NodeType.NOT, edge_type, target_name)
                target_name = _target_name_from_reaction_or_effector(effector)
                _add_information_from_effector_to_graph(effector.expr, EdgeInteractionType.NOT, target_name)

            elif isinstance(effector, OrEffector):
                if effector.name is not None:
                    add_node_and_edge(effector.name, NodeType.OR, edge_type, target_name)
                    target_name = _target_name_from_reaction_or_effector(effector)

                [_add_information_from_effector_to_graph(expr, EdgeInteractionType.OR, target_name) for expr in effector.exprs]

            elif isinstance(effector, AndEffector):
                if effector.name is not None:
                    add_node_and_edge(effector.name, NodeType.AND, edge_type, target_name)
                    target_name = _target_name_from_reaction_or_effector(effector)

                [_add_information_from_effector_to_graph(expr, EdgeInteractionType.AND, target_name) for expr in effector.exprs]
            else:
                raise AssertionError

        for contingency in contingencies:
            _add_information_from_effector_to_graph(contingency.effector, contingency.type,
                                                         _target_name_from_reaction_or_effector(contingency.target))

    def _add_node(self, id: str, label: str, type: NodeType) -> None:
        """
        Adding a node to the graph.

        Args:
            id: Node ID
            label: Node label
            type: Node type

        Mutates:
            The regulatory graph.

        Returns:
            None

        """
        self.regulatory_graph.add_node(self._replace_invalid_chars(id), dict(label=self._replace_invalid_chars(label), type=type.value))

    def _add_edge(self, source: str, target: str, interaction: EdgeInteractionType) -> None:
        """
        Adding a edge to the graph

        Args:
            source: Node ID of the edge source
            target: Node ID of the edge target
            interaction: Type of interaction.

        Mutates:
            The regulatory graph.

        Returns:
            None

        """
        if not self.regulatory_graph.has_edge(self._replace_invalid_chars(source), self._replace_invalid_chars(target)):
            self.regulatory_graph.add_edge(self._replace_invalid_chars(source), self._replace_invalid_chars(target), interaction=interaction.value)

    def _replace_invalid_chars(self, name: str) -> str:
        """
        Replacing invalid characters.

        Args:
            name: String where the invalid characters should be replaced.

        Returns:
            Valid string for xgmml files.

        """
        name = re.sub('[<>]', '', name)
        name = re.sub('[&]', 'AND', name)
        return name



