from typing import Union, List, Optional, Tuple
import re
from enum import Enum, unique
from networkx import DiGraph

from rxncon.venntastic.sets import Union as VennUnion, Complement, ValueSet, Set as VennSet, Intersection
from rxncon.core.contingency import ContingencyType
from rxncon.core.state import State, InteractionState
from rxncon.core.reaction import OutputReaction
from rxncon.core.spec import Spec
from rxncon.core.rxncon_system import RxnConSystem, Reaction, ReactionTerm, Contingency
from rxncon.core.effector import StateEffector, AndEffector, NotEffector, OrEffector, Effector


@unique
class NodeType(Enum):
    reaction  = 'reaction'
    state     = 'state'
    component = 'component'
    AND       = "boolean_and"
    OR        = "boolean_or"
    NOT       = "boolean_not"
    boolean   = "boolean"
    input     = 'input'
    output    = 'output'


@unique
class EdgeInteractionType(Enum):
    produce        = 'produce'
    consume        = 'consume'
    synthesis      = 'synthesis'
    degrade        = 'degrade'
    maybe_degraded = 'maybe_degraded'
    required       = '!'
    no_effect      = '0'
    unknown_effect = '?'
    inhibited      = "x"
    positive       = 'k+'
    negative       = 'k-'
    AND            = 'AND'
    OR             = 'OR'
    NOT            = 'NOT'
    source_state   = 'ss'
    input_state    = 'is'
    mutually_exclusive = 'mutually_exclusive'
    modify = "modifiy"

edge_type_mapping = {ContingencyType.requirement: EdgeInteractionType.required,
                     ContingencyType.inhibition: EdgeInteractionType.inhibited,
                     ContingencyType.positive: EdgeInteractionType.positive,
                     ContingencyType.negative: EdgeInteractionType.negative,
                     ContingencyType.no_effect: EdgeInteractionType.no_effect,
                     ContingencyType.unknown: EdgeInteractionType.unknown_effect}

class SpeciesReactionGraph:
    """
    Definition of the regulatory Graph.

    Args:
        rxncon_system: The rxncon system.
    """
    def __init__(self, rxncon_system: RxnConSystem) -> None:
        self.rxncon_system = rxncon_system
        self.species_reaction_graph = DiGraph()

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
        self.add_synthesised_or_degraded_components()
        return self.species_reaction_graph

    def add_synthesised_or_degraded_components(self) -> None:
        """
        Adding components to the graph, which do not belong to any state of the system but get synthesised or degraded.

        Mutates:
            The regulatory graph, by adding reaction nodes and/or component nodes.

        Returns:
            None

        """
        def calc_components_without_states() -> List[Spec]:
            """
            Calculating the components without states:

            Returns:
                List of components which do not belong to at least one state in the system.

            """
            return [comp for comp in self.rxncon_system.components() if not self.rxncon_system.states_for_component(comp)]

        def add_synthesised_components_and_reactions() -> None:
            """
            Adding components to the graph, which do not belong to at least one state but gets synthesised.

            Mutates:
                The regulatory graph, by adding component nodes and edges from reactions to component nodes.
                components_by_reactions: storing components which are synthesised or degraded but do no belong to any
                state.

            Returns:
                None

            """
            nonlocal components_by_reactions
            for synthesised_component in reaction.synthesised_components:
                if synthesised_component in components_without_states:
                    self._add_node(id=str(synthesised_component), label=str(synthesised_component), type=NodeType.component)
                    self._add_edge(source=str(reaction), target=str(synthesised_component), interaction=EdgeInteractionType.synthesis)
                    components_by_reactions.append(synthesised_component)

        def get_degraded_components_and_reactions() -> None:
            """
            Adding components to the graph, which do not belong to at least one state but gets degraded.

            Mutates:
                The regulatory graph, by adding component nodes and edges from reactions to components.
                components_by_reactions: storing components which are synthesised or degraded but do no belong to any
                state.

            Returns:
                None

            """
            nonlocal components_by_reactions
            for degraded_component in reaction.degraded_components:
                if degraded_component in components_without_states:
                    self._add_node(id=str(degraded_component), label=str(degraded_component), type=NodeType.component)
                    self._add_edge(source=str(reaction), target=str(degraded_component), interaction=EdgeInteractionType.degrade)
                    components_by_reactions.append(degraded_component)

        def connect_components_and_reactions() -> None:
            """
            Connecting components, which do not belong to any state to reactions.

            Note: If a component, which do not belong to any state, is synthesised or degraded the component is
                  mentioned explicitly in the regulatory graph. In this case the component has to be connected to all
                  reactions it is involved in.

            Mutates:
                The regulatory graph, by adding edges.
            Returns:
                None

            """
            for component in components_by_reactions:
                for reaction in self.rxncon_system.reactions:
                    if component in reaction.components_lhs and not component in reaction.components_rhs:
                        self._add_edge(source=str(component), target=str(reaction), interaction=EdgeInteractionType.source_state)
                    elif component in reaction.components_lhs:
                        self._add_edge(source=str(component), target=str(reaction), interaction=EdgeInteractionType.input_state)

        components_by_reactions = []  # type: List[Spec]
        components_without_states = calc_components_without_states()

        for reaction in self.rxncon_system.reactions:
            add_synthesised_components_and_reactions()
            get_degraded_components_and_reactions()

        connect_components_and_reactions()

    def add_degradation_reaction_information_to_graph(self, reaction: Reaction, contingencies: List[Contingency]) -> None:
        """
        Adding degradation information to the graph.

        Args:
            reaction: rxncon reaction.
            contingencies: list of rxncon contingencies.

        Mutates:
            The regulatory Graph.

        Returns:
            None

        """

        def _effector_to_vennset(eff: Effector, con_type: Optional[ContingencyType]=None) -> VennSet:
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

        def _add_interaction_state_for_degradation(state: State, reaction: Reaction, edge_type: EdgeInteractionType) ->None:
            """
            Adds interaction state for degradation.

            Note:
                If we degrade an interaction state, we degrade one component of the state and release the other unbound.
                For visualisation purpose we add an additional boolean node reflecting this interpretation.

            Args:
                state: interaction state
                reaction: reaction of the degradation reaction.
                edge_type: type of the edge. This can be DEGRADED or MAYBE_DEGRADED

            Returns:
                None

            """
            assert isinstance(state, InteractionState)
            assert reaction.degraded_components
            assert edge_type in (EdgeInteractionType.degrade, EdgeInteractionType.maybe_degraded)

            self._add_edge(source=str(reaction), target=str(state), interaction=edge_type)

            if not state.is_homodimer:
                boolean_node_id = '{0}_ON_{1}'.format(str(reaction), str(state))
                self._add_node(id=boolean_node_id, label=' ', type=NodeType.AND)

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
                if isinstance(state, InteractionState):
                    _add_interaction_state_for_degradation(state, reaction, EdgeInteractionType.degrade)
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
            nonlocal reaction
            degraded_states = [x for degraded_component in reaction.degraded_components
                               for x in self.rxncon_system.states_for_component(degraded_component)
                               if x not in positive_states and x not in negative_states]

            for state in degraded_states:
                if not any(positive_state for positive_state in positive_states if state.is_mutually_exclusive_with(positive_state)):
                    if isinstance(state, InteractionState):
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
                    complements = self.rxncon_system.complement_states_for_component(degraded_component, state)
                    # If we have one unique complement of a state we know what gets degraded.
                    if len(complements) == 1:
                        self._add_edge(source=str(reaction), target=str(complements[0]), interaction=EdgeInteractionType.degrade)
                    # If we have more than one complement of a state. We say that its a possible degradation.
                    else:
                        for state in complements:
                            self._add_edge(source=str(reaction), target=str(state), interaction=EdgeInteractionType.maybe_degraded)

        def _update_contingency_information(value_set: ValueSet[State]) -> None:
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
                    if isinstance(state, InteractionState):
                        _add_interaction_state_for_degradation(state, reaction, EdgeInteractionType.degrade)
                    else:
                        self._add_edge(source=str(reaction), target=str(state), interaction=EdgeInteractionType.degrade)

        def _update_contingency_information_for_complement(complement_value: VennSet[State]) -> None:
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

        def _get_positive_and_negative_states(nested_list: List[VennSet[State]], dnf_of_cont: List[List[VennSet[State]]])\
                -> Tuple[List[Complement[State]], List[ValueSet[State]]]:
            """
            Calculating a list of negative states (complements) and positive states (not complements).

            Args:
                nested_list: List if Complements and ValueSets.
                dnf_of_cont: disjunctive normal form of the contingency list belonging to a certain reaction.

            Returns:
                A list of complements (negative_value_set) and positive_value_set.

            """
            negative_value_set = []  # type: List[Complement[State]]
            positive_value_set = []  # type: List[ValueSet[State]]

            for value in nested_list:
                if all(value in nested_list for nested_list in dnf_of_cont):
                    if isinstance(value, Complement):
                        negative_value_set.append(value)
                    elif isinstance(value, ValueSet):
                        positive_value_set.append(value)
                    else:
                        raise AssertionError

            return negative_value_set, positive_value_set

        # First case: reaction with a non-trivial contingency should degrade only the states appearing
        # in the contingency that are connected to the degraded component.
        self._add_node(id=str(reaction), label=str(reaction), type=NodeType.reaction)
        if contingencies:
            self.add_contingency_information_to_graph(contingencies)

            cont = Intersection(*(_effector_to_vennset(contingency.effector, contingency.contingency_type) for contingency in contingencies)).to_simplified_set()
            dnf_of_cont = cont.to_dnf_nested_list()

            for index, nested_list in enumerate(dnf_of_cont):
                negative_value_sets, positive_value_sets = _get_positive_and_negative_states(nested_list, dnf_of_cont)

                for pos_value_set in positive_value_sets:
                    _update_contingency_information(pos_value_set)

                for neg_value_set in negative_value_sets:
                    _update_contingency_information_for_complement(neg_value_set)

                positive_states = [pos_value_set.value for pos_value_set in positive_value_sets]

                negative_states = []
                for neg_value_set in negative_value_sets:
                    assert isinstance(neg_value_set.expr, ValueSet)
                    negative_states.append(neg_value_set.expr.value)

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

        def _add_reactant_states() -> None:
            """
            Adding the states of reactants.

            Returns:
                None

            """
            for reactant_post in reaction.terms_rhs:
                if reaction.synthesised_states:
                    _add_reaction_reactant_to_graph(reaction, reactant_post, EdgeInteractionType.synthesis)
                else:
                    _add_reaction_reactant_to_graph(reaction, reactant_post, EdgeInteractionType.produce)

            for reactant_pre in reaction.terms_lhs:
                _add_reaction_reactant_to_graph(reaction, reactant_pre, EdgeInteractionType.consume)
                _add_reaction_source_state_edges_to_graph(reaction, reactant_pre)
        if isinstance(reaction, OutputReaction):
            self._add_node(id=str(reaction), type=NodeType.output, label=str(reaction))
        else:
            self._add_node(id=str(reaction), type=NodeType.reaction, label=str(reaction))
            _add_reactant_states()



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
                assert target.name is not None
                return self._replace_invalid_chars(target.name)
            else:
                raise AssertionError

        def _add_information_from_effector_to_graph(effector: Effector, edge_type: EdgeInteractionType, target_name: str) -> None:
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

                for expr in effector.exprs:
                    _add_information_from_effector_to_graph(expr, EdgeInteractionType.OR, target_name)
            elif isinstance(effector, AndEffector):
                if effector.name is not None:
                    add_node_and_edge(effector.name, NodeType.AND, edge_type, target_name)
                    target_name = _target_name_from_reaction_or_effector(effector)

                for expr in effector.exprs:
                    _add_information_from_effector_to_graph(expr, EdgeInteractionType.AND, target_name)
            else:
                raise AssertionError

        for contingency in contingencies:
            _add_information_from_effector_to_graph(contingency.effector, edge_type_mapping[contingency.contingency_type],
                                                    _target_name_from_reaction_or_effector(contingency.reaction))

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
        self.species_reaction_graph.add_node(self._replace_invalid_chars(id), label=self._replace_invalid_chars(label), type=type.value)

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
        if not self.species_reaction_graph.has_edge(self._replace_invalid_chars(source), self._replace_invalid_chars(target)):
            self.species_reaction_graph.add_edge(self._replace_invalid_chars(source), self._replace_invalid_chars(target), interaction=interaction.value)

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



class RegulatoryGraph:
    """
    Definition of the regulatory Graph.

    Args:
        rxncon_system: The rxncon system.
    """
    def __init__(self, rxncon_system: RxnConSystem) -> None:
        self.rxncon_system = rxncon_system
        self.regulatory_graph = DiGraph()

    def to_graph(self) -> DiGraph:
        """
        Creates the regulatory graph.

        Returns:
            The regulatory graph as networkx DiGraph.

        """

        for reaction in self.rxncon_system.reactions:
            # if not reaction.degraded_components:
            self.add_reaction_information_to_graph(reaction)
            self.add_contingency_information_to_graph(self.rxncon_system.contingencies_for_reaction(reaction))
            # else:
            #     self.add_degradation_reaction_information_to_graph(reaction, self.rxncon_system.contingencies_for_reaction(reaction))
        self.add_synthesised_or_degraded_components()
        return self.regulatory_graph

    def add_synthesised_or_degraded_components(self) -> None:
        """
        Adding components to the graph, which do not belong to any state of the system but get synthesised or degraded.

        Mutates:
            The regulatory graph, by adding reaction nodes and/or component nodes.

        Returns:
            None

        """

        ######## Start no use in regulatory graph
        def calc_components_without_states() -> List[Spec]:
            """
            Calculating the components without states:

            Returns:
                List of components which do not belong to at least one state in the system.

            """
            return [comp for comp in self.rxncon_system.components() if not self.rxncon_system.states_for_component(comp)]

        ######## END no use in regulatory graph
        def add_synthesised_components_and_reactions() -> None:
            """
            Adding components to the graph, which do not belong to at least one state but gets synthesised.

            Mutates:
                The regulatory graph, by adding component nodes and edges from reactions to component nodes.
                components_by_reactions: storing components which are synthesised or degraded but do no belong to any
                state.

            Returns:
                None

            """
            nonlocal components_by_reactions
            for synthesised_component in reaction.synthesised_components:
                # if synthesised_component in components_without_states:
                self._add_node(id=str(synthesised_component), label=str(synthesised_component), type=NodeType.component)
                self._add_edge(source=str(reaction), target=str(synthesised_component), interaction=EdgeInteractionType.synthesis)
                components_by_reactions.append(synthesised_component)

        def get_degraded_components_and_reactions() -> None:
            """
            Adding components to the graph, which do not belong to at least one state but gets degraded.

            Mutates:
                The regulatory graph, by adding component nodes and edges from reactions to components.
                components_by_reactions: storing components which are synthesised or degraded but do no belong to any
                state.

            Returns:
                None

            """
            nonlocal components_by_reactions
            for degraded_component in reaction.degraded_components:
                # if degraded_component in components_without_states:
                self._add_node(id=str(degraded_component), label=str(degraded_component), type=NodeType.component)
                self._add_edge(source=str(reaction), target=str(degraded_component), interaction=EdgeInteractionType.degrade)
                components_by_reactions.append(degraded_component)

        def component_in_conts(comp, reaction):

            cont = Intersection(*(x.to_venn_set() for x
                                in self.rxncon_system.contingencies_for_reaction(reaction)))

            for state in cont.values:
                if comp in state.to_non_structured().components:
                    return True
            return False

        def get_component_for_neutral_states(reaction):
            nonlocal components_by_reactions
            for state in reaction.consumed_states:
                # if len(state.components) == 1 and state.components[0] == component and state.is_neutral:
                if state.is_neutral and state.components[0] in components_by_reactions:
                    assert len(state.components) == 1
                    if not component_in_conts(state.components[0], reaction):
                        self._add_edge(source=str(state.components[0]), target=str(reaction),
                                   interaction=EdgeInteractionType.input_state)
                    # mutually exclusivity test:
                    complements = self.rxncon_system.complement_states(state)
                    if len(complements) > 1:
                        for state in complements:
                            if state not in reaction.produced_states:
                                self._add_edge(source=str(state), target=str(reaction),
                                               interaction=EdgeInteractionType.mutually_exclusive)

        def _add_modifier_components(reaction):
            modifiers = reaction.modifier_components
            for mod in modifiers:
                if mod in components_by_reactions:
                    self._add_node(id=str(mod), type=NodeType.component, label=str(mod))
                    self._add_edge(source=str(mod), target=str(reaction), interaction=EdgeInteractionType.modify)




        def connect_components_and_reactions() -> None:
            """
            Connecting components, which do not belong to any state to reactions.

            Note: If a component, which do not belong to any state, is synthesised or degraded the component is
                  mentioned explicitly in the regulatory graph. In this case the component has to be connected to all
                  reactions it is involved in.

            Mutates:
                The regulatory graph, by adding edges.
            Returns:
                None

            """
            for component in components_by_reactions:
                for reaction in self.rxncon_system.reactions:
                    if component in reaction.components_lhs and not component in reaction.components_rhs:
                        self._add_edge(source=str(component), target=str(reaction), interaction=EdgeInteractionType.source_state)
                    elif component in reaction.components_lhs :
                        if any(state.is_neutral and component in state.components and not component_in_conts(component, reaction) for state in reaction.consumed_states):
                            self._add_edge(source=str(component), target=str(reaction),
                                           interaction=EdgeInteractionType.input_state)

        components_by_reactions = []  # type: List[Spec]
        components_without_states = calc_components_without_states()

        for reaction in self.rxncon_system.reactions:
            add_synthesised_components_and_reactions()
            get_degraded_components_and_reactions()
        for reaction in self.rxncon_system.reactions:
            get_component_for_neutral_states(reaction)
            _add_modifier_components(reaction)
        connect_components_and_reactions()


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

        def _add_reactant_states() -> None:
            """
            Adding the states of reactants.

            Returns:
                None

            """
            for reactant_post in reaction.terms_rhs:
                if reaction.synthesised_states:
                    if reactant_post.states and not reactant_post.states[0].is_neutral:
                        _add_reaction_reactant_to_graph(reaction, reactant_post, EdgeInteractionType.synthesis)
                else:

                    if reactant_post.states and not reactant_post.states[0].is_neutral:
                        _add_reaction_reactant_to_graph(reaction, reactant_post, EdgeInteractionType.produce)

            for reactant_pre in reaction.terms_lhs:
                if reactant_pre.states and not reactant_pre.states[0].is_neutral:
                    _add_reaction_reactant_to_graph(reaction, reactant_pre, EdgeInteractionType.consume)
                    _add_reaction_source_state_edges_to_graph(reaction, reactant_pre)

        if isinstance(reaction, OutputReaction):
            self._add_node(id=str(reaction), type=NodeType.output, label=str(reaction))
        else:
            self._add_node(id=str(reaction), type=NodeType.reaction, label=str(reaction))
            _add_reactant_states()



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
                assert target.name is not None
                return self._replace_invalid_chars(target.name)
            else:
                raise AssertionError

        def check_state_list_for_neutrals(states):

            output = []
            for state in states:
                assert(isinstance(state, State))
                if state.is_neutral:
                    output.append(state)
            return output

        def _add_information_from_effector_to_graph(effector: Effector, edge_type: EdgeInteractionType, target_name: str) -> None:
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
                neutrals_in_cont = check_state_list_for_neutrals(effector.states)

                if re.match("^\[[/w]*\]?", name):
                    add_node_and_edge(name, NodeType.input, edge_type, target_name)

                elif neutrals_in_cont:
                    add_node_and_edge('<'+name+'>', NodeType.NOT, edge_type, target_name)
                    # check if single

                    complementary_states = self.rxncon_system.complement_states(effector.expr)
                    if len(complementary_states) > 1:
                        add_node_and_edge('not <' + name + '>', NodeType.NOT, edge_type=EdgeInteractionType.NOT,
                                          target_name='<' + name + '>')
                        for state in complementary_states:
                            add_node_and_edge(name=str(state.to_non_structured()), node_type=NodeType.state,
                                              edge_type= EdgeInteractionType.OR, target_name='not <'+name+'>')
                    else:
                        add_node_and_edge(name=str(complementary_states[0].to_non_structured()), node_type=NodeType.state,
                                          edge_type=EdgeInteractionType.NOT, target_name='<' + name + '>')


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

                for expr in effector.exprs:
                    _add_information_from_effector_to_graph(expr, EdgeInteractionType.OR, target_name)
            elif isinstance(effector, AndEffector):
                if effector.name is not None:
                    add_node_and_edge(effector.name, NodeType.AND, edge_type, target_name)
                    target_name = _target_name_from_reaction_or_effector(effector)

                for expr in effector.exprs:
                    _add_information_from_effector_to_graph(expr, EdgeInteractionType.AND, target_name)
            else:
                raise AssertionError

        for contingency in contingencies:
            _add_information_from_effector_to_graph(contingency.effector, edge_type_mapping[contingency.contingency_type],
                                                    _target_name_from_reaction_or_effector(contingency.reaction))

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
        self.regulatory_graph.add_node(n=self._replace_invalid_chars(id), attr_dict=None,
                                       label=self._replace_invalid_chars(label), type=type.value)

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



