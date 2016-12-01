from typing import Union, List
import re
from enum import Enum, unique
from networkx import DiGraph

from rxncon.core.rxncon_system import RxnConSystem
from rxncon.venntastic.sets import Intersection, Union as VennUnion, Complement, ValueSet, UniversalSet, Set as VennSet
from rxncon.core.effector import Effector, AndEffector, OrEffector, NotEffector, StateEffector
from rxncon.core.reaction import ReactionTerm, Reaction
from rxncon.core.contingency import Contingency, ContingencyType
from rxncon.core.state import State

from rxncon.core.rxncon_system import RxnConSystem, Reaction, ReactionTerm, Contingency
from rxncon.core.effector import StateEffector, AndEffector, NotEffector, OrEffector, Effector
from rxncon.venntastic.sets import Intersection
from copy import deepcopy


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

edge_type_mapping = {ContingencyType.requirement: EdgeInteractionType.required,
                     ContingencyType.inhibition: EdgeInteractionType.inhibition,
                     ContingencyType.positive: EdgeInteractionType.positive,
                     ContingencyType.negative: EdgeInteractionType.negative}

class BooleanMemory():

    def __init__(self, state, cont_type: ContingencyType, boolean_target: str):

        self.state = state
        self.boolean_target = boolean_target
        self.cont_type = cont_type

    def __hash__(self):
        return hash(str(self))

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "state: {0} boolean_target: {1}".format(str(self.state), self.boolean_target)


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
            if not reaction.degraded_components:
                self.add_reaction_information_to_graph(reaction)
                self.add_contingency_information_to_graph(self.rxncon_system.contingencies_for_reaction(reaction))
            else:
                self.add_degradation_reaction_information_to_graph(reaction, self.rxncon_system.contingencies_for_reaction(reaction))
        return self.graph


    def _add_node(self, id: str, label: str , type: NodeType):
        self.graph.add_node(self._replace_invalid_chars(id), dict(label=self._replace_invalid_chars(label), type=type.value))

    def _add_edge(self, source: str, target: str, interaction: EdgeInteractionType):
        self.graph.add_edge(self._replace_invalid_chars(source), self._replace_invalid_chars(target), interaction=interaction.value)

    def add_degradation_reaction_information_to_graph(self, reaction, contingencies):


            def parse_effector(eff: Effector, cont_type: ContingencyType, boolean_AND_target=None) -> VennSet:
                if isinstance(eff, StateEffector):
                    return ValueSet(BooleanMemory(eff.expr.to_non_structured_state(), cont_type, boolean_AND_target))
                elif isinstance(eff, NotEffector):
                    return Complement(parse_effector(eff.expr, cont_type))
                elif isinstance(eff, OrEffector):
                    return VennUnion(parse_effector(eff.left_expr, cont_type), parse_effector(eff.right_expr, cont_type))
                elif isinstance(eff, AndEffector):
                    if isinstance(eff.left_expr, StateEffector) and isinstance(eff.right_expr, StateEffector):
                        return Intersection(parse_effector(eff.left_expr, cont_type, eff.name), parse_effector(eff.right_expr, cont_type, eff.name))
                    else:
                        return Intersection(parse_effector(eff.left_expr, cont_type), parse_effector(eff.right_expr, cont_type))
                else:
                    raise AssertionError

            def parse_venn(venn: ValueSet):
                if isinstance(venn, ValueSet):
                    return _update_contingency_information(venn)
                #elif isinstance(venn, Complement):
                #    parse_venn(venn.expr)
                #elif isinstance(venn, VennUnion):
                #    parse_venn(venn.exprs[0])
                #    parse_venn(venn.exprs[1])
                elif isinstance(venn, Intersection):
                    for expr in venn.exprs:
                        parse_venn(expr)
                    return
                else:
                    raise AssertionError

            def _interaction_state_degradation(state: State, reaction_id: str):
                self._add_edge(source=reaction_id, target=str(state), interaction=EdgeInteractionType.consume)
                boolean_node_id = '{0}_AND_{1}'.format(reaction_id, str(state))
                self._add_node(id=boolean_node_id, label=' ', type=NodeType.AND)

                self._add_edge(source=reaction_id, target=boolean_node_id, interaction=EdgeInteractionType.AND)
                self._add_edge(source=str(state), target=boolean_node_id, interaction=EdgeInteractionType.AND)

                produced_neutral_states = [neutral_state for neutral_state in state.neutral_states
                                           if not any(component in reaction.degraded_components
                                                      for component in neutral_state.components)]
                for neutral_state in produced_neutral_states:
                    self._add_edge(source=boolean_node_id, target=str(neutral_state), interaction=EdgeInteractionType.produce)

            def _update_no_contingency_case():
                nonlocal reaction_id
                degraded_states = [x for degraded_component in reaction.degraded_components
                                   for x in self.rxncon_system.states_for_component(degraded_component)]
                for state in degraded_states:
                    if len(state.components) > 1:
                        _interaction_state_degradation(state, reaction_id)
                    else:
                        self._add_edge(source=str(reaction), target=str(state), interaction=EdgeInteractionType.consume)

            def _add_contingency_information(state, boolean_target: str, cont_type: EdgeInteractionType, reaction_id):
                boolean_target = self._replace_invalid_chars(boolean_target)
                self._add_node(boolean_target, label=boolean_target, type=NodeType.AND)

                self._add_edge(source=str(state), target=boolean_target, interaction=EdgeInteractionType.AND)
                self._add_edge(source=boolean_target, target=reaction_id, interaction=edge_type_mapping[cont_type])

            def _update_contingency_information(valueset: ValueSet):
                nonlocal reaction_id

                state = valueset.value.state
                cont_type = valueset.value.cont_type
                boolean_target = valueset.value.boolean_target

                if boolean_target:
                    _add_contingency_information(state, boolean_target, cont_type, reaction_id)
                elif not boolean_target:
                    self._add_edge(source=str(state), target=reaction_id, interaction=edge_type_mapping[cont_type])

                for degraded_component in reaction.degraded_components:
                    if degraded_component in state.components:
                        if len(state.components) > 1:
                            _interaction_state_degradation(state, reaction_id)
                        else:
                            self._add_edge(source=reaction_id, target=str(state), interaction=EdgeInteractionType.consume)

            # First case: reaction with a non-trivial contingency should degrade only the states appearing
            # in the contingency that are connected to the degraded component.
            if contingencies:
                cont = Intersection(*(parse_effector(contingency.effector, contingency.type) for contingency in contingencies)).to_simplified_set()
                for index, effector in enumerate(cont.to_dnf_list()):
                    reaction_id = "{0}#{1}".format(str(reaction), str(index))
                    self._add_node(reaction_id, label=str(reaction), type=NodeType.reaction)

                    parse_venn(effector)

            # Second case: reaction with a trivial contingency should degrade all states for the degraded component.
            else:
                self._add_node(id=str(reaction), label=str(reaction), type=NodeType.reaction)
                reaction_id = str(reaction)
                _update_no_contingency_case()


    def add_degradation_reaction_information_to_graph(self, reaction, contingencies):

            def parse_effector(eff: Effector, cont_type: ContingencyType, boolean_AND_target=None) -> VennSet:
                if isinstance(eff, StateEffector):
                    return ValueSet(BooleanMemory(eff.expr.to_non_structured_state(), cont_type, boolean_AND_target))
                elif isinstance(eff, NotEffector):
                    return Complement(parse_effector(eff.expr, cont_type))
                elif isinstance(eff, OrEffector):
                    return VennUnion(parse_effector(eff.left_expr, cont_type), parse_effector(eff.right_expr, cont_type))
                elif isinstance(eff, AndEffector):
                    if isinstance(eff.left_expr, StateEffector) and isinstance(eff.right_expr, StateEffector):
                        return Intersection(parse_effector(eff.left_expr, cont_type, eff.name), parse_effector(eff.right_expr, cont_type, eff.name))
                    else:
                        return Intersection(parse_effector(eff.left_expr, cont_type), parse_effector(eff.right_expr, cont_type))
                else:
                    raise AssertionError

            def parse_venn(venn: ValueSet):
                if isinstance(venn, ValueSet):
                    return _update_contingency_information(venn)
                #elif isinstance(venn, Complement):
                #    parse_venn(venn.expr)
                #elif isinstance(venn, VennUnion):
                #    parse_venn(venn.exprs[0])
                #    parse_venn(venn.exprs[1])
                elif isinstance(venn, Intersection):
                    for expr in venn.exprs:
                        parse_venn(expr)
                    return
                else:
                    raise AssertionError


            def _interaction_state_degradation(state: State, reaction_id: str):
                self._add_edge(source=reaction_id, target=str(state), interaction=EdgeInteractionType.consume)
                boolean_node_id = '{0}_AND_{1}'.format(reaction_id, str(state))
                self._add_node(id=boolean_node_id, label=' ', type=NodeType.AND)

                self._add_edge(source=reaction_id, target=boolean_node_id, interaction=EdgeInteractionType.AND)
                self._add_edge(source=str(state), target=boolean_node_id, interaction=EdgeInteractionType.AND)

                produced_neutral_states = [neutral_state for neutral_state in state.neutral_states
                                           if not any(component in reaction.degraded_components
                                                      for component in neutral_state.components)]
                for neutral_state in produced_neutral_states:
                    self._add_edge(source=boolean_node_id, target=str(neutral_state), interaction=EdgeInteractionType.produce)

            def _update_no_contingency_case():
                nonlocal reaction_id
                degraded_states = [x for degraded_component in reaction.degraded_components
                                   for x in self.rxncon_system.states_for_component(degraded_component)]
                for state in degraded_states:
                    if len(state.components) > 1:
                        _interaction_state_degradation(state, reaction_id)
                    else:
                        self._add_edge(source=str(reaction), target=str(state), interaction=EdgeInteractionType.consume)

            def _add_contingency_information(state, boolean_target: str, cont_type: EdgeInteractionType, reaction_id):
                boolean_target = self._replace_invalid_chars(boolean_target)
                self._add_node(boolean_target, label=boolean_target, type=NodeType.AND)

                self._add_edge(source=str(state), target=boolean_target, interaction=EdgeInteractionType.AND)
                self._add_edge(source=boolean_target, target=reaction_id, interaction=edge_type_mapping[cont_type])

            def _update_contingency_information(valueset: ValueSet):
                nonlocal reaction_id

                state = valueset.value.state
                cont_type = valueset.value.cont_type
                boolean_target = valueset.value.boolean_target

                if boolean_target:
                    _add_contingency_information(state, boolean_target, cont_type, reaction_id)
                elif not boolean_target:
                    self._add_edge(source=str(state), target=reaction_id, interaction=edge_type_mapping[cont_type])

                for degraded_component in reaction.degraded_components:
                    if degraded_component in state.components:
                        if len(state.components) > 1:
                            _interaction_state_degradation(state, reaction_id)
                        else:
                            self._add_edge(source=reaction_id, target=str(state), interaction=EdgeInteractionType.consume)

            # First case: reaction with a non-trivial contingency should degrade only the states appearing
            # in the contingency that are connected to the degraded component.
            if contingencies:
                cont = Intersection(*(parse_effector(contingency.effector, contingency.type) for contingency in contingencies)).to_simplified_set()
                for index, effector in enumerate(cont.to_dnf_list()):
                    reaction_id = "{0}#{1}".format(str(reaction), str(index))
                    self._add_node(reaction_id, label=str(reaction), type=NodeType.reaction)

                    parse_venn(effector)

            # Second case: reaction with a trivial contingency should degrade all states for the degraded component.
            else:
                self._add_node(id=str(reaction), label=str(reaction), type=NodeType.reaction)
                reaction_id = str(reaction)
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

            for reactant_state in reactants.states:
                self._add_node(id=str(reactant_state), type=NodeType.state, label=str(reactant_state))
                self._add_edge(source=str(reaction), target=str(reactant_state), interaction=edge_type)

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
                self._add_edge(source=str(reactant_state), target=str(reaction), interaction=EdgeInteractionType.source_state)

        self._add_node(id=str(reaction), type=NodeType.reaction, label=str(reaction))

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
            name = str(effector.expr.to_non_structured_state())
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
