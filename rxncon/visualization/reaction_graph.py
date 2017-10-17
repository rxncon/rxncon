import logging
import re
from typing import Union, Type
from enum import Enum
from networkx import DiGraph
from rxncon.core.reaction import Reaction, OutputReaction
from rxncon.core.rxncon_system import RxnConSystem
from rxncon.core.spec import Spec, LocusResolution
from rxncon.core.state import ModificationState, InteractionState, SelfInteractionState, EmptyBindingState

logger = logging.getLogger(__name__)


class EdgeWith(Enum):
    internal = '0'
    external = '10'


class EdgeType(Enum):
    interaction = 'interaction'
    modification = 'modification'
    bimodification = 'bimodification'
    synthesis = 'synthesis'
    degradation = 'degradation'
    unknown = 'unknown'


class NodeType(Enum):
    component = '80'
    domain = '40'
    residue = '20'


STATES = {
    InteractionState: EdgeType.interaction,
    EmptyBindingState: EdgeType.interaction,
    SelfInteractionState: EdgeType.interaction,
    ModificationState: EdgeType.modification
}


class ReactionGraph:
    """
    Stores the reaction graph.

    Args:
        reaction_graph: A DiGraph with reaction information.

    """

    def __init__(self, reaction_graph: DiGraph) -> None:
        self.reaction_graph = reaction_graph
        self._validate_graph()

    def _validate_graph(self) -> None:
        """
        Testing if all nodes are connected through the graph.
        Note:
            The reaction_graph is a directed graph. We are checking both directions source -> target by asking for
            neighbors and target -> source by asking for predecessors. If one of both properties are fulfilled the node
            is connected.

        Returns:
            None

        Raises:
            AssertionError if not all nodes are connected.

        """
        assert all(self.reaction_graph.neighbors(node) != [] or self.reaction_graph.predecessors(node) != []
                   for node in self.reaction_graph.nodes())


class GraphBuilder():
    """
    A class for building the reaction graph.
    """

    def __init__(self) -> None:
        self._reaction_graph = DiGraph()

    def _add_node(self, node_id: str, type: NodeType, label: str) -> None:
        """
        Adding a node if the node is not already in the graph.

        Args:
            node_id: The ID of a graph node e.g. a domain will have the ID A_[dom]
            type: The type of a node e.g. component, domain, residue
            label: The label of the node e.g. the domain name (dom)

        Returns:

        """
        if not self._reaction_graph.has_node(node_id):
            logger.debug('Adding new node node_id: {0} label: {1}, type: {2}'.format(node_id, label, type))
            self._reaction_graph.add_node(node_id, label=label, type=type.value)

    def _add_edge(self, source: str, target: str, interaction: EdgeType, width: EdgeWith) -> None:
        """
        Adding an edge to the graph.

        Note:
             Internal edges are edges within the different levels of an specification.
             External edges are edges between specific resolution levels of two specifications
                e.g. between component and domain, domain and domain and so on.

        Args:
            source: The source of an edge e.g. component node, domain node, residue node.
            target: The target of an edge e.g. component node, domain node, residue node.
            interaction (EdgeType): The type of an interaction.
            width (EdgeWith): The width of an edge.

        Returns:
            None

        """
        if not self._reaction_graph.has_edge(source, target):
            logger.debug('Adding new edge source: {0} target: {1}, interaction: {2}, width: {3}'.format(source,
                                                                                                        target,
                                                                                                        interaction.value,
                                                                                                        width))
            self._reaction_graph.add_edge(source, target, interaction=interaction.value, width=width.value)
        elif width == EdgeWith.external:
            logger.debug(
                'Adding replace inner edge with external edge source: {0} target: {1}, interaction: {2}, width: {3}'.format(
                    source,
                    target,
                    interaction.value,
                    width))
            self._reaction_graph.add_edge(source, target, interaction=interaction.value, width=width.value)

    def add_external_edge(self, source: Spec, target: Spec, type: EdgeType) -> None:
        """
        Adding an external edge.

        Note:
            An external edges is an edge between two specific resolution levels of two specifications respectively
                e.g. between component and domain, domain and domain and so on.
        Args:
            source: The source specification.
            target: The target specification.
            type (EdgeType): The type of this edge e.g. interaction, modification etc.

        Returns:
            None

        """
        logger.info('Adding external edge source: {0} target: {1}, interaction: {2}'.format(
            get_node_id(source, source.resolution),
            get_node_id(target, target.resolution),
            type))
        self._add_edge(get_node_id(source, source.resolution), get_node_id(target, target.resolution),
                       interaction=type, width=EdgeWith.external)

    def add_spec_information(self, specification: Spec) -> None:
        """
        Adding specification information to the reaction graph.

        Args:
            specification: The specification of a reaction reactant

        Returns:
            None

        """

        def _add_spec_nodes() -> None:
            logger.info(
                'Adding component node -> id: {0}, label: {1}'.format(get_node_id(specification, NodeType.component),
                                                                      get_node_label(specification,
                                                                                     NodeType.component)))
            self._add_node(node_id=get_node_id(specification, NodeType.component), type=NodeType.component,
                           label=get_node_label(specification, NodeType.component))

            if specification.locus.domain:
                logger.info(
                    'Adding domain node -> id: {0}, label: {1}'.format(get_node_id(specification, NodeType.domain),
                                                                       get_node_label(specification, NodeType.domain)))
                self._add_node(get_node_id(specification, NodeType.domain), type=NodeType.domain,
                               label=get_node_label(specification, NodeType.domain))
            if specification.locus.residue:
                logger.info(
                    'Adding residue node id: {0}, label: {1}'.format(get_node_id(specification, NodeType.residue),
                                                                     get_node_label(specification, NodeType.residue)))
                self._add_node(get_node_id(specification, NodeType.residue), type=NodeType.residue,
                               label=get_node_label(specification, NodeType.residue))

        def _add_spec_edges() -> None:
            """
            Adding internal edges between nodes.

            Note:
                Internal edges are edges within the different levels of an specification.

            Returns:
                None

            """

            if specification.locus.domain:
                logger.info('Adding internal edge component -> domain source: {} target: {}'.format(
                    get_node_id(specification, NodeType.component),
                    get_node_id(specification, NodeType.domain)))
                self._add_edge(get_node_id(specification, NodeType.component),
                               get_node_id(specification, NodeType.domain), interaction=EdgeType.interaction,
                               width=EdgeWith.internal)

                if specification.locus.residue:
                    logger.info('Adding internal edge domain -> residue source: {} target: {}'.format(
                        get_node_id(specification, NodeType.domain),
                        get_node_id(specification, NodeType.residue)))
                    self._add_edge(get_node_id(specification, NodeType.domain),
                                   get_node_id(specification, NodeType.residue),
                                   interaction=EdgeType.interaction, width=EdgeWith.internal)
            elif specification.locus.residue:
                logger.info('Adding internal edge component -> residue source: {} target: {}'.format(
                    get_node_id(specification, NodeType.component),
                    get_node_id(specification, NodeType.residue)))
                self._add_edge(get_node_id(specification, NodeType.component),
                               get_node_id(specification, NodeType.residue),
                               interaction=EdgeType.interaction, width=EdgeWith.internal)

        logger.info('Adding nodes for {}'.format(specification))
        _add_spec_nodes()
        logger.info('Adding internal edges for {}'.format(specification))
        _add_spec_edges()

    def get_graph(self) -> DiGraph:
        """
        Returning the reaction graph

        Returns:
            The reaction graph (DiGraph).

        """
        return self._reaction_graph


def rxngraph_from_rxncon_system(rxncon_system: RxnConSystem) -> ReactionGraph:
    """
    Creating the reaction graph from the rxncon system.

    Args:
        rxncon_system: The reconstructed rxncon system.

    Returns:
        The reaction graph (ReactionGraph)

    """

    def get_reaction_type(rxn: Reaction) -> EdgeType:
        """
        Getting the type of a reaction.

        Args:
            rxn: Elemental rxncon reaction.

        Returns:
            The type of the reaction (EdgeType)

        """
        if len(rxn.components_lhs) > len(rxn.components_rhs):
            return EdgeType.degradation
        elif len(rxn.components_lhs) < len(rxn.components_rhs):
            return EdgeType.synthesis
        elif len(rxn.produced_states) == 1:
            return STATES[type(rxn.produced_states[0])]  # type: ignore
        elif len(rxn.produced_states) == 2:
            if all(STATES[type(produced_state)] == EdgeType.modification  # type: ignore
                   for produced_state in rxn.produced_states):
                return EdgeType.bimodification
            elif all(STATES[type(produced_state)] == EdgeType.interaction  # type: ignore
                     for produced_state in rxn.produced_states):
                return EdgeType.interaction
            else:
                return EdgeType.unknown
        else:
            return EdgeType.unknown

    def add_reaction_to_graph(rxn: Reaction) -> None:
        """
        Adding a reaction to the reaction graph.

        Args:
            rxn: Elemental rxncon reaction.

        Returns:
            None

        """
        edge_type = get_reaction_type(rxn)
        if edge_type is EdgeType.synthesis:
            for syn_comp in rxn.synthesised_components:
                logger.info('Adding synthesis of {0} -> {1}'.format(rxn['$x'], syn_comp))
                builder.add_spec_information(rxn['$x'])
                builder.add_spec_information(syn_comp)
                builder.add_external_edge(rxn['$x'], syn_comp, edge_type)
        else:
            logger.info('Adding {2} of {0} -> {1}'.format(rxn['$x'], rxn['$y'], get_reaction_type(rxn)))
            builder.add_spec_information(rxn['$x'])
            builder.add_spec_information(rxn['$y'])
            builder.add_external_edge(rxn['$x'], rxn['$y'], get_reaction_type(rxn))

    builder = GraphBuilder()

    for reaction in rxncon_system.reactions:
        if not isinstance(reaction, OutputReaction):
            logger.debug('rxngraph_from_rxncon_system: {}'.format(str(reaction)))
            add_reaction_to_graph(reaction)

    return ReactionGraph(builder.get_graph())


def get_node_id(specification: Spec, node_type: Union[NodeType, LocusResolution]) -> str:
    """
    Building the respective node_id

    Args:
        specification: Specification of an elemental state
        node_type: type of the node to get the proper ID

    Returns:
        The string ID of the node corresponding to asked conditions.

    Raises:
        AssertionError if non of the conditions are fulfilled.

    """

    if node_type in [NodeType.component, LocusResolution.component]:
        return str(specification.to_component_spec())

    if node_type in [NodeType.domain, LocusResolution.domain] and specification.locus.domain:
        return '{0}_[{1}]'.format(specification.name, specification.locus.domain)

    if node_type in [NodeType.residue, LocusResolution.residue] and specification.locus.residue:
        return str(specification)

    raise AssertionError


def get_node_label(specification: Spec, node_type: NodeType) -> str:
    """
    Asks for a label a specific specification resolution.

    Args:
        specification:
        node_type:

    Returns:

    """
    if node_type == NodeType.component:
        return re.match('[a-zA-Z0-9]+', str(specification)).group()

    if node_type == NodeType.domain and specification.locus.domain:
        return specification.locus.domain

    if node_type == NodeType.residue and specification.locus.residue:
        return specification.locus.residue

    raise AssertionError
