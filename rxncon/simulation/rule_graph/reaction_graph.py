from typing import Union
from enum import Enum
from networkx import DiGraph
import re
from rxncon.core.reaction import Reaction
from rxncon.core.rxncon_system import RxnConSystem
from rxncon.core.spec import Spec, Locus, LocusResolution


class EdgeWith(Enum):
    internal = '0'
    external = '10'


class EdgeType(Enum):
    interaction    = 'interaction'
    modification   = 'modification'
    bimodification = 'bimodification'
    synthesis      = 'synthesis'
    degradation    = 'degradation'
    unknown        = 'unknown'


class NodeType(Enum):
    component = '80'
    domain    = '40'
    residue   = '20'


STATES = {'$x--$y': EdgeType.interaction,
          '$x--[$y]': EdgeType.interaction,
          '$x-{$y}': EdgeType.modification}



class ReactionGraph:
    def __init__(self, reaction_graph: DiGraph) -> None:
        self.reaction_graph = reaction_graph
        self._validate_graph()

    def _validate_graph(self) -> None:
        pass


class GraphBuilder():
    def __init__(self) -> None:
        self._reaction_graph = DiGraph()

    def _add_node(self, node_id: str, type: NodeType, label: str) -> None:
        if not self._reaction_graph.has_node(node_id):
            self._reaction_graph.add_node(node_id, label=label, type=type.value)

    def _add_edge(self, source: str, target: str, interaction: EdgeType, width: EdgeWith) -> None:
        if not self._reaction_graph.has_edge(source, target):
            self._reaction_graph.add_edge(source, target, interaction=interaction.value, width=width.value)
        elif width == EdgeWith.external:
            self._reaction_graph.add_edge(source, target, interaction=interaction.value, width=width.value)

    def add_external_edge(self, source: Spec, target: Spec, type: EdgeType) -> None:

        self._add_edge(get_node_id(source, source.resolution), get_node_id(target, target.resolution),
                       interaction=type, width=EdgeWith.external)

    def add_spec_information(self, specification: Spec) -> None:

        def _add_spec_nodes() -> None:
            self._add_node(node_id=get_node_id(specification, NodeType.component), type=NodeType.component,
                           label=get_node_label(specification, NodeType.component))

            if specification.locus.domain:
                self._add_node(get_node_id(specification, NodeType.domain), type=NodeType.domain,
                               label=get_node_label(specification, NodeType.domain))
            if specification.locus.residue:
                self._add_node(get_node_id(specification, NodeType.residue), type=NodeType.residue,
                               label=get_node_label(specification, NodeType.residue))

        def _add_spec_edges() -> None:

            if specification.locus.domain:
                self._add_edge(get_node_id(specification, NodeType.component),
                               get_node_id(specification, NodeType.domain), interaction=EdgeType.interaction,
                               width=EdgeWith.internal)
                if specification.locus.residue:
                    self._add_edge(get_node_id(specification, NodeType.domain), get_node_id(specification, NodeType.residue),
                                   interaction=EdgeType.interaction, width=EdgeWith.internal)
            elif specification.locus.residue:
                self._add_edge(get_node_id(specification, NodeType.component), get_node_id(specification, NodeType.residue),
                           interaction=EdgeType.interaction, width=EdgeWith.internal)

        _add_spec_nodes()
        _add_spec_edges()

    def get_graph(self) -> DiGraph:
        return self._reaction_graph


def rxngraph_from_rxncon_system(rxncon_system: RxnConSystem) -> ReactionGraph:
    def get_reaction_type(rxn: Reaction) -> EdgeType:
        if len(rxn.components_lhs) > len(rxn.components_rhs):
            return EdgeType.degradation
        elif len(rxn.components_lhs) < len(rxn.components_rhs):
            return EdgeType.synthesis
        elif len(rxn.produced_states) == 1:
            return STATES[rxn.produced_states[0].repr_def]
        elif len(rxn.produced_states) == 2:
            if all(STATES[produced_state.repr_def] == EdgeType.modification for produced_state in rxn.produced_states):
                return EdgeType.bimodification
            else:
                return EdgeType.unknown
        else:
            return EdgeType.unknown

    def add_reaction_to_graph(rxn: Reaction) -> None:
        edge_type = get_reaction_type(rxn)
        if edge_type is EdgeType.synthesis:
            for syn_comp in rxn.synthesised_components:
                builder.add_spec_information(rxn['$x'])
                builder.add_spec_information(syn_comp)
                builder.add_external_edge(rxn['$x'], syn_comp, edge_type)
        else:
            builder.add_spec_information(rxn['$x'])
            builder.add_spec_information(rxn['$y'])
            builder.add_external_edge(rxn['$x'], rxn['$y'], get_reaction_type(rxn))

    builder = GraphBuilder()

    for reaction in rxncon_system.reactions:
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
        return re.match('[a-zA-Z0-9]+', str(specification)).group()

    if node_type in [NodeType.domain, LocusResolution.domain] and specification.locus.domain:
        return '{0}_[{1}]'.format(specification.component_name, specification.locus.domain)

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
