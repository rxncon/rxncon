from enum import Enum
from networkx import DiGraph
from rxncon.core.reaction import Reaction
from rxncon.core.rxncon_system import RxnConSystem
from rxncon.core.spec import Spec


class EdgeWith(Enum):
    internal = '0'
    external = '10'


class EdgeType(Enum):
    interaction  = 'interaction'
    modification = 'modification'
    synthesis    = 'synthesis'
    degradation  = 'degradation'


class NodeType(Enum):
    component = '80'
    domain    = '40'
    residue   = '20'


def get_node_id(specification: Spec, node_type: NodeType):
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

    if node_type == NodeType.component:
        return specification.component_name

    if node_type == NodeType.domain and specification.locus.domain:
        return '{0}_[{1}]'.format(specification.component_name, specification.locus.domain)

    if node_type == NodeType.residue and specification.locus.residue:
            return str(specification)

    raise AssertionError

def get_node_label(specification: Spec, node_type: NodeType):
    """
    Asks for a label a specific specification resolution.

    Args:
        specification:
        node_type:

    Returns:

    """
    if node_type == NodeType.component:
        return specification.component_name

    if node_type == NodeType.domain and specification.locus.domain:
        return specification.locus.domain

    if node_type == NodeType.residue and specification.locus.residue:
            return specification.locus.residue

    raise AssertionError

INTERACTION_STATES = {
    # Interaction state.
    '$x--$y': [
        lambda state, builder: builder.add_spec_nodes(state.specs[0]),
        lambda state, builder: builder.add_spec_nodes(state.specs[1]),
        lambda state, builder: builder.add_spec_edges(state.specs[0]),
        lambda state, builder: builder.add_spec_edges(state.specs[1]),
        lambda state, builder: builder.add_external_edge(source=get_node_id(state.specs[0], NodeType.domain),
                                                         target=get_node_id(state.specs[1], NodeType.domain),
                                                         type=EdgeType.interaction)
    ],
    # Self-interaction state.
    '$x--[$y]': [
        lambda state, builder: builder.add_node(node_id=get_node_id(state.specs[0], NodeType.component), type=NodeType.component, label=state.specs[0].component_name),
        lambda state, builder: builder.add_node(node_id=state.specs[1], type=NodeType.domain, label=state.specs[1]),
        lambda state, builder: builder.add_spec_edges(state.specs[0]),
        lambda state, builder: builder.add_spec_edges(state.specs[1]),
        lambda state, builder: builder.add_external_edge(source=get_node_id(state.specs[0], NodeType.domain), target=state.specs[1], type=EdgeType.interaction)
    ],
}

class ReactionGraph:
    def __init__(self, rxncon_system: RxnConSystem):
        self.rxncon_system = rxncon_system
        self.builder = GraphBuilder()
        #self.reaction_graph = DiGraph()

    def to_graph(self):
        for reaction in self.rxncon_system.reactions:
            self.add_reaction_information_to_graph(reaction)
        return self.builder.get_graph()

    def add_reaction_information_to_graph(self, reaction: Reaction):

        def _add_trans_reaction_to_graph(reaction: Reaction):
            for product_states in reaction.produced_states:
                for func in INTERACTION_STATES[product_states.repr_def]:
                    func(product_states, self.builder)

        def _add_cis_reaction_to_graph(reaction: Reaction):
            pass


        if self.is_trans(reaction):
            _add_trans_reaction_to_graph(reaction)
        elif self.is_cis(reaction):
            _add_cis_reaction_to_graph(reaction)
        else:
            raise AssertionError

    def is_cis(self, reaction: Reaction) -> bool:
        return len(reaction.components_lhs) == 1 and len(reaction.components_rhs) == 1

    def is_trans(self, reaction: Reaction) -> bool:
        return not self.is_cis(reaction)


class GraphBuilder():
    def __init__(self):
        self._reaction_graph = DiGraph()

    def _add_node(self, node_id: str, type: NodeType, label: str):
        self._reaction_graph.add_node(node_id, label=label, type=type.value)

    def _add_edge(self, source: str, target: str, interaction: EdgeType, width: EdgeWith):
        self._reaction_graph.add_edge(source, target, interaction=interaction.value, width=width.value)

    def add_external_edge(self, source: str, target: str, type: EdgeType):
        self._add_edge(source, target, interaction=type, width=EdgeWith.external)

    def add_spec_nodes(self, specification: Spec):
        self._add_node(node_id=get_node_id(specification, NodeType.component), type=NodeType.component,
                       label=get_node_label(specification, NodeType.component))

        if specification.locus.domain:
            self._add_node(get_node_id(specification, NodeType.domain), type=NodeType.domain,
                           label=get_node_label(specification, NodeType.domain))
        if specification.locus.residue:
            self._add_node(get_node_id(specification, NodeType.residue), type=NodeType.residue,
                           label=get_node_label(specification, NodeType.residue))

    def add_spec_edges(self, specification: Spec):

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

    def get_graph(self):
        return self._reaction_graph
