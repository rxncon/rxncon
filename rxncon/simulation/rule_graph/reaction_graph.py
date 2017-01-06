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
        lambda state, builder: builder.add_site(state['$x']) if builder.name == state['$x'].component_name \
                               else builder.add_site(state['$y']),
    ],
    # Self-interaction state.
    '$x--[$y]': [
        lambda state, builder: builder.add_node(node_id=get_node_id(state.specs[0], NodeType.component), type=NodeType.component, label=state.specs[0].component_name),
        lambda state, builder: builder.add_node(node_id=state.specs[1], type=NodeType.domain, label=state.specs[1]),
        lambda state, builder: builder.add_edge()
    ],
}

class ReactionGraph:
    def __init__(self, rxncon_system: RxnConSystem):
        self.rxncon_system = rxncon_system
        #self.reaction_graph = DiGraph()

    def to_graph(self):
        for reaction in self.rxncon_system.reactions:
            self.add_reaction_information_to_graph(reaction)

    def add_reaction_information_to_graph(self, reaction: Reaction):

        def _add_trans_reaction_to_graph(reaction: Reaction):
            pass

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

    def add_node(self, node_id: str, type: NodeType, label: str):
        self.reaction_graph.add_node(node_id, label=label, type=type)

    def add_edge(self, source: str, target: str, width: EdgeWith, type: EdgeType):
        self.reaction_graph.add_edge(source, target, interaction=type, width=width)

    def get_graph(self):
        return self._reaction_graph
