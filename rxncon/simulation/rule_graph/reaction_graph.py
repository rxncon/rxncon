from enum import Enum
from networkx import DiGraph
from rxncon.core.reaction import Reaction
from rxncon.core.rxncon_system import RxnConSystem


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

INTERACTION_STATES = {
    # Interaction state.
    '$x--$y': [
        lambda state, builder: builder.add_site(state['$x']) if builder.name == state['$x'].component_name \
                               else builder.add_site(state['$y']),
    ],
    # Self-interaction state.
    '$x--[$y]': [
        lambda state, builder: builder.add_node(state.specs[0]),
        lambda state, builder: builder.add_node(state.specs[1])
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
