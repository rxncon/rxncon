from typing import Union, List, Dict
import re
import os
import copy
from enum import Enum, unique
from networkx import DiGraph
from rxncon.core.rxncon_system import RxnConSystem, Reaction, ReactionTerm, Contingency
from rxncon.core.effector import StateEffector, AndEffector, NotEffector, OrEffector, Effector
from xml.dom import minidom

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


def layout2regulatory_graph(no_layout_str: str, template_file_str: str) -> str:
    def check_filepath(file_path):
        if not os.path.exists(file_path):
            assert "File path {} does not exist!".format(file_path)

    def apply_template_layout(no_layout_str: str, template_file_str: str) -> str:
        #_get_node_information(no_layout_str)

        pass

    def _get_labels_and_coordinates_dict(xmldoc) -> Dict:

        graphics_list = xmldoc.getElementsByTagName('graphics')
        coordinates_dict = {graphic.parentNode.getAttribute('label'): {"x": graphic.getAttribute('x'),
                                                                       "y": graphic.getAttribute('y'),
                                                                       "z": graphic.getAttribute('z')}
                            for graphic in graphics_list if graphic.attributes.values() and graphic.parentNode.tagName == "node"}

        return coordinates_dict

    check_filepath(template_file_str)
    xmldoc_no_layout = minidom.parseString(no_layout_str)
    node_list_no_layout = xmldoc_no_layout.getElementsByTagName('node')

    xmldoc_template = minidom.parse(template_file_str)
    node_list_template = xmldoc_template.getElementsByTagName('node')
    template_coordinates = _get_labels_and_coordinates_dict(xmldoc_template)

    for no_layout_node in node_list_no_layout:
        if no_layout_node.getAttribute('label') in template_coordinates:
            node_name = no_layout_node.getAttribute('label')
            element = xmldoc_no_layout.createElement("graphics")
            element.setAttribute("x", template_coordinates[node_name]["x"])
            element.setAttribute("y", template_coordinates[node_name]["y"])
            element.setAttribute("z", template_coordinates[node_name]["z"])
            element.appendChild(xmldoc_no_layout.createTextNode(''))
            no_layout_node.appendChild(element)

    return xmldoc_no_layout.toprettyxml()



# def insert_graphics_elements_into_graph(graph_data_unlayouted, graph_data_layouted, node_names_layouted):
#     coordinates_dict = _get_labels_and_coordinates_dict(graph_data_layouted["xmldoc"])
#     for item in graph_data_unlayouted["node_names_dicts"]:
#         if item["name"] in node_names_layouted and item["name"] in coordinates_dict:
#             element = graph_data_unlayouted["xmldoc"].createElement("graphics")
#             element.setAttribute("x", coordinates_dict[item["name"]]["x"])
#             element.setAttribute("y", coordinates_dict[item["name"]]["y"])
#             element.setAttribute("z", coordinates_dict[item["name"]]["z"])
#             element.appendChild(graph_data_unlayouted["xmldoc"].createTextNode(''))
#             target_node = graph_data_unlayouted["node_list"].item(item["index"])
#             target_node.appendChild(element)
#     return graph_data_unlayouted
#

# def apply_template_layout(request, graph_file_path):
#     template_file = request.FILES.get('template')
#     graph_data_unlayouted = get_graph_nodes_and_lables_from_file(graph_file_path)
#     graph_data_layouted = get_graph_nodes_and_lables_from_file(template_file)
#     node_names_layouted = [d["name"] for d in graph_data_layouted["node_names_dicts"]]
#     graph_data_now_layouted = insert_graphics_elements_into_graph(graph_data_unlayouted, graph_data_layouted, node_names_layouted)
#     graph_file = open(graph_file_path, "w")
#     graph_data_now_layouted["xmldoc"].writexml(graph_file)
#     graph_string = graph_data_now_layouted["xmldoc"].toprettyxml()
#     graph_file.close()
#
#     return graph_file, graph_string
#