import os
from typing import Dict, Union
from networkx import DiGraph
from xml.dom import minidom


class XGMML:
    """
    Definition of the XGMML format.

    Args:
        graph (DiGraph): A graph which will be written in XGMML format.
        graph_name: The name of the graph.

    """

    def __init__(self, graph: DiGraph, graph_name: str) -> None:
        self.graph = graph
        self.graph_name = graph_name

    def to_string(self) -> str:
        """
        Translating the graph information into a xgmml string representation.

        Returns:
            A xgmml string representation.

        """
        xgmml = [self._header_string(), self._nodes_string(), self._edges_string(), self._footer_string()]
        return "\n".join(xgmml)

    def to_file(self, file_path: str, force: bool=False) -> None:

        """
        Writes the xgmml graph into a file.

        Args:
            file_path: path to the file the xgmml graph is written to.
            force:     overwrite file if it already exists.
        Returns:
            None

        Raises:
            FileExistsError: If the file already exists.
            NotADirectoryError: If the path to the file does not exists.

        """
        path, file = os.path.split(file_path)
        if path and os.path.exists(path):
            if force:
                self._write_to_file(file_path)
            elif not os.path.isfile(file_path):
                self._write_to_file(file_path)
            else:
                raise FileExistsError("{0} exists! remove file and run again".format(file_path))
        elif not path:
            if not os.path.isfile(file):
                self._write_to_file(file_path)
            else:
                if force:
                    self._write_to_file(file_path)
                print(os.path.dirname(file_path))
                raise FileExistsError("{0} exists! remove file and run again".format(file_path))
        elif path and not os.path.exists(path):
            raise NotADirectoryError("Path {0} does not exists.".format(path))

    def _write_to_file(self, file_path: str) -> None:
        with open(file_path, mode='w') as writehandle:
            writehandle.write(self.to_string())

    def _header_string(self) -> str:
        """
        Defining the header of the xgmml file.

        Returns:
            The header of the xgmml file.

        """
        return """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
    <graph directed="1"  xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns="http://www.cs.rpi.edu/XGMML">
    <att name="selected" value="1" type="boolean" />
    <att name="name" value="{0}" type="string"/>
    <att name="shared name" value="{0}" type="string"/>
    """.format(self.graph_name)

    def _footer_string(self) -> str:
        """
        Defining the footer of the xgmml file.

        Returns:
            Footer of the xgmml file.

        """
        return '</graph>'

    def _nodes_string(self) -> str:
        """
        Translating the graph node information into xgmml node information.

        Returns:
            String of all node information in xgmml format.

        """
        nodes = []
        for graph_node in self.graph.nodes(data=True):
            id = graph_node[0]
            attr = dict(graph_node[1])

            if 'label' in attr:
                label = attr['label']
                del attr['label']
            else:
                label = id

            node = '<node id="{id}" label="{label}">'.format(id=id, label=label)
            node += '<att name="{}" value="{}" type="string"/>'.format("rxnconID", id)

            for name, value in attr.items():
                node += self._format_attribute(name, value)

            node += '</node>'
            nodes.append(node)
        return "\n".join(nodes)

    def _edges_string(self) -> str:
        """
        Translating the graph edge information into xgmml format.

        Returns:
            String of all edge information in xgmml format.

        """
        edges = []

        for graph_edge in self.graph.edges(data=True):
            edge = '<edge source="{}" target="{}">'.format(graph_edge[0], graph_edge[1])
            for name, value in graph_edge[2].items():
                edge += self._format_attribute(name, value)

            edge += '</edge>'
            edges.append(edge)
        return "\n".join(edges)

    def _format_attribute(self, name: str, value: Union[float, int, str]) -> str:
        if isinstance(value, float):
            return '<att name="{}" value="{}" type="double"/>'.format(name, value)
        elif isinstance(value, int):
            return '<att name="{}" value="{}" type="integer"/>'.format(name, value)
        else:
            return '<att name="{}" value="{}" type="string"/>'.format(name, value)


def map_layout2xgmml(no_layout_graph_str: str, layout_template: str, str_template: bool=False) -> str:
    """
    Mapping the layout information of a xgmml file to a xgmml string.

    Args:
        no_layout_graph_str: xgmml graph with no layout information.
        layout_template: xgmml graph with layout information
        str_template: @todo

    Returns:
        XGMML String with mapped layout information.

    """

    def _check_filepath(file_path: str) -> None:
        """
        checks if file path/file already exists

        Args:
            file_path: Path to a file.

        Returns:
            None

        Raises:
            FileExistsError: If file exists.
            NotADirectoryError: If path does not exists.

        """

        path, file = os.path.split(file_path)
        if path and os.path.exists(path) and not os.path.isfile(file_path):
            raise FileNotFoundError("{0} does not exists!".format(file_path))
        elif not path and not os.path.isfile(file):
            raise FileNotFoundError("{0} does not exists!".format(file_path))
        elif path and not os.path.exists(path):
            raise NotADirectoryError("Path {0} does not exists.".format(path))

    def _apply_template_layout() -> None:
        """
        Writes the template layout information to the xml nodes with no layout information.

        Returns:
            None

        """
        nonlocal template_coordinates
        nonlocal node_list_no_layout
        nonlocal xmldoc_no_layout

        for no_layout_node in node_list_no_layout:
            rxnconID = _get_rxnconID(no_layout_node)
            if rxnconID in template_coordinates:
                element = xmldoc_no_layout.createElement("graphics")
                element.setAttribute("x", template_coordinates[rxnconID]["x"])
                element.setAttribute("y", template_coordinates[rxnconID]["y"])
                element.setAttribute("z", template_coordinates[rxnconID]["z"])
                element.appendChild(xmldoc_no_layout.createTextNode(''))
                no_layout_node.appendChild(element)
    if not str_template:
        _check_filepath(layout_template)
    xmldoc_no_layout = minidom.parseString(no_layout_graph_str)
    node_list_no_layout = xmldoc_no_layout.getElementsByTagName('node')
    template_coordinates = _get_labels_and_coordinates_dict(minidom.parse(layout_template))

    _apply_template_layout()

    return xmldoc_no_layout.toprettyxml()


def _get_rxnconID(element: minidom.Element) -> str:
    for child in element.childNodes:
        if child.attributes and child.getAttribute('name') == "rxnconID":
            return child.getAttribute('value')
    raise AssertionError('Could not find rxnconID for element {}'.format(element))


def _get_labels_and_coordinates_dict(xmldoc: minidom.Document) -> Dict[str, Dict[str, str]]:
    """
    Creates a mapping of node names and their coordinates.

    Args:
        xmldoc: xml information.

    Returns:
        Dict[str, Dict[str, str]]
        Dictionary of nodes and coordinates.

    """

    template_coordinates = {}  # type: Dict[str, Dict[str, str]]
    for graphic in xmldoc.getElementsByTagName('graphics'):
        if graphic.attributes.values() and graphic.parentNode.tagName == "node":
            rxnconID = _get_rxnconID(graphic.parentNode)
            if rxnconID not in template_coordinates:
                template_coordinates[rxnconID] = {"x": graphic.getAttribute('x'),
                                                  "y": graphic.getAttribute('y'),
                                                  "z": graphic.getAttribute('z')}
            else:
                raise AssertionError
    return template_coordinates
