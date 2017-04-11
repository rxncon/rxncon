import os
from xml.dom import minidom

from rxncon.input.excel_book.excel_book import ExcelBook
from rxncon.visualization.graphML import XGMML, map_layout2xgmml, _get_labels_and_coordinates_dict
from rxncon.visualization.regulatory_graph import SpeciesReactionGraph


PHEROMONE_XLS   = os.path.join(os.path.dirname(__file__), 'pheromone.xls')
PHEROMONE_XGMML = os.path.join(os.path.dirname(__file__), 'pheromone_layout.xgmml')


def test_regulatory_graph_layout_remains() -> None:
    """
    Testing if the new graph gets the correct layout of the old graph.

    Returns:
        None

    Raises:
        AssertionError: If the node coordinates of the new graph differ from expected coordinates.

    """

    book = ExcelBook(PHEROMONE_XLS)
    reg_graph = SpeciesReactionGraph(book.rxncon_system)
    gml_system = XGMML(reg_graph.to_graph(), "pheromone_layout")
    mapped_layout = map_layout2xgmml(gml_system.to_string(), PHEROMONE_XGMML)
    xmldoc_no_layout = minidom.parseString(mapped_layout)
    xmldoc_no_layout_info = _get_labels_and_coordinates_dict(xmldoc_no_layout)
    xmldoc_expected_layout = minidom.parse(PHEROMONE_XGMML)
    xmldoc_expected_layout_info = _get_labels_and_coordinates_dict(xmldoc_expected_layout)
    assert all(xmldoc_no_layout_info[no_layout_node] == xmldoc_expected_layout_info[no_layout_node]
               for no_layout_node in xmldoc_no_layout_info)
    assert all(xmldoc_no_layout_info[expected_layout_node] == xmldoc_expected_layout_info[expected_layout_node]
               for expected_layout_node in xmldoc_expected_layout_info)
