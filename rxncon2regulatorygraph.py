from typing import Tuple
import os
import click

from rxncon.core.rxncon_system import RxnConSystem
from rxncon.input.excel_book.excel_book import ExcelBook
from rxncon.simulation.rule_graph.regulatory_graph import RegulatoryGraph
from rxncon.simulation.rule_graph.graphML import XGMML


def write_xgmml(excel_filename: str, base_name=None):
    if not base_name:
        base_name = os.path.splitext(os.path.basename(excel_filename))[0]

    base_path = os.path.dirname(excel_filename)

    graph_filename = os.path.join(base_path, '{0}.xgmml'.format(base_name))

    print('Reading in Excel file [{}] ...'.format(excel_filename))
    excel_book = ExcelBook(excel_filename)

    rxncon_system = excel_book.rxncon_system
    print('Constructed rxncon system: [{} reactions], [{} contingencies]'
          .format(len(rxncon_system.reactions), len(rxncon_system.contingencies)))

    print('Generating regulatory graph output...')
    reg_system = RegulatoryGraph(rxncon_system)
    graph = reg_system.to_graph()


    print('Writing regulatory graph file [{}] ...'.format(graph_filename))
    gml_system = XGMML(graph, "{}".format(base_name))
    gml_system.to_file(graph_filename)


@click.command()
@click.option('--output', default=None,
              help='Base name for output files. Default: \'fn\' for input file \'fn.xls\'')
@click.argument('excel_file')
def run(output, excel_file):
    write_xgmml(excel_file, output)

if __name__ == '__main__':
    run()
