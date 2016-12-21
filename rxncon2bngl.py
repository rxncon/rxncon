import os
import click

from rxncon.core.rxncon_system import RxnConSystem
from rxncon.input.excel_book.excel_book import ExcelBook
from rxncon.simulation.rule_based.rule_based_model import rule_based_model_from_rxncon
from rxncon.simulation.rule_based.bngl_from_rule_based_model import bngl_from_rule_based_model


def bngl_str_from_rxncon(rxncon: RxnConSystem) -> str:
    rbm = rule_based_model_from_rxncon(rxncon)
    return bngl_from_rule_based_model(rbm)


def write_bngl(excel_filename: str, base_name=None):
    if not base_name:
        base_name = os.path.splitext(os.path.basename(excel_filename))[0]

    base_path = os.path.dirname(excel_filename)

    bngl_model_filename = os.path.join(base_path, '{0}.bngl'.format(base_name))

    print('Reading in Excel file [{}] ...'.format(excel_filename))
    excel_book = ExcelBook(excel_filename)

    rxncon_system = excel_book.rxncon_system
    print('Constructed rxncon system: [{} reactions], [{} contingencies]'
          .format(len(rxncon_system.reactions), len(rxncon_system.contingencies)))

    print('Generating BNGL output ...')
    model_str = bngl_str_from_rxncon(rxncon_system)

    print('Writing BNGL model file [{}] ...'.format(bngl_model_filename))
    with open(bngl_model_filename, mode='w') as f:
        f.write(model_str)


@click.command()
@click.option('--output', default=None,
              help='Base name for output files. Default: \'fn\' for input file \'fn.xls\'')
@click.argument('excel_file')
def run(output, excel_file):
    write_bngl(excel_file, output)

if __name__ == '__main__':
    run()
