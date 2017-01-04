#!/usr/bin/python3

from typing import Tuple
import os
import click

from rxncon.core.rxncon_system import RxnConSystem
from rxncon.input.excel_book.excel_book import ExcelBook
from rxncon.simulation.boolean.boolean_model import boolnet_from_boolean_model, boolean_model_from_rxncon, \
    SmoothingStrategy


def boolnet_strs_from_rxncon(rxncon: RxnConSystem, smoothing_strategy: SmoothingStrategy) -> Tuple[str, str, str]:
    def sort_key(key_val_pair):
        k, v = key_val_pair
        return k[0], int(k[1:])

    model_str, symbol_dict, initial_val_dict = \
        boolnet_from_boolean_model(boolean_model_from_rxncon(rxncon, smoothing_strategy))

    symbol_str      = '\n'.join('{0}, {1}'.format(boolnet_sym, rxncon_sym) for boolnet_sym, rxncon_sym
                                in sorted(symbol_dict.items(), key=sort_key)) + '\n'

    initial_val_str = '\n'.join('{0}, {1: <5}  , #  {2}'.format(boolnet_sym, initial_val, symbol_dict[boolnet_sym])
                                for boolnet_sym, initial_val in sorted(initial_val_dict.items(), key=sort_key)) + '\n'

    return model_str, symbol_str, initial_val_str


def write_boolnet(excel_filename: str, smoothing_strategy: SmoothingStrategy, base_name=None):
    if not base_name:
        base_name = os.path.splitext(os.path.basename(excel_filename))[0]

    base_path = os.path.dirname(excel_filename)

    boolnet_model_filename = os.path.join(base_path, '{0}.boolnet'.format(base_name))
    boolnet_symbol_filename = os.path.join(base_path, '{0}_symbols.csv'.format(base_name))
    boolnet_initial_val_filename = os.path.join(base_path, '{0}_initial_vals.csv'.format(base_name))

    print('Reading in Excel file [{}] ...'.format(excel_filename))
    excel_book = ExcelBook(excel_filename)

    rxncon_system = excel_book.rxncon_system
    print('Constructed rxncon system: [{} reactions], [{} contingencies]'
          .format(len(rxncon_system.reactions), len(rxncon_system.contingencies)))

    print('Generating BoolNet output using smoothing strategy [{}] ...'.format(smoothing_strategy.name))
    model_str, symbol_str, initial_val_str = boolnet_strs_from_rxncon(rxncon_system, smoothing_strategy)

    print('Writing BoolNet model file [{}] ...'.format(boolnet_model_filename))
    with open(boolnet_model_filename, mode='w') as f:
        f.write(model_str)

    print('Writing BoolNet symbol file [{}] ...'.format(boolnet_symbol_filename))
    with open(boolnet_symbol_filename, mode='w') as f:
        f.write(symbol_str)

    print('Writing BoolBet initial value file [{}] ...'.format(boolnet_initial_val_filename))
    with open(boolnet_initial_val_filename, mode='w') as f:
        f.write(initial_val_str)


valid_strategies = [strategy.value for strategy in SmoothingStrategy.__members__.values()]


def validate_smoothing_strategy(ctx, param, value):
    try:
        strategy = SmoothingStrategy(value)
        return value
    except ValueError:
        raise click.BadParameter('Valid strategies are: {}'.format(', '.join(valid_strategies)))


@click.command()
@click.option('--smoothing', default='no_smoothing',
              help='Smoothing strategy. Default: no_smoothing. Choices: {}'.format(', '.join(valid_strategies)),
              callback=validate_smoothing_strategy)
@click.option('--output', default=None,
              help='Base name for output files. Default: \'fn\' for input file \'fn.xls\'')
@click.argument('excel_file')
def run(smoothing, output, excel_file):
    smoothing_strategy = SmoothingStrategy(smoothing)
    write_boolnet(excel_file, smoothing_strategy, output)

if __name__ == '__main__':
    run()
