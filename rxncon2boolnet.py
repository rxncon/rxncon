#!/usr/bin/python3

from typing import Tuple, Optional
import os, sys
import click
import logging
import click_log
import colorama
from enum import Enum

from rxncon.core.rxncon_system import RxnConSystem
from rxncon.input.excel_book.excel_book import ExcelBook
from rxncon.simulation.boolean.boolean_model import boolean_model_from_rxncon, \
    SmoothingStrategy, KnockoutStrategy, OverexpressionStrategy
from rxncon.simulation.boolean.boolnet_from_boolean_model import boolnet_from_boolean_model


colorama.init()
LOGGER = logging.getLogger(__name__)


class QuantitativeContingencyStrategy(Enum):
    strict = 'strict'
    ignore = 'ignore'


def boolnet_strs_from_rxncon(rxncon: RxnConSystem, smoothing_strategy: SmoothingStrategy, knockout_strategy: KnockoutStrategy,
                             overexpression_strategy: OverexpressionStrategy, k_plus_strategy: QuantitativeContingencyStrategy,
                             k_minus_strategy: QuantitativeContingencyStrategy) \
        -> Tuple[str, str, str]:
    def sort_key(key_val_pair):
        k, v = key_val_pair
        return k[0], int(k[1:])

    if k_plus_strategy == QuantitativeContingencyStrategy.strict:
        k_plus_strict = True
    elif k_plus_strategy == QuantitativeContingencyStrategy.ignore:
        k_plus_strict = False
    else:
        raise AssertionError

    if k_minus_strategy == QuantitativeContingencyStrategy.strict:
        k_minus_strict = True
    elif k_minus_strategy == QuantitativeContingencyStrategy.ignore:
        k_minus_strict = False
    else:
        raise AssertionError

    model_str, symbol_dict, initial_val_dict = \
        boolnet_from_boolean_model(boolean_model_from_rxncon(rxncon, smoothing_strategy=smoothing_strategy,
                                                             knockout_strategy=knockout_strategy,
                                                             overexpression_strategy=overexpression_strategy,
                                                             k_plus_strict=k_plus_strict, k_minus_strict=k_minus_strict))

    symbol_str      = '\n'.join('{0}, {1}'.format(boolnet_sym, rxncon_sym) for boolnet_sym, rxncon_sym
                                in sorted(symbol_dict.items(), key=sort_key)) + '\n'

    initial_val_str = '\n'.join('{0}, {1: <5}  , #  {2}'.format(boolnet_sym, initial_val, symbol_dict[boolnet_sym])
                                for boolnet_sym, initial_val in sorted(initial_val_dict.items(), key=sort_key)) + '\n'

    return model_str, symbol_str, initial_val_str


def write_boolnet(excel_filename: str, smoothing_strategy: SmoothingStrategy, knockout_strategy: KnockoutStrategy,
                  overexpression_strategy: OverexpressionStrategy, k_plus_strategy: QuantitativeContingencyStrategy,
                  k_minus_strategy: QuantitativeContingencyStrategy, base_name: Optional[str]=None):
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
    model_str, symbol_str, initial_val_str = boolnet_strs_from_rxncon(rxncon_system, smoothing_strategy, knockout_strategy,
                                                                      overexpression_strategy, k_plus_strategy, k_minus_strategy)

    print('Writing BoolNet model file [{}] ...'.format(boolnet_model_filename))
    with open(boolnet_model_filename, mode='w') as f:
        f.write(model_str)

    print('Writing BoolNet symbol file [{}] ...'.format(boolnet_symbol_filename))
    with open(boolnet_symbol_filename, mode='w') as f:
        f.write(symbol_str)

    print('Writing BoolBet initial value file [{}] ...'.format(boolnet_initial_val_filename))
    with open(boolnet_initial_val_filename, mode='w') as f:
        f.write(initial_val_str)


valid_smoothing_strategies      = [strategy.value for strategy in SmoothingStrategy.__members__.values()]
valid_knockout_strategies       = [strategy.value for strategy in KnockoutStrategy.__members__.values()]
valid_overexpression_strategies = [strategy.value for strategy in OverexpressionStrategy.__members__.values()]
valid_quantitative_contingency_strategies = [strategy.value for strategy in QuantitativeContingencyStrategy.__members__.values()]

def validate_smoothing_strategy(ctx, param, value):
    try:
        SmoothingStrategy(value)
        return value
    except ValueError:
        raise click.BadParameter('Valid strategies are: {}'.format(', '.join(valid_smoothing_strategies)))


def validate_knockout_strategy(ctx, param, value):
    try:
        KnockoutStrategy(value)
        return value
    except ValueError:
        raise click.BadParameter('Valid strategies are: {}'.format(', '.join(valid_knockout_strategies)))


def validate_overexpression_strategy(ctx, param, value):
    try:
        OverexpressionStrategy(value)
        return value
    except ValueError:
        raise click.BadParameter('Valid strategies are: {}'.format(', '.join(valid_overexpression_strategies)))


def validate_quantitative_contingency_strategy(ctx, param, value):
    try:
        QuantitativeContingencyStrategy(value)
        return value
    except ValueError:
        raise click.BadParameter('Valid strategies are: {}'.format(', '.join(valid_quantitative_contingency_strategies)))


@click.command()
@click.option('--smoothing', default='no_smoothing',
              help='Smoothing strategy. Default: no_smoothing. Choices: {}'.format(', '.join(valid_smoothing_strategies)),
              callback=validate_smoothing_strategy)
@click.option('--knockout', default='no_knockout',
              help='Generate knockouts. Default: no_knockout. Choices: {}'.format(', '.join(valid_knockout_strategies)),
              callback=validate_knockout_strategy)
@click.option('--overexpression', default='no_overexpression',
              help='Generate overexpressions. Default: no_overexpression. Choices: {}'.format(', '.join(valid_overexpression_strategies)),
              callback=validate_overexpression_strategy)
@click.option('--k_plus', default='strict',
              help='Strategy for handling k+ contingencies. Default: strict. Choices: {}'.format(', '.join(valid_quantitative_contingency_strategies)),
              callback=validate_quantitative_contingency_strategy)
@click.option('--k_minus', default='strict',
              help='Strategy for handling k- contingencies. Default: strict. Choices: {}'.format(', '.join(valid_quantitative_contingency_strategies)),
              callback=validate_quantitative_contingency_strategy)
@click.option('--output', default=None,
              help='Base name for output files. Default: \'fn\' for input file \'fn.xls\'')
@click.argument('excel_file')
@click_log.simple_verbosity_option(default='WARNING')
@click_log.init()
def run(overexpression, knockout, smoothing, output, excel_file, k_plus, k_minus):
    smoothing_strategy = SmoothingStrategy(smoothing)
    knockout_strategy = KnockoutStrategy(knockout)
    overexpression_strategy = OverexpressionStrategy(overexpression)
    k_plus_strategy = QuantitativeContingencyStrategy(k_plus)
    k_minus_strategy = QuantitativeContingencyStrategy(k_minus)
    write_boolnet(excel_file, smoothing_strategy, knockout_strategy, overexpression_strategy,
                  k_plus_strategy, k_minus_strategy, output)


def setup_logging_colors():
    click_log.ColorFormatter.colors = {
        'error': dict(fg='red'),
        'exception': dict(fg='red'),
        'critical': dict(fg='red'),
        'debug': dict(fg='yellow'),
        'warning': dict(fg='yellow'),
        'info': dict(fg='yellow')
    }

    def format(self, record):
        if not record.exc_info:
            level = record.levelname.lower()
            if level in self.colors:
                padding_size = 7  # Assume just INFO / DEBUG entries.

                prefix = click.style('{}: '.format(level).ljust(padding_size),
                                     **self.colors[level])

                prefix += click.style('{} '.format(record.name), fg='blue')

                msg = record.msg
                if isinstance(msg, bytes):
                    msg = msg.decode(sys.getfilesystemencoding(),
                                     'replace')
                elif not isinstance(msg, str):
                    msg = str(msg)
                record.msg = '\n'.join(prefix + x for x in msg.splitlines())

        return logging.Formatter.format(self, record)

    click_log.ColorFormatter.format = format

if __name__ == '__main__':
    try:
        setup_logging_colors()
        run()
    except Exception as e:
        print('ERROR: {}\n{}\nPlease re-run this command with the \'-v DEBUG\' option.'.format(type(e), e))


