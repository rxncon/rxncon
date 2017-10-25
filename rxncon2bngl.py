#!/usr/bin/python3

import os, logging
import click
import click_log
import colorama
import sys

from rxncon.input.excel_book.excel_book import ExcelBook
from rxncon.simulation.rule_based.rule_based_model import rule_based_model_from_rxncon
from rxncon.simulation.rule_based.bngl_from_rule_based_model import bngl_from_rule_based_model


colorama.init()
LOGGER = logging.getLogger(__name__)


def write_bngl(excel_filename: str, base_name=None):
    if not base_name:
        base_name = os.path.splitext(os.path.basename(excel_filename))[0]

    base_path = os.path.dirname(excel_filename)

    bngl_model_filename = os.path.join(base_path, '{0}.bngl'.format(base_name))

    print('Reading in Excel file [{}] ...'.format(excel_filename))
    excel_book = ExcelBook(excel_filename)

    rxncon_system = excel_book.rxncon_system
    print('Constructed rxncon system: [{} reactions], [{} contingencies], [{} components], [{} elemental states]'
          .format(len(rxncon_system.reactions), len(rxncon_system.contingencies), len(rxncon_system.components()),
                  len(rxncon_system.states)))

    print('Generating BNGL output ...')
    rbm = rule_based_model_from_rxncon(rxncon_system)
    print('Constructed rule-based model: [{} molecule types], [{} rules], [{} observables]'
          .format(len(rbm.mol_defs), len(rbm.rules), len(rbm.observables)))
    model_str = bngl_from_rule_based_model(rbm)

    print('Writing BNGL model file [{}] ...'.format(bngl_model_filename))
    with open(bngl_model_filename, mode='w') as f:
        f.write(model_str)


@click.command()
@click.option('--output', default=None,
              help='Base name for output files. Default: \'fn\' for input file \'fn.xls\'')
@click.argument('excel_file')
@click_log.simple_verbosity_option(default='WARNING')
@click_log.init()
def run(output, excel_file):
    write_bngl(excel_file, output)

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
