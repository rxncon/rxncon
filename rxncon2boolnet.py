from typing import Tuple, Optional
import os
import sys

from rxncon.core.rxncon_system import RxnConSystem
from rxncon.input.excel_book.excel_book import ExcelBook
from rxncon.simulation.boolean.boolean_model import boolnet_from_boolean_model, boolean_model_from_rxncon


def boolnet_strs_from_rxncon(rxncon: RxnConSystem) -> Tuple[str, str, str]:
    def sort_key(key_val_pair):
        k, v = key_val_pair
        return k[0], int(k[1:])

    model_str, symbol_dict, initial_val_dict = boolnet_from_boolean_model(boolean_model_from_rxncon(rxncon))

    symbol_str       = '\n'.join('{0}, {1}'.format(boolnet_sym, rxncon_sym) for boolnet_sym, rxncon_sym
                                 in sorted(symbol_dict.items(), key=sort_key)) + '\n'

    initial_val_str  = ', '.join(sym for sym, _ in sorted(initial_val_dict.items(), key=sort_key)) + '\n'
    initial_val_str += ', '.join('1' if val else '0' for _, val in sorted(initial_val_dict.items(), key=sort_key)) + '\n'

    return model_str, symbol_str, initial_val_str


def write_boolnet(excel_filename: str, boolnet_model_filename: Optional[str]=None, boolnet_symbol_filename: Optional[str]=None,
                  boolnet_initial_val_filename: Optional[str]=None):
    base_name = os.path.splitext(os.path.basename(excel_filename))[0]
    base_path = os.path.dirname(excel_filename)

    if not boolnet_model_filename:
        boolnet_model_filename = os.path.join(base_path, '{0}.boolnet'.format(base_name))
    if not boolnet_symbol_filename:
        boolnet_symbol_filename = os.path.join(base_path, '{0}_symbols.csv'.format(base_name))
    if not boolnet_initial_val_filename:
        boolnet_initial_val_filename = os.path.join(base_path, '{0}_initial_vals.csv'.format(base_name))

    excel_book = ExcelBook(excel_filename)
    model_str, symbol_str, initial_val_str = boolnet_strs_from_rxncon(excel_book.rxncon_system)

    with open(boolnet_model_filename, mode='w') as f:
        f.write(model_str)

    with open(boolnet_symbol_filename, mode='w') as f:
        f.write(symbol_str)

    with open(boolnet_initial_val_filename, mode='w') as f:
        f.write(initial_val_str)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SyntaxError('Please provide the name of the Excel file as an argument.')
    excel_file = sys.argv[1]
    write_boolnet(excel_file)
