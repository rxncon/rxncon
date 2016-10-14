import os

from rxncon.input.excel_book.excel_book import ExcelBook

CELL_CYCLE_XLS = os.path.join(os.path.dirname(__file__), 'cell_cycle_toy_model.xls')
PHEROMONE_XLS  = os.path.join(os.path.dirname(__file__), 'pheromone.xls')


# def test_cell_cycle_toy_model():
#     book = ExcelBook(CELL_CYCLE_XLS)
#     print('hallo')


def test_pheromone():
    book = ExcelBook(PHEROMONE_XLS)
    system = book.rxncon_system
    print('hallo')
