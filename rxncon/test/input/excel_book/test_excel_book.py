import os

from rxncon.input.excel_book.excel_book import ExcelBook
from rxncon.simulation.boolean.boolean_model import boolean_model_from_rxncon


CELL_CYCLE_XLS = os.path.join(os.path.dirname(__file__), 'cell_cycle_toy_model.xls')
PHEROMONE_XLS  = os.path.join(os.path.dirname(__file__), 'pheromone_structured.xls')
SPS_XLS        = os.path.join(os.path.dirname(__file__), 'sps.xls')
CIRCADIAN_XLS  = os.path.join(os.path.dirname(__file__), 'CircadianClock.xls')


def test_pheromone():
    book = ExcelBook(PHEROMONE_XLS)
    system = book.rxncon_system
    print('hallo')
