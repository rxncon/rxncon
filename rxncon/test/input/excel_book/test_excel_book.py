import os

from rxncon.input.excel_book.excel_book import ExcelBook
from rxncon.core.reaction import reaction_from_str

CELL_CYCLE_XLS = os.path.join(os.path.dirname(__file__), 'cell_cycle_toy_model.xls')
PHEROMONE_XLS  = os.path.join(os.path.dirname(__file__), 'pheromone_structured.xls')
SPS_XLS        = os.path.join(os.path.dirname(__file__), 'sps.xls')


def test_pheromone():
    book = ExcelBook(PHEROMONE_XLS)
    system = book.rxncon_system
    cs = system.contingencies_for_reaction(reaction_from_str('Ste4_[Ste5]_ppi+_Ste5_[Ste4]'))

    for c in cs:
        print()
        print(c)

