import os
import pytest

from rxncon.input.excel_book.excel_book import ExcelBook
from rxncon.core.reaction import reaction_from_str

MISSING_NECESSARY_SHEET_XLS   = os.path.join(os.path.dirname(__file__), 'missing_necessary_sheet.xls')
MISSING_UNNECESSARY_SHEET_XLS = os.path.join(os.path.dirname(__file__), 'missing_unnecessary_sheet.xls')
SHUFFLED_COLUMNS_XLS          = os.path.join(os.path.dirname(__file__), 'shuffled_columns.xls')


def test_missing_necessary_sheet() -> None:
    with pytest.raises(SyntaxError):
        ExcelBook(MISSING_NECESSARY_SHEET_XLS)


def test_missing_unnecessary_sheet() -> None:
    rxncon_system = ExcelBook(MISSING_UNNECESSARY_SHEET_XLS).rxncon_system

    expected_reactions = ['A_[x]_ppi+_B_[y]', 'A_[x]_ppi-_B_[y]', 'C_p+_A_[(z)]']

    for rxn in expected_reactions:
        assert reaction_from_str(rxn) in rxncon_system.reactions


def test_shuffled_columns() -> None:
    rxncon_system = ExcelBook(SHUFFLED_COLUMNS_XLS).rxncon_system

    expected_reactions = ['A_[x]_ppi+_B_[y]', 'A_[x]_ppi-_B_[y]', 'C_p+_A_[(z)]']

    for rxn in expected_reactions:
        assert reaction_from_str(rxn) in rxncon_system.reactions
