import os
import pytest

from rxncon.input.excel_book.excel_book import ExcelBook
from rxncon.core.reaction import reaction_from_str, initialize_reaction_defs
from rxncon.core.state import state_from_str, initialize_state_modifiers


MISSING_NECESSARY_SHEET_XLS   = os.path.join(os.path.dirname(__file__), 'missing_necessary_sheet.xls')
MISSING_UNNECESSARY_SHEET_XLS = os.path.join(os.path.dirname(__file__), 'missing_unnecessary_sheet.xls')
SHUFFLED_COLUMNS_XLS          = os.path.join(os.path.dirname(__file__), 'shuffled_columns.xls')
ADDITIONAL_MODIFIERS_XLS      = os.path.join(os.path.dirname(__file__), 'additional_modifiers.xls')
ADDITIONAL_REACTIONS_XLS      = os.path.join(os.path.dirname(__file__), 'additional_rxns.xls')


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


def test_additional_modifiers() -> None:
    with pytest.raises(ValueError):
        state_from_str('A_[(x)]-{bladiebla}')

    rxncon_system = ExcelBook(ADDITIONAL_MODIFIERS_XLS).rxncon_system
    assert str(state_from_str('A_[(x)]-{bladiebla}')) == 'A_[(x)]-{bladiebla}'

    initialize_state_modifiers()

    with pytest.raises(ValueError):
        state_from_str('A_[(x)]-{bladiebla}')


def test_additional_reactions() -> None:
    with pytest.raises(SyntaxError):
        reaction_from_str('A_smurf+_B_[(r)]')

    rxncon_system = ExcelBook(ADDITIONAL_REACTIONS_XLS).rxncon_system
    rxn = reaction_from_str('A_smurf+_B_[(r)]')

    assert rxn.produced_states == [state_from_str('B_[(r)]-{smurf}')]
    assert rxn.consumed_states == [state_from_str('B_[(r)]-{0}')]

    initialize_reaction_defs()

    with pytest.raises(SyntaxError):
        reaction_from_str('A_smurf+_B_[(r)]')
