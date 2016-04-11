import os

import rxncon.input.excel_book.excel_book as exc

TIGER_FILENAME = 'Tiger_et_al_TableS1.xls'
TIGER_PATH = os.path.join(os.path.dirname(__file__), TIGER_FILENAME)

RXNCON_NEW_TEMPLATE_FILENAME = 'rxncon_template_2_0_no_relocalisation.xls'
RXNCON_NEW_TEMPLATE_PATH = os.path.join(os.path.dirname(__file__), RXNCON_NEW_TEMPLATE_FILENAME)

def test_excel_book_tiger_network():
    excel_book = exc.ExcelBookWithoutReactionType(TIGER_PATH)

    assert isinstance(excel_book, exc.ExcelBook)
    assert isinstance(excel_book, exc.ExcelBookWithoutReactionType)
    assert len(excel_book._reactions) == 222
    assert len(excel_book._contingencies) == 279

    print(len(excel_book.rxncon_system.reactions))


def test_new_excel_book_input():
    excel_book = exc.ExcelBookWithReactionType(RXNCON_NEW_TEMPLATE_PATH)

    assert isinstance(excel_book, exc.ExcelBook)
    assert isinstance(excel_book, exc.ExcelBookWithReactionType)
    assert len(excel_book._reactions) == 16
    assert len(excel_book._contingencies) == 2
