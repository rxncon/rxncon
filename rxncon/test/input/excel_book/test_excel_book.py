import os

import rxncon.input.excel_book.excel_book as exc

TIGER_FILENAME = 'Tiger_et_al_TableS1.xls'
TIGER_PATH = os.path.join(os.path.dirname(__file__), TIGER_FILENAME)


# def test_excel_book_tiger_network():
#     excel_book = exc.ExcelBook(TIGER_PATH)
#
#     assert isinstance(excel_book, exc.ExcelBook)
#     assert len(excel_book._reactions) == 222
#     assert len(excel_book._contingencies) == 279
#
#     print(len(excel_book.rxncon_system.reactions))
