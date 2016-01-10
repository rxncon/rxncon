from typing import List
import os.path

import xlrd

import rxncon.core.error as err
import rxncon.core.rxncon_system as rxs
import rxncon.input.rxncon_input as inp
import rxncon.input.shared.contingency_list as cli
import rxncon.syntax.rxncon_from_string as fst

NOT_APPLICABLE = 'N/A'

SHEET_REACTION_LIST = '(I) Reaction list'
SHEET_METABOLIC_REACTION_LIST = '(II) Metabolic reaction list'
SHEET_CONTINGENCY_LIST = '(III) Contingency list'
SHEET_REACTION_DEFINITION = '(IV) Reaction definition'

REACTION_LIST_COLUMN_FULL_NAME = 1
REACTION_LIST_COLUMN_SOURCE_STATE = 2
REACTION_LIST_COLUMN_PRODUCT_STATE = 3

CONTINGENCY_LIST_COLUMN_TARGET = 1
CONTINGENCY_LIST_COLUMN_TYPE = 2
CONTINGENCY_LIST_COLUMN_MODIFIER = 3


class ExcelBook(inp.RxnConInput):
    def __init__(self, filename: str):
        self.filename = filename
        self._xlrd_book = None               # type: xlrd.Book
        self._reactions = []                 # type: List[rxn.Reaction]
        self._contingency_list_entries = []  # type: List[cli.ContingencyListEntry]
        self._contingencies = []             # type: List[con.Contingency]
        self._rxncon_system = None           # type: rxs.RxnConSystem

        self._open_file()
        self._validate_book()
        self._load_reaction_list()
        self._load_contingency_list_entries()
        self._construct_contingencies()
        self._construct_rxncon_system()

    @property
    def rxncon_system(self) -> rxs.RxnConSystem:
        return self._rxncon_system

    def _open_file(self):
        if not os.path.isfile(self.filename):
            raise err.RxnConInputError('Could not find file {}'.format(self.filename))

        self._xlrd_book = xlrd.open_workbook(self.filename)

    def _validate_book(self):
        expected_sheets = [SHEET_REACTION_LIST, SHEET_METABOLIC_REACTION_LIST, SHEET_CONTINGENCY_LIST, SHEET_REACTION_DEFINITION]

        if not all([sheet in self._xlrd_book.sheet_names() for sheet in expected_sheets]):
            raise err.RxnConParseError('Excel book does not contain expected sheets')

    def _load_reaction_list(self):
        sheet = self._xlrd_book.sheet_by_name(SHEET_REACTION_LIST)
        reaction_rows = [row for row in sheet.get_rows()][1:]

        for row in reaction_rows:
            reaction = fst.reaction_from_string(row[REACTION_LIST_COLUMN_FULL_NAME].value)
            self._reactions.append(reaction)

    def _load_contingency_list_entries(self):
        sheet = self._xlrd_book.sheet_by_name(SHEET_CONTINGENCY_LIST)
        contingency_rows = [row for row in sheet.get_rows()][1:]

        for row in contingency_rows:
            entry = cli.contingency_list_entry_from_subject_predicate_agent_strings(
                row[CONTINGENCY_LIST_COLUMN_TARGET].value, row[CONTINGENCY_LIST_COLUMN_TYPE].value, row[CONTINGENCY_LIST_COLUMN_MODIFIER].value
            )
            self._contingency_list_entries.append(entry)

    def _construct_contingencies(self):
        self._contingencies = cli.contingencies_from_contingency_list_entries(self._contingency_list_entries)

    def _construct_rxncon_system(self):
        pass

