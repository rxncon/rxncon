from typing import List
from abc import abstractmethod

import os.path

import xlrd

import rxncon.core.error as err
import rxncon.core.rxncon_system as rxs
import rxncon.input.rxncon_input as inp
import rxncon.input.shared.contingency_list as cli
import rxncon.syntax.rxncon_from_string as fst

NOT_APPLICABLE = 'N/A'

class ExcelBook(inp.RxnConInput):

    def __init__(self):
        raise AssertionError

    @property
    def rxncon_system(self) -> rxs.RxnConSystem:
        return self._rxncon_system

    def _open_file(self):
        if not os.path.isfile(self.filename):
            raise err.RxnConInputError('Could not find file {}'.format(self.filename))

        self._xlrd_book = xlrd.open_workbook(self.filename)

    @abstractmethod
    def _validate_book(self):
        pass

    @abstractmethod
    def _load_reaction_list(self):
        pass

    @abstractmethod
    def _load_contingency_list_entries(self):
        pass

    def _construct_contingencies(self):
        self._contingencies = cli.contingencies_from_contingency_list_entries(self._contingency_list_entries)

    def _construct_rxncon_system(self):
        self._rxncon_system = rxs.RxnConSystem(self._reactions, self._contingencies)


class ExcelBookWithoutReactionType(ExcelBook):

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

    def __init__(self, filename: str):
        self.filename = filename
        self._xlrd_book = None               # type: xlrd.Book
        self._reactions = []                 # type: List[rxn.Reaction]
        self._contingency_list_entries = []  # type: List[cli.ContingencyListEntry]
        self._contingencies = []             # type: List[con.Contingency]
        self._rxncon_system = None           # type: rxs.RxnConSystem

        super()._open_file()
        self._validate_book()
        self._load_reaction_list()
        self._load_contingency_list_entries()
        super()._construct_contingencies()
        super()._construct_rxncon_system()

    def _validate_book(self):
        expected_sheets = [self.SHEET_REACTION_LIST, self.SHEET_METABOLIC_REACTION_LIST,
                           self.SHEET_CONTINGENCY_LIST, self.SHEET_REACTION_DEFINITION]

        if not all([sheet in self._xlrd_book.sheet_names() for sheet in expected_sheets]):
            raise err.RxnConParseError('Excel book does not contain expected sheets')

    def _load_reaction_list(self):
        sheet = self._xlrd_book.sheet_by_name(self.SHEET_REACTION_LIST)
        reaction_rows = [row for row in sheet.get_rows()][1:]

        for row in reaction_rows:
            reaction = fst.reaction_from_string(row[self.REACTION_LIST_COLUMN_FULL_NAME].value)
            self._reactions.append(reaction)

    def _load_contingency_list_entries(self):
        sheet = self._xlrd_book.sheet_by_name(self.SHEET_CONTINGENCY_LIST)
        contingency_rows = [row for row in sheet.get_rows()][1:]

        for row in contingency_rows:
            entry = cli.contingency_list_entry_from_subject_predicate_agent_strings(
                row[self.CONTINGENCY_LIST_COLUMN_TARGET].value, row[self.CONTINGENCY_LIST_COLUMN_TYPE].value,
                row[self.CONTINGENCY_LIST_COLUMN_MODIFIER].value
            )
            self._contingency_list_entries.append(entry)


class ExcelBookWithReactionType(ExcelBook):

    SHEET_COMPONENT_LIST = 'ComponentList'
    SHEET_REACTION_DEFINITION = 'ReactionDefinition'
    SHEET_REACTION_LIST = 'ReactionList'
    SHEET_CONTINGENCY_LIST = 'ContingencyList'

    REACTION_LIST_COLUMN_UID = 0
    REACTION_LIST_COLUMN_SOURCE_STATE = 1
    REACTION_LIST_COLUMN_PRODUCT_STATE = 2

    CONTINGENCY_LIST_COLUMN_TARGET = 1
    CONTINGENCY_LIST_COLUMN_TYPE = 2
    CONTINGENCY_LIST_COLUMN_MODIFIER = 3


    def __init__(self, filename: str):
        self.filename = filename
        self._xlrd_book = None               # type: xlrd.Book
        self._reactions = []                 # type: List[rxn.Reaction]
        self._contingency_list_entries = []  # type: List[cli.ContingencyListEntry]
        self._contingencies = []             # type: List[con.Contingency]
        self._rxncon_system = None           # type: rxs.RxnConSystem

        super()._open_file()
        self._validate_book()
        self._load_reaction_list()
        self._load_contingency_list_entries()
        super()._construct_contingencies()
        super()._construct_rxncon_system()

    def _validate_book(self):
        expected_sheets = [self.SHEET_COMPONENT_LIST, self.SHEET_REACTION_DEFINITION,
                           self.SHEET_REACTION_LIST, self.SHEET_CONTINGENCY_LIST]

        if not all([sheet in self._xlrd_book.sheet_names() for sheet in expected_sheets]):
            raise err.RxnConParseError('Excel book does not contain expected sheets')

    def _load_reaction_list(self):
        sheet = self._xlrd_book.sheet_by_name(self.SHEET_REACTION_LIST)
        reaction_rows = [row for row in sheet.get_rows()][2:]

        for row in reaction_rows:
            reaction = fst.reaction_from_string(row[self.REACTION_LIST_COLUMN_UID].value)
            self._reactions.append(reaction)

    def _load_contingency_list_entries(self):
        sheet = self._xlrd_book.sheet_by_name(self.SHEET_CONTINGENCY_LIST)
        contingency_rows = [row for row in sheet.get_rows()][2:]

        for row in contingency_rows:
            entry = cli.contingency_list_entry_from_subject_predicate_agent_strings(
                row[self.CONTINGENCY_LIST_COLUMN_TARGET].value, row[self.CONTINGENCY_LIST_COLUMN_TYPE].value,
                row[self.CONTINGENCY_LIST_COLUMN_MODIFIER].value
            )
            self._contingency_list_entries.append(entry)

    def _load_reaction_definition(self):
        pass