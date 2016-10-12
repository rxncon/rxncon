from typing import List
import os.path
import xlrd

from rxncon.core.rxncon_system import RxnConSystem
from rxncon.input.shared.contingency_list import contingencies_from_contingency_list_entries, \
    contingency_list_entry_from_subject_verb_object_strings, ContingencyListEntry
from rxncon.input.shared.reaction_preprocess import preprocessed_reaction_strs
from rxncon.core.reaction import Reaction, reaction_from_str
from rxncon.core.contingency import Contingency


NOT_APPLICABLE = 'N/A'

SHEET_COMPONENT_LIST = 'ComponentList'
SHEET_REACTION_DEFINITION = 'ReactionDefinition'
SHEET_REACTION_LIST = 'ReactionList'
SHEET_CONTINGENCY_LIST = 'ContingencyList'

REACTION_LIST_COLUMN_FULL_NAME = 0
REACTION_LIST_FIRST_ROW        = 2

CONTINGENCY_LIST_COLUMN_TARGET   = 1
CONTINGENCY_LIST_COLUMN_TYPE     = 2
CONTINGENCY_LIST_COLUMN_MODIFIER = 3
CONTINGENCY_LIST_FIRST_ROW       = 2

class ExcelBook:
    def __init__(self, filename: str):
        self.filename = filename
        self._xlrd_book = None               # type: xlrd.Book
        self._reactions = []                 # type: List[Reaction]
        self._contingency_list_entries = []  # type: List[ContingencyListEntry]
        self._contingencies = []             # type: List[Contingency]
        self._rxncon_system = None           # type: RxnConSystem

        self._open_file()
        self._validate_book()
        self._load_reaction_list()
        self._load_contingency_list_entries()
        self._construct_contingencies()
        self._construct_rxncon_system()

    @property
    def rxncon_system(self) -> RxnConSystem:
        return self._rxncon_system

    def _open_file(self):
        if not os.path.isfile(self.filename):
            raise IOError('Could not find file {}'.format(self.filename))

        self._xlrd_book = xlrd.open_workbook(self.filename)

    def _validate_book(self):
        expected_sheets = [SHEET_COMPONENT_LIST, SHEET_REACTION_DEFINITION, SHEET_REACTION_LIST, SHEET_CONTINGENCY_LIST]

        if not all([sheet in self._xlrd_book.sheet_names() for sheet in expected_sheets]):
            raise SyntaxError('Excel book does not contain expected sheets')

    def _load_reaction_list(self):
        sheet = self._xlrd_book.sheet_by_name(SHEET_REACTION_LIST)
        reaction_rows = [row for row in sheet.get_rows()][REACTION_LIST_FIRST_ROW:]

        for row in reaction_rows:
            # When a verb such as 'ppi' is encountered, the function 'preprocessed_reaction_strs'
            # will split it into 'ppi+' and 'ppi-'.
            reaction_strs = preprocessed_reaction_strs(row[REACTION_LIST_COLUMN_FULL_NAME].value)
            self._reactions += [reaction_from_str(x) for x in reaction_strs]

    def _load_contingency_list_entries(self):
        sheet = self._xlrd_book.sheet_by_name(SHEET_CONTINGENCY_LIST)
        contingency_rows = [row for row in sheet.get_rows()][CONTINGENCY_LIST_FIRST_ROW:]

        for row in contingency_rows:
            # When a reaction with a verb such as 'ppi' is encountered, we only apply the contingency to the
            # positive reaction.
            target = row[CONTINGENCY_LIST_COLUMN_TARGET].value
            if preprocessed_reaction_strs(target) != [target]:
                target = preprocessed_reaction_strs(target)[0]

            entry = contingency_list_entry_from_subject_verb_object_strings(
                target, row[CONTINGENCY_LIST_COLUMN_TYPE].value,
                row[CONTINGENCY_LIST_COLUMN_MODIFIER].value
            )
            self._contingency_list_entries.append(entry)

    def _construct_contingencies(self):
        self._contingencies = contingencies_from_contingency_list_entries(self._contingency_list_entries)

    def _construct_rxncon_system(self):
        self._rxncon_system = RxnConSystem(self._reactions, self._contingencies)
