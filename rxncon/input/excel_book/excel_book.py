from typing import List
import os.path
import xlrd

from rxncon.core.rxncon_system import RxnConSystem
from rxncon.input.shared.contingency_list import contingencies_from_contingency_list_entries, \
    contingency_list_entry_from_strs, ContingencyListEntry
from rxncon.input.shared.reaction_preprocess import split_bidirectional_reaction_str
from rxncon.core.reaction import Reaction, reaction_from_str, OutputReaction
from rxncon.core.contingency import Contingency


NOT_APPLICABLE = 'N/A'

SHEET_COMPONENT_LIST      = 'ComponentList'
SHEET_REACTION_DEFINITION = 'ReactionDefinition'
SHEET_REACTION_LIST       = 'ReactionList'
SHEET_CONTINGENCY_LIST    = 'ContingencyList'

NECESSARY_SHEETS = [SHEET_REACTION_LIST, SHEET_CONTINGENCY_LIST]

HEADER_ROW = 1
DATA_ROW   = 2

REACTION_LIST_COLUMN_FULL_NAME = '!UID:Reaction'

CONTINGENCY_LIST_COLUMN_TARGET   = '!Target'
CONTINGENCY_LIST_COLUMN_TYPE     = '!Contingency'
CONTINGENCY_LIST_COLUMN_MODIFIER = '!Modifier'


class ExcelBook:
    def __init__(self, filename: str):
        self.filename           = filename
        self._xlrd_book         = None   # type: xlrd.Book
        self._reactions         = []     # type: List[Reaction]
        self._cont_list_entries = []     # type: List[ContingencyListEntry]
        self._contingencies     = []     # type: List[Contingency]
        self._rxncon_system     = None   # type: RxnConSystem

        self._column_reaction_full_name   = None
        self._column_contingency_target   = None
        self._column_contingency_type     = None
        self._column_contingency_modifier = None

        self._open_file()
        self._validate_book()
        self._determine_column_numbers()
        self._load_reaction_list()
        self._load_contingency_list_entries()
        self._construct_output_reactions()
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
        if not all(sheet in self._xlrd_book.sheet_names() for sheet in NECESSARY_SHEETS):
            raise SyntaxError('Excel book does not contain expected sheets')

    def _determine_column_numbers(self):
        sheet = self._xlrd_book.sheet_by_name(SHEET_REACTION_LIST)
        row   = list(sheet.get_rows())[HEADER_ROW]
        for num, header in enumerate(row):
            if header.value == REACTION_LIST_COLUMN_FULL_NAME:
                self._column_reaction_full_name = num

        sheet = self._xlrd_book.sheet_by_name(SHEET_CONTINGENCY_LIST)
        row   = list(sheet.get_rows())[HEADER_ROW]
        for num, header in enumerate(row):
            if header.value == CONTINGENCY_LIST_COLUMN_TARGET:
                self._column_contingency_target = num
            elif header.value == CONTINGENCY_LIST_COLUMN_TYPE:
                self._column_contingency_type = num
            elif header.value == CONTINGENCY_LIST_COLUMN_MODIFIER:
                self._column_contingency_modifier = num

        for col_num in (self._column_reaction_full_name, self._column_contingency_target,
                        self._column_contingency_type, self._column_contingency_modifier):
            assert col_num is not None

    def _load_reaction_list(self):
        sheet = self._xlrd_book.sheet_by_name(SHEET_REACTION_LIST)
        reaction_rows = [row for row in sheet.get_rows()][DATA_ROW:]

        for row in reaction_rows:
            if not row[self._column_reaction_full_name].value:
                continue

            # When a verb such as 'ppi' is encountered, the function 'preprocessed_reaction_strs'
            # will split it into 'ppi+' and 'ppi-'.
            reaction_strs = split_bidirectional_reaction_str(row[self._column_reaction_full_name].value)
            self._reactions += [reaction_from_str(x) for x in reaction_strs]

    def _load_contingency_list_entries(self):
        sheet = self._xlrd_book.sheet_by_name(SHEET_CONTINGENCY_LIST)
        contingency_rows = [row for row in sheet.get_rows()][DATA_ROW:]

        for row in contingency_rows:
            if not row[self._column_contingency_target].value:
                continue

            # When a reaction with a verb such as 'ppi' is encountered, we only apply the contingency to the
            # positive reaction.
            target = row[self._column_contingency_target].value
            if split_bidirectional_reaction_str(target) != [target]:
                target = split_bidirectional_reaction_str(target)[0]

            entry = contingency_list_entry_from_strs(
                target, row[self._column_contingency_type].value,
                row[self._column_contingency_modifier].value
            )
            self._cont_list_entries.append(entry)

    def _construct_output_reactions(self):
        self._reactions += [cle.subject for cle in self._cont_list_entries
                            if isinstance(cle.subject, OutputReaction)]

    def _construct_contingencies(self):
        self._contingencies = contingencies_from_contingency_list_entries(self._cont_list_entries)

    def _construct_rxncon_system(self):
        self._rxncon_system = RxnConSystem(self._reactions, self._contingencies)
