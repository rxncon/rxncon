"""Module containing the class ExcelBook, which is used to read in a rxncon system from an Excel file."""

from typing import List, Optional
import os.path
import xlrd
import logging

from rxncon.core.rxncon_system import RxnConSystem
from rxncon.input.shared.contingency_list import contingencies_from_contingency_list_entries, \
    contingency_list_entry_from_strs, ContingencyListEntry
from rxncon.input.shared.reaction_preprocess import split_bidirectional_reaction_str
from rxncon.core.reaction import Reaction, reaction_from_str, OutputReaction, initialize_reaction_defs
from rxncon.core.state import initialize_state_modifiers
from rxncon.core.contingency import Contingency

NOT_APPLICABLE = 'N/A'

SHEET_COMPONENT_LIST = 'ComponentList'
SHEET_REACTION_TYPE_LIST = 'ReactionTypeDefinition'
SHEET_REACTION_LIST = 'ReactionList'
SHEET_CONTINGENCY_LIST = 'ContingencyList'
SHEET_MODIFICATION_TYPE_LIST = 'ModificationTypeDefinition'

NECESSARY_SHEETS = [SHEET_REACTION_LIST, SHEET_CONTINGENCY_LIST]

HEADER_ROW = 1
DATA_ROW = 2

REACTION_LIST_COLUMN_FULL_NAME = '!UID:Reaction'

CONTINGENCY_LIST_COLUMN_TARGET = '!Target'
CONTINGENCY_LIST_COLUMN_TYPE = '!Contingency'
CONTINGENCY_LIST_COLUMN_MODIFIER = '!Modifier'

MODIFICATION_TYPE_LIST_COLUMN_TYPE = '!UID:ModificationType'
MODIFICATION_TYPE_LIST_COLUMN_LABEL = '!UID:ModificationLabel'

RXN_DEF_LIST_COLUMN_REACTION = '!UID:Reaction'
RXN_DEF_LIST_COLUMN_REACTION_KEY = '!UID:ReactionKey'
RXN_DEF_LIST_COLUMN_BIDIRECTIONAL_VERB = '!BidirectionalVerb'
RXN_DEF_LIST_COLUMN_MOLTYPE_X = '!MolTypeX'
RXN_DEF_LIST_COLUMN_RESOLUTION_X = '!ResolutionX'
RXN_DEF_LIST_COLUMN_MOLTYPE_Y = '!MolTypeY'
RXN_DEF_LIST_COLUMN_RESOLUTION_Y = '!ResolutionY'
RXN_DEF_LIST_COLUMN_RULE = '!SkeletonRule'

logger = logging.getLogger(__name__)


class ExcelBook:
    def __init__(self, filename: str) -> None:
        self.filename = filename
        self._xlrd_book = None  # type: Optional[xlrd.Book]
        self._reactions = []  # type: List[Reaction]
        self._cont_list_entries = []  # type: List[ContingencyListEntry]
        self._contingencies = []  # type: List[Contingency]
        self._rxncon_system = None  # type: Optional[RxnConSystem]

        self._column_modification_type = None  # type: Optional[int]
        self._column_modification_label = None  # type: Optional[int]
        self._column_reaction_full_name = None  # type: Optional[int]
        self._column_contingency_target = None  # type: Optional[int]
        self._column_contingency_type = None  # type: Optional[int]
        self._column_contingency_modifier = None  # type: Optional[int]
        self._column_rxn_def_reaction = None  # type: Optional[int]
        self._column_rxn_def_reaction_key = None  # type: Optional[int]
        self._column_rxn_def_bi_verb = None  # type: Optional[int]
        self._column_rxn_def_moltype_x = None  # type: Optional[int]
        self._column_rxn_def_resolution_x = None  # type: Optional[int]
        self._column_rxn_def_moltype_y = None  # type: Optional[int]
        self._column_rxn_def_resolution_y = None  # type: Optional[int]
        self._column_rule = None  # type: Optional[int]

        self._open_file()
        self._validate_book()
        self._determine_column_numbers()
        self._initialize_modification_types()
        self._initialize_reaction_defs()
        self._load_reaction_list()
        self._load_contingency_list_entries()
        self._construct_output_reactions()
        self._construct_contingencies()
        self._construct_rxncon_system()

    @property
    def rxncon_system(self) -> RxnConSystem:
        assert self._rxncon_system is not None, 'Could not construct rxncon system!'
        return self._rxncon_system

    def _open_file(self) -> None:
        if not os.path.isfile(self.filename):
            raise IOError('Could not find file {}'.format(self.filename))

        self._xlrd_book = xlrd.open_workbook(self.filename)

    def _validate_book(self) -> None:
        assert isinstance(self._xlrd_book, xlrd.Book), 'Could not instantiate Excel workbook!'
        if not all(sheet in self._xlrd_book.sheet_names() for sheet in NECESSARY_SHEETS):  # type: ignore
            raise SyntaxError('Excel book does not contain expected sheets')

    def _determine_column_numbers(self) -> None:
        sheet = self._xlrd_book.sheet_by_name(SHEET_REACTION_LIST)  # type: ignore
        row = list(sheet.get_rows())[HEADER_ROW]
        for num, header in enumerate(row):
            if header.value == REACTION_LIST_COLUMN_FULL_NAME:
                self._column_reaction_full_name = num

        sheet = self._xlrd_book.sheet_by_name(SHEET_CONTINGENCY_LIST)  # type: ignore
        row = list(sheet.get_rows())[HEADER_ROW]
        for num, header in enumerate(row):
            if header.value == CONTINGENCY_LIST_COLUMN_TARGET:
                self._column_contingency_target = num
            elif header.value == CONTINGENCY_LIST_COLUMN_TYPE:
                self._column_contingency_type = num
            elif header.value == CONTINGENCY_LIST_COLUMN_MODIFIER:
                self._column_contingency_modifier = num

        for col_num in (self._column_reaction_full_name, self._column_contingency_target,
                        self._column_contingency_type, self._column_contingency_modifier):
            assert col_num is not None, 'Missing column!'

        try:
            sheet = self._xlrd_book.sheet_by_name(SHEET_MODIFICATION_TYPE_LIST)  # type: ignore
            row = list(sheet.get_rows())[HEADER_ROW]
            for num, header in enumerate(row):
                if header.value == MODIFICATION_TYPE_LIST_COLUMN_TYPE:
                    self._column_modification_type = num
                elif header.value == MODIFICATION_TYPE_LIST_COLUMN_LABEL:
                    self._column_modification_label = num
        except xlrd.XLRDError:
            # Skip empty rows.
            pass

        try:
            sheet = self._xlrd_book.sheet_by_name(SHEET_REACTION_TYPE_LIST)  # type: ignore
            row = list(sheet.get_rows())[HEADER_ROW]
            for num, header in enumerate(row):
                if header.value == RXN_DEF_LIST_COLUMN_REACTION:
                    self._column_rxn_def_reaction = num
                elif header.value == RXN_DEF_LIST_COLUMN_REACTION_KEY:
                    self._column_rxn_def_reaction_key = num
                elif header.value == RXN_DEF_LIST_COLUMN_BIDIRECTIONAL_VERB:
                    self._column_rxn_def_bi_verb = num
                elif header.value == RXN_DEF_LIST_COLUMN_MOLTYPE_X:
                    self._column_rxn_def_moltype_x = num
                elif header.value == RXN_DEF_LIST_COLUMN_RESOLUTION_X:
                    self._column_rxn_def_resolution_x = num
                elif header.value == RXN_DEF_LIST_COLUMN_MOLTYPE_Y:
                    self._column_rxn_def_moltype_y = num
                elif header.value == RXN_DEF_LIST_COLUMN_RESOLUTION_Y:
                    self._column_rxn_def_resolution_y = num
                elif header.value == RXN_DEF_LIST_COLUMN_RULE:
                    self._column_rule = num
        except xlrd.XLRDError:
            # Skip empty rows.
            pass

    def _initialize_modification_types(self) -> None:
        if self._column_modification_type is None or self._column_modification_label is None:
            return

        sheet = self._xlrd_book.sheet_by_name(SHEET_MODIFICATION_TYPE_LIST)  # type: ignore
        modifier_rows = [row for row in sheet.get_rows()][DATA_ROW:]

        modifiers = {}

        for row in modifier_rows:
            k = row[self._column_modification_type].value
            v = row[self._column_modification_label].value

            if not isinstance(k, str):
                raise SyntaxError('Modification type {} needs to be a str.'.format(k))

            if isinstance(v, float) or isinstance(v, int):
                v = str(int(v))

            if not isinstance(v, str):
                raise SyntaxError('Modification label {} needs to be str or number.'.format(v))

            modifiers[k] = v

        # The new modifiers are loaded using this function.
        initialize_state_modifiers(modifiers)

    def _initialize_reaction_defs(self) -> None:
        if self._column_rxn_def_reaction is None:
            return

        sheet = self._xlrd_book.sheet_by_name(SHEET_REACTION_TYPE_LIST)  # type: ignore
        rxn_def_rows = [row for row in sheet.get_rows()][DATA_ROW:]

        rxn_defs = []

        for row in rxn_def_rows:
            rxn_defs.append({
                RXN_DEF_LIST_COLUMN_REACTION: row[self._column_rxn_def_reaction].value,
                RXN_DEF_LIST_COLUMN_REACTION_KEY: row[self._column_rxn_def_reaction_key].value,
                RXN_DEF_LIST_COLUMN_BIDIRECTIONAL_VERB: row[self._column_rxn_def_bi_verb].value,
                RXN_DEF_LIST_COLUMN_MOLTYPE_X: row[self._column_rxn_def_moltype_x].value,
                RXN_DEF_LIST_COLUMN_RESOLUTION_X: row[self._column_rxn_def_resolution_x].value,
                RXN_DEF_LIST_COLUMN_MOLTYPE_Y: row[self._column_rxn_def_moltype_y].value,
                RXN_DEF_LIST_COLUMN_RESOLUTION_Y: row[self._column_rxn_def_resolution_y].value,
                RXN_DEF_LIST_COLUMN_RULE: row[self._column_rule].value
            })

        initialize_reaction_defs(rxn_defs)

    def _load_reaction_list(self) -> None:
        sheet = self._xlrd_book.sheet_by_name(SHEET_REACTION_LIST)  # type: ignore
        reaction_rows = [row for row in sheet.get_rows()][DATA_ROW:]

        for row in reaction_rows:
            if not row[self._column_reaction_full_name].value:
                logger.debug('_load_reaction_list: Empty row')
                continue

            logger.debug('_load_reaction_list: {}'.format(row[self._column_reaction_full_name].value))

            # When a verb such as 'ppi' is encountered, the function 'preprocessed_reaction_strs'
            # will split it into 'ppi+' and 'ppi-'.
            reaction_strs = split_bidirectional_reaction_str(row[self._column_reaction_full_name].value)
            self._reactions += [reaction_from_str(x) for x in reaction_strs]

    def _load_contingency_list_entries(self) -> None:
        sheet = self._xlrd_book.sheet_by_name(SHEET_CONTINGENCY_LIST)  # type: ignore
        contingency_rows = [row for row in sheet.get_rows()][DATA_ROW:]

        for row in contingency_rows:
            if not row[self._column_contingency_target].value:
                continue

            # When a reaction with a verb such as 'ppi' is encountered, we only apply the contingency to the
            # positive reaction.
            target = row[self._column_contingency_target].value
            if split_bidirectional_reaction_str(target) != [target]:
                target = split_bidirectional_reaction_str(target)[0]
                logger.info('_load_contingency_list_entries: applying contingency {} to forward target {}'
                            .format(row[self._column_contingency_modifier].value, target))

            entry = contingency_list_entry_from_strs(
                target, row[self._column_contingency_type].value,
                row[self._column_contingency_modifier].value
            )
            self._cont_list_entries.append(entry)

    def _construct_output_reactions(self) -> None:
        self._reactions += [cle.subj for cle in self._cont_list_entries
                            if isinstance(cle.subj, OutputReaction)]

    def _construct_contingencies(self) -> None:
        self._contingencies = contingencies_from_contingency_list_entries(self._cont_list_entries)

    def _construct_rxncon_system(self) -> None:
        self._rxncon_system = RxnConSystem(self._reactions, self._contingencies)
