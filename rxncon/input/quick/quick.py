"""Module containing the class Quick, a text-based format to read a rxncon system."""

import re
from typing import List, Optional

from rxncon.input.shared.contingency_list import contingencies_from_contingency_list_entries, \
    contingency_list_entry_from_strs, ContingencyListEntry
from rxncon.core.reaction import reaction_from_str
from rxncon.core.rxncon_system import RxnConSystem
from rxncon.core.reaction import Reaction
from rxncon.core.contingency import Contingency
from rxncon.input.shared.reaction_preprocess import split_bidirectional_reaction_str


class Quick:
    def __init__(self, rxncon_str: str) -> None:
        self.quick_input = rxncon_str.split('\n')
        self._rxncon_system = None              # type: Optional[RxnConSystem]
        self._reactions = []                    # type: List[Reaction]
        self._contingencies = []                # type: List[Contingency]
        self._contingency_list_entries = []     # type: List[ContingencyListEntry]
        self._parse_str()
        self._construct_contingencies()
        self._construct_rxncon_system()
        assert self._rxncon_system is not None

    @property
    def rxncon_system(self) -> RxnConSystem:
        assert self._rxncon_system is not None
        return self._rxncon_system

    def _parse_str(self) -> None:
        BOOL_REGEX = '^\<.+?\>$'

        for line in self.quick_input:
            reaction_string = line.split(';')[0].strip()
            contingency_strings = line.split(';')[1:]
            if reaction_string:
                if not re.match(BOOL_REGEX, reaction_string):
                    self._add_reaction_from_string(reaction_string)
                self._add_contingency_list_entries(contingency_strings, reaction_string)

    def _add_reaction_from_string(self, reaction_str: str) -> None:
        reaction_strs = split_bidirectional_reaction_str(reaction_str)
        for rxn in reaction_strs:
            reaction = reaction_from_str(rxn)
            self._reactions.append(reaction)

    def _add_contingency_list_entries(self, contingency_strs: List[str], reaction_str: str) -> None:
        for cont in contingency_strs:
            cont = cont.strip()
            cont_type = cont.split()[0]
            modifier = cont.split()[-1]

            #  If the verb is bidirectional, only apply the contingency to the forward direction.
            reaction_strs = split_bidirectional_reaction_str(reaction_str)
            entry = contingency_list_entry_from_strs(reaction_strs[0], cont_type, modifier)
            self._contingency_list_entries.append(entry)

    def _construct_contingencies(self) -> None:
        self._contingencies = contingencies_from_contingency_list_entries(self._contingency_list_entries)

    def _construct_rxncon_system(self) -> None:
        self._rxncon_system = RxnConSystem(self._reactions, self._contingencies)