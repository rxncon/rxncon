import re
from typing import List
from typecheck import typecheck

from rxncon.input.shared.contingency_list import contingencies_from_contingency_list_entries, \
    contingency_list_entry_from_subject_predicate_agent_strings
from rxncon.core.reaction import reaction_from_string, REACTION_DEFINITIONS
from rxncon.core.rxncon_system import RxnConSystem
from rxncon.core.reaction import Reaction
from rxncon.core.contingency import Contingency
from rxncon.input.shared.contingency_list import ContingencyListEntry

class Quick:
    @typecheck
    def __init__(self, rxncon_str: str):
        self.quick_input = rxncon_str.split('\n')
        self.rxncon_system = None               # type: RxnConSystem
        self._reactions = []                    # type: List[Reaction]
        self._contingencies = []                # type: List[Contingency]
        self._contingency_list_entries = []     # type: List[ContingencyListEntry]

        self._parse_str()
        self._construct_contingencies()
        self._construct_rxncon_system()

    def _parse_str(self):
        INPUT_REGEX = '^\[.+?\]$'
        BOOL_REGEX = '^\<.+?\>$'

        for line in self.quick_input:
            reaction_string = line.split(';')[0].strip()
            contingency_strings = line.split(';')[1:]
            if reaction_string:
                if not re.match(BOOL_REGEX, reaction_string) and not re.match(INPUT_REGEX, reaction_string):
                    self._add_reaction_from_string(reaction_string)
                self._add_contingency_list_entries(contingency_strings, reaction_string)

    @typecheck
    def _add_reaction_from_string(self, full_reaction_str: str):
        reaction = reaction_from_string(REACTION_DEFINITIONS, full_reaction_str)
        self._reactions.append(reaction)

    @typecheck
    def _add_contingency_list_entries(self, contingency_strings: List[str], reaction_string: str):
        for cont in contingency_strings:
            cont = cont.strip()
            cont_type = cont.split()[0]
            modifier = cont.split()[-1]
            entry = contingency_list_entry_from_subject_predicate_agent_strings(reaction_string, cont_type, modifier)
            self._contingency_list_entries.append(entry)

    def _construct_contingencies(self):
        self._contingencies = contingencies_from_contingency_list_entries(self._contingency_list_entries)

    def _construct_rxncon_system(self):
        self.rxncon_system = RxnConSystem(self._reactions, self._contingencies)
