from typing import List
import typecheck as tc
import re

import rxncon.input.rxncon_input as inp
import rxncon.core.rxncon_system as rxs
import rxncon.input.shared.contingency_list as cli
import rxncon.syntax.rxncon_from_string as fst


class Quick(inp.RxnConInput):
    @tc.typecheck
    def __init__(self, rxncon_string: str):
        self.quick_input = rxncon_string.split("\n")
        self._reactions = []                    # type: List[rxn.Reaction]
        self._contingencies = []                # type: List[con.Contingency]
        self._contingency_list_entries = []     # type: List[cli.ContingencyListEntry]
        self._rxncon_system = None              # type: rxs.RxnConSystem

        self._parse_str()
        self._construct_contingencies()
        self._construct_rxncon_system()

    @property
    def rxncon_system(self) -> rxs.RxnConSystem:
        return self._rxncon_system

    def _parse_str(self):
        INPUT_REGEX = '^\[.+?\]$'
        BOOL_REGEX = '^\<.+?\>$'
        for line in self.quick_input:
            reaction_string = line.split(";")[0].strip()
            contingency_strings = line.split(";")[1:]
            if reaction_string:
                if not re.match(BOOL_REGEX, reaction_string) and not re.match(INPUT_REGEX, reaction_string):
                    self._add_reaction_from_string(reaction_string)
                self._add_contingency_list_entries(contingency_strings, reaction_string)

    @tc.typecheck
    def _add_reaction_from_string(self, full_reaction_str: str):
        reaction = fst.reaction_from_string(full_reaction_str)
        self._reactions.append(reaction)

    @tc.typecheck
    def _add_contingency_list_entries(self, contingency_strings: List[str], reaction_string: str):
        for cont in contingency_strings:
            cont = cont.strip()
            cont_type = cont.split()[0]
            modifier = cont.split()[-1]
            entry = cli.contingency_list_entry_from_subject_predicate_agent_strings(reaction_string, cont_type, modifier)
            self._contingency_list_entries.append(entry)

    def _construct_contingencies(self):
        self._contingencies = cli.contingencies_from_contingency_list_entries(self._contingency_list_entries)

    def _construct_rxncon_system(self):
        self._rxncon_system = rxs.RxnConSystem(self._reactions, self._contingencies)
