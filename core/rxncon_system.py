from typing import List

import core.reaction as rxn
import core.contingency as con
import core.component as com


class RxnConSystem:
    def __init__(self, reactions: List[rxn.Reaction], contingencies: List[con.Contingency]):
        self.reactions = reactions
        self.contingencies = contingencies
        self.states = None
        self._generate_state_list()
        self._assert_consistency()

    @property
    def constant_components(self) -> List[com.Component]:
        # @todo implement this?
        return []

    def _generate_state_list(self):
        for reaction in self.reactions:
            if reaction.product and reaction.product not in self.states:
                self.states.append(reaction.product)

            if reaction.source and reaction.source not in self.states:
                self.states.append(reaction.source)

    def _assert_consistency(self):
        pass
