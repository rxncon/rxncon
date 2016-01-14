from typing import List

import rxncon.core.component as com
import rxncon.core.contingency as con
import rxncon.core.reaction as rxn


class RxnConSystem:
    def __init__(self, reactions: List[rxn.Reaction], contingencies: List[con.Contingency]):
        self.reactions = reactions
        self.implicit_reactions = []  # type: List[rxn.Reaction]
        self.contingencies = contingencies
        self.implicit_contingencies = []  # type: List[con.Contingency]
        self.states = None

        self._generate_state_list()
        self._generate_implicit_reactions_contingencies()
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

    def _generate_implicit_reactions_contingencies(self):
        # generate reaction-contingency pairs implied by the explicitly passed reactions and contingencies.
        # e.g. if A_ppi_B requires A-{P}, and there exists a dephosphorylation reaction for A,
        # we add a negative A_ppi_B that has as its contingency NOT A-{P}.
        pass

    def _assert_consistency(self):
        pass
