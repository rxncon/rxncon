from typing import List
import typecheck as tc

import rxncon.core.contingency as con
import rxncon.core.reaction as rxn


class RxnConSystem:
    @tc.typecheck
    def __init__(self, reactions: List[rxn.Reaction], contingencies: List[con.Contingency]):
        self.reactions = reactions
        self.implicit_reactions = []  # type: List[rxn.Reaction]
        self.contingencies = contingencies
        self.implicit_contingencies = []  # type: List[con.Contingency]

        self._generate_implicit_reactions_contingencies()
        self._assert_consistency()

    @property
    def source_states(self):
        return list({reaction.source for reaction in self.reactions if reaction.source})

    @property
    def product_states(self):
        return list({reaction.product for reaction in self.reactions if reaction.product})

    @property
    def effector_states(self):
        return list({state for contingency in self.contingencies for state in contingency.effector.states})

    @property
    def states(self):
        return list({state for state in self.source_states + self.product_states + self.effector_states})

    def _generate_implicit_reactions_contingencies(self):
        # generate reaction-contingency pairs implied by the explicitly passed reactions and contingencies.
        # e.g. if A_ppi_B requires A-{P}, and there exists a dephosphorylation reaction for A,
        # we add a negative A_ppi_B that has as its contingency NOT A-{P}.
        pass

    def _assert_consistency(self):
        pass
