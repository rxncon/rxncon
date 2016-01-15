from typing import List

import rxncon.core.contingency as con
import rxncon.core.reaction as rxn


class RxnConSystem:
    def __init__(self, contingencies: List[con.Contingency]):
        self.implicit_reactions = []  # type: List[rxn.Reaction]
        self.contingencies = contingencies
        self.implicit_contingencies = []  # type: List[con.Contingency]

        self._generate_implicit_reactions_contingencies()
        self._assert_consistency()

    @property
    def reactions(self):
        rxns = []
        for contingency in self.contingencies:
            if contingency.target not in rxns:
                rxns.append(contingency.target)

        return rxns

    @property
    def source_states(self):
        states = []
        for reaction in self.reactions:
            if reaction.source and reaction.source not in states:
                states.append(reaction.source)

        return states

    @property
    def product_states(self):
        states = []
        for reaction in self.reactions:
            if reaction.product and reaction.product not in states:
                states.append(reaction.source)

        return states

    @property
    def effector_states(self):
        states = []
        for contingency in self.contingencies:
            effector_states = contingency.effector.states

            for state in effector_states:
                if state not in states:
                    states.append(state)

        return states

    @property
    def states(self):
        return self.source_states + self.product_states + self.effector_states

    def _generate_implicit_reactions_contingencies(self):
        # generate reaction-contingency pairs implied by the explicitly passed reactions and contingencies.
        # e.g. if A_ppi_B requires A-{P}, and there exists a dephosphorylation reaction for A,
        # we add a negative A_ppi_B that has as its contingency NOT A-{P}.
        pass

    def _assert_consistency(self):
        pass
