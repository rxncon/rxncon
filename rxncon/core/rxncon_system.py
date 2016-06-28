from typing import List
import typecheck as tc

from rxncon.core.contingency import ContingencyType, Contingency
from rxncon.core.reaction import Reaction


class RxnConSystem:
    @tc.typecheck
    def __init__(self, reactions: List[Reaction], contingencies: List[Contingency]):
        self.reactions = reactions
        self.contingencies = contingencies

        self._assert_consistency()

    def quantitative_contingencies_for_reaction(self, reaction: Reaction) -> List[Contingency]:
        return [x for x in self.contingencies if x.target == reaction and x.type
                in [ContingencyType.positive, ContingencyType.negative]]

    def strict_contingencies_for_reaction(self, reaction: Reaction) -> List[Contingency]:
        return [x for x in self.contingencies if x.target == reaction and x.type
                in [ContingencyType.requirement, ContingencyType.inhibition]]

    @property
    def product_states(self):
        states = []
        for reaction in self.reactions:
            states += reaction.products

        return list(set(states))

    @property
    def source_states(self):
        states = []
        for reaction in self.reactions:
            states += reaction.sources

        return list(set(states))

    def _assert_consistency(self):
        required_states = []
        for contingency in self.contingencies:
            required_states += contingency.effector.states

        reacting_states = self.product_states + self.source_states

        for state in required_states:
            assert state in reacting_states, \
                'State {0} appears in contingencies, but is neither produced or consumed'.format(str(state))
