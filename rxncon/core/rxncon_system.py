from typing import List
from typecheck import typecheck

from rxncon.core.contingency import ContingencyType, Contingency
from rxncon.core.reaction import Reaction
from rxncon.core.state import State


class RxnConSystem:
    @typecheck
    def __init__(self, reactions: List[Reaction], contingencies: List[Contingency]):
        self.reactions = reactions
        self.contingencies = contingencies

        self._assert_consistency()

    @typecheck
    def quantitative_contingencies_for_reaction(self, reaction: Reaction) -> List[Contingency]:
        return [x for x in self.contingencies if x.target == reaction and x.type
                in [ContingencyType.positive, ContingencyType.negative]]

    @typecheck
    def strict_contingencies_for_reaction(self, reaction: Reaction) -> List[Contingency]:
        return [x for x in self.contingencies if x.target == reaction and x.type
                in [ContingencyType.requirement, ContingencyType.inhibition]]

    @property
    @typecheck
    def product_states(self) -> List[State]:
        states = []
        for reaction in self.reactions:
            states += reaction.products

        return sorted(list(set(states)))

    @property
    @typecheck
    def source_states(self) -> List[State]:
        states = []
        for reaction in self.reactions:
            states += reaction.sources

        return sorted(list(set(states)))

    def _assert_consistency(self):
        required_states = []
        for contingency in self.contingencies:
            required_states += contingency.effector.states

        reacting_states = self.product_states + self.source_states

        for state in required_states:
            assert state in reacting_states, \
                'State {0} appears in contingencies, but is neither produced or consumed'.format(str(state))
