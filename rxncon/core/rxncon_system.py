from typing import List
import typecheck as tc

import rxncon.core.contingency as con
import rxncon.core.reaction as rxn
import rxncon.core.component as com


class RxnConSystem:
    @tc.typecheck
    def __init__(self, reactions: List[rxn.Reaction], contingencies: List[con.Contingency]):
        self.reactions = reactions
        self.implicit_reactions = []  # type: List[rxn.Reaction]
        self.contingencies = contingencies
        self.source_contingencies = []  # type: List[con.Contingency]

        self._assert_consistency()

    def quantitative_contingencies_for_reaction(self, reaction: rxn.Reaction) -> List[con.Contingency]:
        return [x for x in self.contingencies if x.target == reaction and x.type
                in [con.ContingencyType.positive, con.ContingencyType.negative]]

    def strict_contingencies_for_reaction(self, reaction: rxn.Reaction) -> List[con.Contingency]:
        return [x for x in self.contingencies if x.target == reaction and x.type
                in [con.ContingencyType.requirement, con.ContingencyType.inhibition]]

    def components_for_reaction(self, reaction: rxn.Reaction) -> List[com.Component]:
        components = reaction.components

        states = []

        for contingency in self.strict_contingencies_for_reaction(reaction) + self.quantitative_contingencies_for_reaction(reaction):
            states += contingency.effector.states

        for state in states:
            components += state.components

        return components


    def _assert_consistency(self):
        # @todo check that every contingency is produced.
        pass
