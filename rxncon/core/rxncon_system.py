from typing import List
import typecheck as tc

import rxncon.core.contingency as con
import rxncon.core.reaction as rxn
import rxncon.core.effector as eff


class RxnConSystem:
    @tc.typecheck
    def __init__(self, reactions: List[rxn.Reaction], contingencies: List[con.Contingency]):
        self.reactions = reactions
        self.implicit_reactions = []  # type: List[rxn.Reaction]
        self.contingencies = contingencies
        self.source_contingencies = []  # type: List[con.Contingency]

        self._generate_source_contingencies()
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

    def quantitative_contingencies_for_reaction(self, reaction: rxn.Reaction) -> List[con.Contingency]:
        return [x for x in self.contingencies if x.target == reaction and x.type in [con.ContingencyType.positive, con.ContingencyType.negative]]

    def strict_contingencies_for_reaction(self, reaction: rxn.Reaction) -> List[con.Contingency]:
        return [x for x in self.contingencies if x.target == reaction and x.type in [con.ContingencyType.requirement, con.ContingencyType.inhibition]]

    def source_contingencies_for_reaction(self, reaction: rxn.Reaction) -> List[con.Contingency]:
        return [x for x in self.source_contingencies if x.target == reaction]

    def _generate_source_contingencies(self):
        # generate reaction-contingency pairs implied by the explicitly passed reactions and contingencies.
        # e.g. if A_ppi_B requires A-{P}, and there exists a dephosphorylation reaction for A,
        # we add a negative A_ppi_B that has as its contingency NOT A-{P}.
        for reaction in self.reactions:
            # @todo bidirectionality???
            if reaction.influence == rxn.Influence.positive or reaction.influence == rxn.Influence.bidirectional:
                self.source_contingencies.append(con.Contingency(reaction,
                                                                 con.ContingencyType.inhibition,
                                                                 eff.StateEffector(reaction.product)))
            elif reaction.influence == rxn.Influence.negative:
                self.source_contingencies.append(con.Contingency(reaction,
                                                                 con.ContingencyType.requirement,
                                                                 eff.StateEffector(reaction.source)))

            elif reaction.influence == rxn.Influence.transfer:
                self.source_contingencies.append(con.Contingency(reaction,
                                                                 con.ContingencyType.requirement,
                                                                 eff.StateEffector(reaction.source)))

                self.source_contingencies.append(con.Contingency(reaction,
                                                                 con.ContingencyType.inhibition,
                                                                 eff.StateEffector(reaction.product)))


    def _assert_consistency(self):
        pass
