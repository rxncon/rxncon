from typing import List
from typecheck import typecheck

from rxncon.core.contingency import ContingencyType, Contingency
from rxncon.core.reaction import Reaction, Reactant
from rxncon.core.state import State
from rxncon.core.spec import MolSpec


class RxnConSystem:
    @typecheck
    def __init__(self, reactions: List[Reaction], contingencies: List[Contingency]):
        self.reactions = reactions
        self.contingencies = contingencies

        self._expand_fully_neutral_states()
        self._assert_consistency()

    @typecheck
    def quant_contingencies(self, reaction: Reaction) -> List[Contingency]:
        return [x for x in self.contingencies if x.target == reaction and x.type
                in [ContingencyType.positive, ContingencyType.negative]]

    @typecheck
    def strict_contingencies(self, reaction: Reaction) -> List[Contingency]:
        return [x for x in self.contingencies if x.target == reaction and x.type
                in [ContingencyType.requirement, ContingencyType.inhibition]]

    @property
    @typecheck
    def produced_states(self) -> List[State]:
        states = []
        for reaction in self.reactions:
            states += reaction.produced_states

        return sorted(list(set(states)))

    @property
    @typecheck
    def consumed_states(self) -> List[State]:
        states = []
        for reaction in self.reactions:
            states += reaction.consumed_states

        return sorted(list(set(states)))

    @property
    @typecheck
    def states(self) -> List[State]:
        return sorted(list(set(self.produced_states + self.consumed_states)))

    @typecheck
    def states_for_component(self, component: MolSpec) -> List[State]:
        assert component.is_component_spec
        return [x for x in self.states if component in x.components]

    def _expand_fully_neutral_states(self):
        for reaction in self.reactions:
            self._expand(reaction.reactants_lhs)
            self._expand(reaction.reactants_rhs)

    def _expand(self, reactants: List[Reactant]):
        for reactant in reactants:
            if reactant.is_molecule_reactant and reactant.is_fully_neutral:
                reactant.value = [x for x in self.states_for_component(reactant.spec) if x.is_neutral]

    def _assert_consistency(self):
        required_states = []
        for contingency in self.contingencies:
            required_states += contingency.effector.states

        for state in required_states:
            assert state in self.states, \
                'State {0} appears in contingencies, but is neither produced nor consumed'.format(str(state))
