from typing import List, FrozenSet
from typecheck import typecheck
from collections import defaultdict

from rxncon.core.contingency import ContingencyType, Contingency
from rxncon.core.reaction import Reaction, ReactionTerm, MoleculeReactionTerm
from rxncon.core.state import State
from rxncon.core.spec import MolSpec, BondSpec


class RxnConSystem:
    @typecheck
    def __init__(self, reactions: List[Reaction], contingencies: List[Contingency]):
        self.reactions = reactions
        self.contingencies = contingencies

        self._expand_fully_neutral_states()
        self._assert_consistency()

    @typecheck
    def contingencies_for_reaction(self, reaction: Reaction) -> List[Contingency]:
        return [x for x in self.contingencies if x.target == reaction]

    @typecheck
    def q_contingencies_for_reaction(self, reaction: Reaction) -> List[Contingency]:
        return [x for x in self.contingencies if x.target == reaction and x.type
                in [ContingencyType.positive, ContingencyType.negative]]

    @typecheck
    def s_contingencies_for_reaction(self, reaction: Reaction) -> List[Contingency]:
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
    def synthesised_states(self) -> List[State]:
        states = []
        for reaction in self.reactions:
            states += reaction.synthesised_states

        return sorted(list(set(states)))

    @property
    @typecheck
    def states(self) -> List[State]:
        return sorted(list(set(self.produced_states + self.consumed_states + self.synthesised_states)))

    @typecheck
    def states_for_component(self, component: MolSpec) -> List[State]:
        assert component.is_component_spec
        return [x for x in self.states if component in x.components]

    @typecheck
    def states_for_component_grouped(self, component: MolSpec) -> List[List[State]]:
        states = self.states_for_component(component)
        grouped = defaultdict(list)

        while states:
            state = states.pop()
            if isinstance(state.target, MolSpec):
                grouped[(state.target, state.definition)].append(state)
            elif isinstance(state.target, BondSpec):
                if state.target.first.to_component_spec() == component:
                    grouped[(state.target.first, state.definition)].append(state)
                if state.target.second.to_component_spec() == component:
                    grouped[(state.target.second, state.definition)].append(state)
            else:
                raise Exception('State target is neither MolSpec nor BondSpec')

        return list(grouped.values())

    def _expand_fully_neutral_states(self):
        for reaction in self.reactions:
            self._expand_reaction_terms(reaction.terms_lhs)
            self._expand_reaction_terms(reaction.terms_rhs)

    def _expand_reaction_terms(self, terms: List[ReactionTerm]):
        for term in terms:
            if isinstance(term, MoleculeReactionTerm) and term.is_fully_neutral:
                term.states = [x for x in self.states_for_component(term.spec) if x.is_neutral]

    def _assert_consistency(self):
        required_states = []
        for contingency in self.contingencies:
            required_states += contingency.effector.states

        for state in required_states:
            assert state in self.states, \
                'State {0} appears in contingencies, but is neither produced nor consumed'.format(str(state))


