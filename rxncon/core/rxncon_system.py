from typing import List, Dict, Tuple
from typecheck import typecheck
from collections import defaultdict

from rxncon.core.contingency import ContingencyType, Contingency
from rxncon.core.reaction import Reaction, ReactionTerm, MoleculeReactionTerm
from rxncon.core.state import State, StateDef
from rxncon.core.spec import MolSpec, BondSpec


class RxnConSystem:
    @typecheck
    def __init__(self, reactions: List[Reaction], contingencies: List[Contingency]):
        self.reactions = reactions
        self.contingencies = contingencies

        self._components         = []
        self._states             = []
        self._produced_states    = []
        self._consumed_states    = []
        self._synthesised_states = []

        self._expand_fully_neutral_states()
        self._assert_consistency()

    @typecheck
    def components(self) -> List[MolSpec]:
        if not self._components:
            self._calculate_components()
        return self._components

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
    def states(self) -> List[State]:
        if not self._states:
            self._calculate_states()
        return self._states

    @property
    @typecheck
    def produced_states(self) -> List[State]:
        if not self._produced_states:
            self._calculate_produced_states()
        return self._produced_states

    @property
    @typecheck
    def consumed_states(self) -> List[State]:
        if not self._consumed_states:
            self._calculate_consumed_states()
        return self._consumed_states

    @property
    @typecheck
    def synthesised_states(self) -> List[State]:
        if not self._synthesised_states:
            self._calculate_synthesised_states()
        return self._synthesised_states

    @typecheck
    def states_for_component(self, component: MolSpec) -> List[State]:
        assert component.is_component_spec
        return [x for x in self.states if component.to_non_struct_spec() in x.components]

    @typecheck
    def states_for_component_grouped(self, component: MolSpec) -> Dict[Tuple[MolSpec, StateDef], List[State]]:
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

        return grouped

    def _calculate_components(self):
        components = []
        for reaction in self.reactions:
            components += [spec.to_non_struct_spec() for spec in reaction.components_lhs] + \
                          [spec.to_non_struct_spec() for spec in reaction.components_rhs]

        self._components = list(set(components))

    def _calculate_states(self):
        self._states = list(set(self.produced_states + self.consumed_states + self.synthesised_states))

    def _calculate_produced_states(self):
        states = []
        for reaction in self.reactions:
            states += [state.to_non_struct_state() for state in reaction.produced_states]

        self._produced_states = list(set(states))

    def _calculate_consumed_states(self):
        states = []
        for reaction in self.reactions:
            states += [state.to_non_struct_state() for state in reaction.consumed_states]

        self._consumed_states = list(set(states))

    def _calculate_synthesised_states(self):
        states = []
        for reaction in self.reactions:
            states += [state.to_non_struct_state() for state in reaction.synthesised_states]

        self._synthesised_states = list(set(states))

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
            required_states += [state.to_non_struct_state() for state in contingency.effector.states]

        for state in required_states:
            assert state in self.states, \
                'State {0} appears in contingencies, but is neither produced nor consumed'.format(str(state))
