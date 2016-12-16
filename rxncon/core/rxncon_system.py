from typing import List, Dict, Tuple
from collections import defaultdict

from rxncon.core.contingency import ContingencyType, Contingency
from rxncon.core.effector import Effector, AndEffector, OrEffector, NotEffector, StateEffector
from rxncon.core.reaction import Reaction, ReactionTerm
from rxncon.core.state import State, StateDef, FullyNeutralState
from rxncon.core.spec import Spec


class RxnConSystem:
    def __init__(self, reactions: List[Reaction], contingencies: List[Contingency]):
        self.reactions     = reactions
        self.contingencies = contingencies

        self._components         = []  # type: List[Spec]
        self._states             = []  # type: List[State]
        self._produced_states    = []  # type: List[State]
        self._consumed_states    = []  # type: List[State]
        self._synthesised_states = []  # type: List[State]
        self._global_states      = []  # type: List[State]

        self._expand_fully_neutral_states()
        self._calculate_produced_states()
        self._calculate_consumed_states()
        self._calculate_synthesised_states()
        self._calculate_global_states()
        self._calculate_states()

        self._expand_non_elemental_contingencies()
        self._structure_contingencies()

        self.validate()

    def validate(self):
        if not self.reactions:
            raise AssertionError('No reactions, boring!')

        missing_states = self._missing_states()

        if missing_states:
            raise AssertionError('State(s) {0} appear(s) in contingencies, but is not produced or consumed'
                                 .format(', '.join(str(state) for state in missing_states)))

    def components(self) -> List[Spec]:
        if not self._components:
            self._calculate_components()
        return self._components

    def contingencies_for_reaction(self, reaction: Reaction) -> List[Contingency]:
        return [x for x in self.contingencies if x.target == reaction]

    def q_contingencies_for_reaction(self, reaction: Reaction) -> List[Contingency]:
        return [x for x in self.contingencies if x.target == reaction and x.type
                in [ContingencyType.positive, ContingencyType.negative]]

    def s_contingencies_for_reaction(self, reaction: Reaction) -> List[Contingency]:
        return [x for x in self.contingencies if x.target == reaction and x.type
                in [ContingencyType.requirement, ContingencyType.inhibition]]

    @property
    def states(self) -> List[State]:
        if not self._states:
            self._calculate_states()
        return self._states

    @property
    def produced_states(self) -> List[State]:
        if not self._produced_states:
            self._calculate_produced_states()
        return self._produced_states

    @property
    def consumed_states(self) -> List[State]:
        if not self._consumed_states:
            self._calculate_consumed_states()
        return self._consumed_states

    @property
    def synthesised_states(self) -> List[State]:
        if not self._synthesised_states:
            self._calculate_synthesised_states()
        return self._synthesised_states

    @property
    def global_states(self):
        if not self._global_states:
            self._calculate_global_states()
        return self._global_states

    def states_for_component(self, component: Spec) -> List[State]:
        assert component.is_component_spec
        return [x for x in self.states if component.to_non_struct_spec() in x.components]

    def states_for_component_grouped(self, component: Spec) -> Dict[Spec, List[State]]:
        states = self.states_for_component(component)
        grouped = defaultdict(list)

        while states:
            state = states.pop()

            for spec in state.specs:
                if spec.to_component_spec() != component:
                    continue

                grouped[spec].append(state)

        return grouped

    def complementary_states(self, state: State) -> List[State]:
        states = []
        for component in state.components:
            states += self.complementary_states_for_component(component, state)

        return states

    def complementary_states_for_component(self, component: Spec, state: State) -> List[State]:
        for group in self.states_for_component_grouped(component.to_non_struct_spec()).values():
            if state.to_non_structured_state() in group:
                complements = [x for x in group if x != state.to_non_structured_state()]
                return complements

    def _calculate_components(self):
        components = []
        for reaction in self.reactions:
            components += [spec.to_non_struct_spec() for spec in reaction.components_lhs] + \
                          [spec.to_non_struct_spec() for spec in reaction.components_rhs]

        self._components = list(set(components))

    def _calculate_states(self):
        self._states = list(set(self.produced_states + self.consumed_states + self.synthesised_states + self.global_states))

    def _calculate_produced_states(self):
        states = []
        for reaction in self.reactions:
            states += [state.to_non_structured_state() for state in reaction.produced_states]

        self._produced_states = list(set(states))

    def _calculate_consumed_states(self):
        states = []
        for reaction in self.reactions:
            states += [state.to_non_structured_state() for state in reaction.consumed_states]

        self._consumed_states = list(set(states))

    def _calculate_synthesised_states(self):
        states = []
        for reaction in self.reactions:
            states += [state.to_non_structured_state() for state in reaction.synthesised_states]

        self._synthesised_states = list(set(states))

    def _calculate_global_states(self):
        states = []
        for contingency in self.contingencies:
            states += [state for state in contingency.effector.states if state.is_global]

        self._global_states = list(set(states))

    def _expand_fully_neutral_states(self):
        for reaction in self.reactions:
            self._expand_reaction_terms(reaction.terms_lhs)
            self._expand_reaction_terms(reaction.terms_rhs)
            reaction.invalidate_state_cache()

    def _expand_reaction_terms(self, terms: List[ReactionTerm]):
        for term in terms:
            if term.is_fully_neutral:
                term.states = [state for component in term.specs for state in self.states_for_component(component)
                               if state.is_neutral and state != FullyNeutralState()]

    def _expand_non_elemental_contingencies(self):
        def expanded_effector(effector: Effector):
            def MultiOrEffector(state_list: List[State]):
                if len(state_list) == 1:
                    return StateEffector(state_list[0])
                else:
                    return OrEffector(MultiOrEffector(state_list[1:]), StateEffector(state_list[0]))

            if isinstance(effector, StateEffector):
                if effector.expr.is_elemental:
                    return effector
                else:
                    elemental_states = [state for state in self.states if state.is_subset_of(effector.expr)]
                    assert elemental_states, 'Could not find elemental states which are subset of the non-elemental ' \
                                             'state {}'.format(effector.expr)
                    assert all(state.is_elemental for state in elemental_states)
                    multi_or_effector = MultiOrEffector(elemental_states)
                    multi_or_effector.name = effector.name
                    return multi_or_effector
            elif isinstance(effector, AndEffector):
                and_effector = AndEffector(*(expanded_effector(x) for x in effector.exprs))
                and_effector.name = effector.name
                return and_effector
            elif isinstance(effector, OrEffector):
                or_effector = OrEffector(*(expanded_effector(x) for x in effector.exprs))
                or_effector.name = effector.name
                return or_effector
            elif isinstance(effector, NotEffector):
                not_effector = NotEffector(expanded_effector(effector.expr))
                not_effector.name = effector.name
                return not_effector
            else:
                raise AssertionError

        for contingency in self.contingencies:
            contingency.effector = expanded_effector(contingency.effector)

    def _structure_contingencies(self):
        self.contingencies = [c.to_structured() for c in self.contingencies]

    def _missing_states(self):
        required_states = []
        for contingency in self.contingencies:
            required_states += [state.to_non_structured_state() for state in contingency.effector.states if not state.is_global]

        return [state for state in required_states if state not in self.states]
