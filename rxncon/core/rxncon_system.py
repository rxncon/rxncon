"""Module containing the class RxnconSystem, which contains all information (reactions, contingencies)
pertaining to a rxncon system. Some information is inherently non-local, e.g. the neutral states that
get created by a synthesis reaction; this information is determined here. Also some validation is performed,
e.g. a check that all States appearing in contingencies are actually appearing in reactions."""


from typing import List, Dict, Tuple
from itertools import product
from collections import defaultdict, Counter
import logging

from rxncon.core.contingency import ContingencyType, Contingency
from rxncon.core.effector import Effector, AndEffector, OrEffector, NotEffector, StateEffector
from rxncon.core.reaction import Reaction, ReactionTerm
from rxncon.core.state import State, FullyNeutralState
from rxncon.core.spec import Spec
from rxncon.venntastic.sets import UniversalSet, Intersection, Set  # pylint: disable=unused-import

LOGGER = logging.getLogger(__name__)


class RxnConSystem:  # pylint: disable=too-many-instance-attributes
    """RxnConSystem holds all reactions and contingencies pertaining to a rxncon system."""
    def __init__(self, reactions: List[Reaction], contingencies: List[Contingency]) -> None:
        self.reactions = reactions
        self.contingencies = contingencies

        self._components = []  # type: List[Spec]
        self._states = []  # type: List[State]
        self._produced_states = []  # type: List[State]
        self._consumed_states = []  # type: List[State]
        self._synthesised_states = []  # type: List[State]
        self._global_states = []  # type: List[State]

        # Components get synthesised in a 'fully neutral' state. What this state is depends on all reactions
        # in the system, so this 'fully neutral' state has to be expanded into the combination of all neutral
        # states pertaining to the component.
        self._expand_fully_neutral_states()
        # These states are lazily generated.
        self._calculate_produced_states()
        self._calculate_consumed_states()
        self._calculate_synthesised_states()
        self._calculate_global_states()
        self._calculate_states()

        # Contingencies containing states at a non-elemental resolution are expanded to an OR of elemental
        # contingencies.
        self._expand_non_elemental_states()
        self._structure_contingencies()

        self.validate()

    def reaction_number(self, reaction: Reaction) -> int:
        return self.reactions.index(reaction)

    def validate(self) -> None:
        if not self.reactions:
            raise AssertionError('No reactions, boring!')

        missing_states = self._missing_states()
        if missing_states:
            raise AssertionError('State(s) {0} appear(s) in contingencies, but is not produced or consumed'
                                 .format(', '.join(str(state) for state in missing_states)))

        missing_reactions = self._missing_reactions()
        if missing_reactions:
            raise AssertionError('Reactions(s) {0} appear(s) in contingencies, but are not defined in reaction list'
                                 .format(', '.join(str(reaction) for reaction in missing_reactions)))

        unsatisfiable_contingencies = self._unsatisfiable_contingencies()
        if unsatisfiable_contingencies:
            reason_str = '\n'.join('{} : {}'.format(rxn, reason) for rxn, reason in unsatisfiable_contingencies)
            raise AssertionError('Unsatisfiable reaction contingencies:\n{}'.format(reason_str))

    def components(self) -> List[Spec]:
        if not self._components:
            self._calculate_components()
        return self._components

    def contingencies_for_reaction(self, reaction: Reaction) -> List[Contingency]:
        assert reaction in self.reactions
        return [x for x in self.contingencies if x.reaction == reaction]

    def q_contingencies_for_reaction(self, reaction: Reaction) -> List[Contingency]:
        assert reaction in self.reactions
        return [x for x in self.contingencies if x.reaction == reaction and x.contingency_type
                in [ContingencyType.positive, ContingencyType.negative]]

    def s_contingencies_for_reaction(self, reaction: Reaction) -> List[Contingency]:
        assert reaction in self.reactions
        return [x for x in self.contingencies if x.reaction == reaction and x.contingency_type
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
    def global_states(self) -> List[State]:
        if not self._global_states:
            self._calculate_global_states()
        return self._global_states

    def states_for_component(self, component: Spec) -> List[State]:
        """Returns all the States that live on a certain Component."""
        assert component.is_component_spec
        return [x for x in self.states if component.to_non_struct_spec() in x.components]

    def states_for_component_grouped(self, component: Spec) -> Dict[Spec, List[State]]:
        """Returns all the States that live on a certain component, grouped by their Spec.
        Within a group, each State is mutually exclusive with each other State."""
        states = self.states_for_component(component)  # type: List[State]
        grouped = defaultdict(list)  # type: Dict[Spec, List[State]]

        while states:
            state = states.pop()

            for spec in state.specs:
                if spec.to_component_spec() != component:
                    continue

                grouped[spec].append(state)

        return grouped

    def complement_states(self, state: State) -> List[State]:
        """Returns all States mutually exclusive with the State given. For multi-component States,
        all components are included."""
        states = []  # type: List[State]
        for component in state.components:
            states += self.complement_states_for_component(component, state)

        return states

    def complement_states_for_component(self, component: Spec, state: State) -> List[State]:
        """Returns all States mutually exclusive with the State given that live on the Component given."""
        if not state.is_structured:
            for group in self.states_for_component_grouped(component).values():
                if state in group:
                    return [x for x in group if x != state]
            raise AssertionError
        else:
            # The structure information is first thrown away to determine the mutually exclusive
            # states, and then applied again (as far as possible).
            complements = self.complement_states_for_component(component.to_non_struct_spec(),
                                                               state.to_non_structured())
            return [x.to_structured_from_state(state) for x in complements]

    def _calculate_components(self) -> None:
        """Determines all components in the system, and stores it in self._components."""
        components = []  # type: List[Spec]
        for reaction in self.reactions:
            components += [spec.to_non_struct_spec() for spec in reaction.components_lhs] + \
                          [spec.to_non_struct_spec() for spec in reaction.components_rhs]

        self._components = list(set(components))

    def _calculate_states(self) -> None:
        self._states = list(
            set(self.produced_states + self.consumed_states + self.synthesised_states + self.global_states))

    def _calculate_produced_states(self) -> None:
        states = []  # type: List[State]
        for reaction in self.reactions:
            states += [state.to_non_structured() for state in reaction.produced_states]

        self._produced_states = list(set(states))

    def _calculate_consumed_states(self) -> None:
        states = []  # type: List[State]
        for reaction in self.reactions:
            states += [state.to_non_structured() for state in reaction.consumed_states]

        self._consumed_states = list(set(states))

    def _calculate_synthesised_states(self) -> None:
        states = []  # type: List[State]
        for reaction in self.reactions:
            states += [state.to_non_structured() for state in reaction.synthesised_states]

        self._synthesised_states = list(set(states))

    def _calculate_global_states(self) -> None:
        states = []  # type: List[State]
        for contingency in self.contingencies:
            states += [state for state in contingency.effector.states if state.is_global]

        self._global_states = list(set(states))

    def _expand_fully_neutral_states(self) -> None:
        for reaction in self.reactions:
            self._expand_reaction_terms(reaction.terms_lhs)
            self._expand_reaction_terms(reaction.terms_rhs)
            reaction.invalidate_state_cache()

    def _expand_reaction_terms(self, terms: List[ReactionTerm]) -> None:
        """Expands the reaction rules that contain a FullyNeutralState, which appear in synthesis reactions.
        This state is expanded in a combination of all neutral states for the component."""
        for term in terms:
            if FullyNeutralState() in term.states:
                existing_states = [state for state in term.states if state != FullyNeutralState()]
                new_states = [state for component in term.specs for state in self.states_for_component(component)
                              if state.is_neutral and state != FullyNeutralState() and state not in existing_states and
                              not any(state.is_mutually_exclusive_with(existing) for existing in existing_states)]

                term.states = existing_states + new_states

    def _expand_non_elemental_states(self) -> None:
        """Expands the non-elemental States appearing in contingencies, these get expanded into a Union of
        the elemental States of which the non-elemental State is a superset."""
        def expanded_effector(effector: Effector) -> Effector:
            if isinstance(effector, StateEffector):
                if effector.expr.is_elemental:
                    return effector
                else:
                    elemental_states = [state.to_structured_from_state(effector.expr)
                                        for state in self.states if state.is_subset_of(effector.expr)]
                    assert elemental_states, 'Could not find elemental states which are subset of the non-elemental ' \
                                             'state {}'.format(effector.expr)
                    assert all(state.is_elemental for state in elemental_states)

                    LOGGER.info('expanded_effector: {} -> {}'
                                .format(str(effector.expr), ' | '.join(str(x) for x in elemental_states)))
                    return OrEffector(*(StateEffector(x) for x in elemental_states), name=str(effector.expr))
            elif isinstance(effector, AndEffector):
                return AndEffector(*(expanded_effector(x) for x in effector.exprs), name=effector.name,
                                   equivs=effector.equivs)
            elif isinstance(effector, OrEffector):
                return OrEffector(*(expanded_effector(x) for x in effector.exprs), name=effector.name,
                                  equivs=effector.equivs)
            elif isinstance(effector, NotEffector):
                return NotEffector(expanded_effector(effector.expr), name=effector.name)
            else:
                raise AssertionError

        for contingency in self.contingencies:
            contingency.effector = expanded_effector(contingency.effector)

    def _structure_contingencies(self) -> None:
        """Structures (i.e. augment with topological data) all contingencies."""
        structured_conts = []

        for rxn in self.reactions:
            counter_start = 2
            for con in self.contingencies_for_reaction(rxn):
                structured_conts.append(con.to_structured(counter_start))
                counter_start += 100

        self.contingencies = structured_conts

    def _missing_states(self) -> List[State]:
        """Returns the States that appear in contingencies, but which are not produced, consumed, synthesised
        or degraded by the Reactions."""
        required_states = []  # type: List[State]
        for contingency in self.contingencies:
            required_states += [state.to_non_structured() for state in contingency.effector.states if
                                not state.is_global]

        return [state for state in required_states if state not in self.states]

    def _missing_reactions(self) -> List[Reaction]:
        """Returns the Reactions that appear in the contingency list, but which are not in the Reaction list."""
        required_reactions = []
        for contingency in self.contingencies:
            required_reactions.append(contingency.reaction)

        return [reaction for reaction in required_reactions if reaction not in self.reactions]

    def _unsatisfiable_contingencies(self) -> List[Tuple[Reaction, str]]:
        """Determines the contingencies that are not satisfiable, returns a list of (Reaction, str) where
        the str contains the human-readable reason for the contingency not being satisfiable."""
        unsatisfiable = []

        for reaction in self.reactions:
            contingencies = self.s_contingencies_for_reaction(reaction)

            # Make sure the contingency does not contain the states produced / consumed by the reaction.
            states = (state for contingency in contingencies
                      for state in contingency.effector.states)

            for state in states:
                # We need to be talking about the states mentioning the reactants.
                if not all(spec.struct_index in (0, 1) for spec in state.specs):
                    continue

                if any(count > 1 for count in Counter(reaction.components_rhs).values()):
                    continue

                # States appear non-structured in the '.produced_states' etc. properties.
                state = state.to_non_structured()

                if state in reaction.produced_states and state not in reaction.synthesised_states:
                    unsatisfiable.append((reaction, 'Produced state {} appears in contingencies'.format(str(state))))
                if state in reaction.consumed_states and state not in reaction.degraded_states:
                    unsatisfiable.append((reaction, 'Consumed state {} appears in contingencies'.format(str(state))))

            # Make sure at least one solution is there (this might still contain mutually exclusive states)
            total_set = UniversalSet()  # type: Set[State]
            for contingency in contingencies:
                total_set = Intersection(total_set,
                                         contingency.to_venn_set())  # pylint: disable=redefined-variable-type

            solutions = total_set.calc_solutions()
            if len(solutions) == 0:
                unsatisfiable.append((reaction, 'Zero consistent solutions found.'))

            # Make sure that at least one solution doesn't contain mutually exclusive states.
            local_unsatisfiables = []
            at_least_one_consistent_soln = False

            for solution in solutions:
                trues = [state for state, val in solution.items() if val]
                if any(state.is_mutually_exclusive_with(other) for state, other in product(trues, trues)):
                    state, other = next((state, other) for state, other in product(trues, trues) if
                                        state.is_mutually_exclusive_with(other))
                    local_unsatisfiables.append(
                        (reaction, 'State {} mutually exclusive with {}.'.format(str(state), str(other))))
                else:
                    at_least_one_consistent_soln = True

            if not at_least_one_consistent_soln:
                unsatisfiable += local_unsatisfiables

        return unsatisfiable
