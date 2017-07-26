from abc import ABCMeta
from copy import deepcopy
from enum import Enum
from itertools import product
from typing import List, Dict, Tuple, Optional

from rxncon.core.reaction import Reaction, OutputReaction
from rxncon.core.rxncon_system import RxnConSystem
from rxncon.core.spec import Spec
from rxncon.core.state import State, InteractionState
from rxncon.venntastic.sets import Set as VennSet, ValueSet, Intersection, Union, Complement, UniversalSet, EmptySet

MAX_STEADY_STATE_ITERS = 20


class BooleanModel:
    """Holds all data describing a Boolean model: a list of targets, a list of update rules and
    a list of initial conditions."""
    def __init__(self, targets: List['Target'], update_rules: List['UpdateRule'],
                 initial_conditions: 'BooleanModelState') -> None:
        self.update_rules = sorted(update_rules)
        self.initial_conditions = initial_conditions
        self._state_targets = {str(x): x for x in targets if isinstance(x, StateTarget)}
        self._reaction_targets = {str(x): x for x in targets if isinstance(x, ReactionTarget)}
        self._knockout_targets = {str(x): x for x in targets if isinstance(x, KnockoutTarget)}
        self._overexpression_targets = {str(x): x for x in targets if isinstance(x, OverexpressionTarget)}
        self._validate_update_rules()
        self._validate_initial_conditions()

        self.current_state = None  # type: Optional[BooleanModelState]

    def set_initial_condition(self, target: 'Target', value: bool) -> None:
        self.initial_conditions.set_target(target, value)

    def update_rule_by_target(self, target: 'Target') -> 'UpdateRule':
        for rule in self.update_rules:
            if rule.target == target:
                return rule

        raise KeyError

    def state_target_by_name(self, name: str) -> 'StateTarget':
        return self._state_targets[name]

    def reaction_target_by_name(self, name: str) -> 'ReactionTarget':
        return self._reaction_targets[name]

    def knockout_target_by_name(self, name: str) -> 'KnockoutTarget':
        return self._knockout_targets[name]

    def overexpression_target_by_name(self, name: str) -> 'OverexpressionTarget':
        return self._overexpression_targets[name]

    def step(self) -> None:
        """Takes one timestep in the Boolean model. This is rather inefficient, but not meant for
        actual simulations, this is only used in the unit tests that test all different motifs
        and their desired steady states."""
        if not self.current_state:
            self.current_state = deepcopy(self.initial_conditions)
        else:
            new_state = dict()
            for rule in self.update_rules:
                new_state[rule.target] = rule.factor.eval_boolean_func(self.current_state.target_to_value)

            self.current_state = BooleanModelState(new_state)

    def calc_steady_state(self) -> 'BooleanModelState':
        """Calculates the steady state by taking max MAX_STEADY_STATE_ITERS steps. If no steady state
        found, raises."""
        iters = 0

        while iters < MAX_STEADY_STATE_ITERS:
            prev = deepcopy(self.current_state)
            self.step()
            if prev == self.current_state:
                assert isinstance(prev, BooleanModelState)
                return prev

            iters += 1

        raise AssertionError('Could not find steady state.')

    def _validate_update_rules(self) -> None:
        """Assert that all targets appearing on the RHS in an update rule have their own LHS."""
        all_lhs_targets = []  # type: List[Target]
        all_rhs_targets = []  # type: List[Target]
        for rule in self.update_rules:
            all_lhs_targets.append(rule.target)
            all_rhs_targets += rule.factor_targets

        assert all(x in all_lhs_targets for x in all_rhs_targets)

    def _validate_initial_conditions(self) -> None:
        self.initial_conditions.validate_by_model(self)


class BooleanModelState:
    def __init__(self, target_to_value: Dict['Target', bool]) -> None:
        self.target_to_value = target_to_value

    def __eq__(self, other):
        if not isinstance(other, BooleanModelState):
            return NotImplemented
        else:
            return self.target_to_value == other.target_to_value

    def __getitem__(self, item):
        return self.target_to_value[item]

    def __str__(self):
        return str(self.target_to_value)

    def __repr__(self):
        return str(self)

    def set_target(self, target: 'Target', value: bool) -> None:
        self.target_to_value[target] = value

    def validate_by_model(self, model: BooleanModel) -> None:
        """Assert that all targets appearing in the model have a Boolean value assigned."""
        model_targets = [rule.target for rule in model.update_rules]
        config_targets = self.target_to_value.keys()

        assert set(model_targets) == set(config_targets) and len(model_targets) == len(config_targets)


class Target(metaclass=ABCMeta):
    """Abstract base class for the different targets."""
    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)


class ReactionTarget(Target):
    """Reaction target of the boolean model. For all non-degrading reactions the relation between
    rxncon reactions and Boolean targets is 1:1. The relation for degrading reactions is more difficult
    since (1) the contingencies determine what the reaction degrades (which obviously becomes problematic
    in the case of a logical disjunction), and (2) the degradation of bonds should produce empty binding
    partners. We refer to our paper."""
    def __init__(self, reaction_parent: Reaction, contingency_variant: Optional[int]=None,
                 interaction_variant: Optional[int] = None, contingency_factor: VennSet['StateTarget']=None) -> None:
        self.reaction_parent = reaction_parent  # type: Reaction
        self.produced_targets = [StateTarget(x) for x in reaction_parent.produced_states]  # type: List[StateTarget]
        self.consumed_targets = [StateTarget(x) for x in reaction_parent.consumed_states]  # type: List[StateTarget]
        self.synthesised_targets = [StateTarget(x) for x in
                                    reaction_parent.synthesised_states]  # type: List[StateTarget]
        self.degraded_targets = [StateTarget(x) for x in reaction_parent.degraded_states]  # type: List[StateTarget]

        self.contingency_variant_index = contingency_variant
        self.interaction_variant_index = interaction_variant

        if contingency_factor is None:
            self.contingency_factor = UniversalSet()  # type: VennSet[StateTarget]
        else:
            self.contingency_factor = contingency_factor  # type: VennSet[StateTarget]

    def __hash__(self) -> int:
        return hash(str(self))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Target):
            return NotImplemented
        return isinstance(other, ReactionTarget) and self.reaction_parent == other.reaction_parent and \
            self.contingency_variant_index == other.contingency_variant_index and \
            self.interaction_variant_index == other.interaction_variant_index

    def __str__(self) -> str:
        suffix = ''
        if self.interaction_variant_index is not None and self.contingency_variant_index is not None:
            suffix = '#c{}/i{}'.format(self.contingency_variant_index, self.interaction_variant_index)
        elif self.contingency_variant_index is not None and self.interaction_variant_index is None:
            suffix = '#c{}'.format(self.contingency_variant_index)
        elif self.interaction_variant_index is not None and self.contingency_variant_index is None:
            suffix = '#i{}'.format(self.interaction_variant_index)

        return str(self.reaction_parent) + suffix

    def __repr__(self) -> str:
        return str(self)

    def produces(self, state_target: 'StateTarget') -> bool:
        return state_target in self.produced_targets

    def consumes(self, state_target: 'StateTarget') -> bool:
        return state_target in self.consumed_targets

    def synthesises(self, state_target: 'StateTarget') -> bool:
        return state_target in self.synthesised_targets

    def degrades(self, state_target: 'StateTarget') -> bool:
        return state_target in self.degraded_targets

    @property
    def components_lhs(self) -> List[Spec]:
        return self.reaction_parent.components_lhs

    @property
    def components_rhs(self) -> List[Spec]:
        return self.reaction_parent.components_rhs

    @property
    def degraded_components(self) -> List[Spec]:
        return [component for component in self.components_lhs if component not in self.components_rhs]

    @property
    def synthesised_components(self) -> List[Spec]:
        return [component for component in self.components_rhs if component not in self.components_lhs]

    def degrades_component(self, spec: Spec) -> bool:
        assert spec.is_component_spec
        return spec in self.degraded_components

    def synthesises_component(self, spec: Spec) -> bool:
        assert spec.is_component_spec
        return spec in self.synthesised_components

    def is_output(self) -> bool:
        return isinstance(self.reaction_parent, OutputReaction)


class StateTarget(Target):
    """State target of the Boolean model. The relation between rxncon states and Boolean state targets
    is generally 1:1, but not quite: the components that carry no internal state are assigned so-called
    ComponentStateTargets (see next class.)"""
    def __init__(self, state_parent: State) -> None:
        self.state_parent = state_parent

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return str(self.state_parent)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Target):
            return NotImplemented
        return isinstance(other, StateTarget) and self.state_parent == other.state_parent

    def is_produced_by(self, reaction_target: ReactionTarget) -> bool:
        return reaction_target.produces(self)

    def is_consumed_by(self, reaction_target: ReactionTarget) -> bool:
        return reaction_target.consumes(self)

    def is_synthesised_by(self, reaction_target: ReactionTarget) -> bool:
        return reaction_target.synthesises(self)

    def is_degraded_by(self, reaction_target: ReactionTarget) -> bool:
        return reaction_target.degrades(self)

    def is_input(self) -> bool:
        return self.state_parent.is_global

    @property
    def components(self) -> List[Spec]:
        return self.state_parent.components

    def shares_component_with(self, other_target: 'StateTarget') -> bool:
        return any(x in other_target.components for x in self.components)

    @property
    def is_neutral(self) -> bool:
        return self.state_parent.is_neutral

    @property
    def is_homodimer(self) -> bool:
        return self.state_parent.is_homodimer

    @property
    def neutral_targets(self) -> List['StateTarget']:
        return [StateTarget(x) for x in self.state_parent.neutral_states]

    @property
    def is_interaction(self) -> bool:
        return isinstance(self.state_parent, InteractionState)

    def is_mutually_exclusive_with(self, other: 'StateTarget') -> bool:
        return self.state_parent.is_mutually_exclusive_with(other.state_parent)

    def complementary_state_targets(self, rxnconsys: RxnConSystem, component: Spec) -> List['StateTarget']:
        others = rxnconsys.complement_states_for_component(component, self.state_parent)
        return [StateTarget(x) for x in others]


class ComponentStateTarget(StateTarget):
    """ComponentStateTarget describes a rxncon component that carries no states."""
    def __init__(self, component: Spec) -> None:
        self.component = component

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Target):
            return NotImplemented
        return isinstance(other, ComponentStateTarget) and self.component == other.component

    def __str__(self) -> str:
        return str(self.component)

    def __repr__(self) -> str:
        return str(self)

    def __hash__(self) -> int:
        return hash(str(self))

    @property
    def components(self) -> List[Spec]:
        return [self.component]

    @property
    def is_neutral(self) -> bool:
        return True

    @property
    def is_interaction(self) -> bool:
        return False


class KnockoutTarget(ComponentStateTarget):
    """When enabled, KnockoutTargets are ANDed into the update rule for each target, one KnockoutTarget
    per component. The KnockoutTarget itself has a trivial update rule. By changing the initial conditions
    for the KnockoutTarget we can knock out all states for a particular component, making them always FALSE."""
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Target):
            return NotImplemented
        return isinstance(other, KnockoutTarget) and self.component == other.component

    def __str__(self) -> str:
        return 'Knockout<{}>'.format(str(self.component))

    def __repr__(self) -> str:
        return str(self)

    def __hash__(self) -> int:
        return hash(str(self))


class OverexpressionTarget(ComponentStateTarget):
    """Similar to KnockoutTarget, but now (a set of) states per component can be made always TRUE."""
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Target):
            return NotImplemented
        return isinstance(other, OverexpressionTarget) and self.component == other.component

    def __str__(self) -> str:
        return 'Overexpression<{}>'.format(str(self.component))

    def __repr__(self) -> str:
        return str(self)

    def __hash__(self) -> int:
        return hash(str(self))


class UpdateRule:
    """The UpdateRule for a Target gives the Boolean value of the target at time t+1 when we
    evaluate the factor at time t."""
    def __init__(self, target: Target, factor: VennSet[Target]) -> None:
        self.target = target
        self.factor = factor

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, UpdateRule):
            return NotImplemented
        else:
            return self.target == other.target and self.factor.is_equivalent_to(other.factor)

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, UpdateRule):
            return NotImplemented
        else:
            return str(self.target) < str(other.target)

    def __str__(self) -> str:
        return "target: {0}, factors: {1}".format(self.target, self.factor)

    @property
    def factor_targets(self) -> List[Target]:
        return self.factor.values


class SmoothingStrategy(Enum):
    """To overcome non biological oscillatory behaviour during the simulation we introduced
    smoothings."""
    no_smoothing = 'no_smoothing'
    smooth_production_sources = 'smooth_production_sources'


class KnockoutStrategy(Enum):
    """For which states should knockout rules be generated: none, the neutral ones, or all."""
    no_knockout = 'no_knockout'
    knockout_neutral_states = 'knockout_neutral_states'
    knockout_all_states = 'knockout_all_states'


class OverexpressionStrategy(Enum):
    """For which states should overexpression rules be generated: none, the neutral ones, or all."""
    no_overexpression = 'no_overexpression'
    overexpress_neutral_states = 'overexpress_neutral_states'
    overexpress_all_states = 'overexpress_all_states'


def boolean_model_from_rxncon(rxncon_sys: RxnConSystem,
                              smoothing_strategy: SmoothingStrategy=SmoothingStrategy.no_smoothing,
                              knockout_strategy: KnockoutStrategy=KnockoutStrategy.no_knockout,
                              overexpression_strategy: OverexpressionStrategy=OverexpressionStrategy.no_overexpression,
                              k_plus_strict: bool=True, k_minus_strict: bool=True) -> BooleanModel:

    def initial_conditions(reaction_targets: List[ReactionTarget], state_targets: List[StateTarget],
                           knockout_targets: List[KnockoutTarget], overexpression_targets: List[OverexpressionTarget]) \
            -> BooleanModelState:
        """As default all the neutral state targets are set to True. All other state targets as well
        as all reaction targets are set to False."""
        conds = {}  # type: Dict[Target, bool]

        for reaction_target in reaction_targets:
            conds[reaction_target] = False

        for state_target in state_targets:
            # Neutral state targets are True.
            if state_target.is_neutral:
                conds[state_target] = True
            # All reaction targets and non-neutral state targets are False.
            else:
                conds[state_target] = False

        for knockout_target in knockout_targets:
            conds[knockout_target] = False

        for overexpression_target in overexpression_targets:
            conds[overexpression_target] = False

        return BooleanModelState(conds)

    def calc_component_presence_factors() -> Tuple[Dict[Spec, VennSet[StateTarget]], List[ComponentStateTarget]]:
        """The form of the component presence factor is:
            (state_a1 | ... | state_an) & (state_b1 | ... | state_bm) & ...

        Mutually exclusive states are combined by boolean OR (state_a1 ... state_an , state_b1 ... state_bm).
        These ORs are then combines with ANDs.

        If a component does not carry states, this will be a ComponentStateTarget.
        """
        component_state_targets = []  # type: List[ComponentStateTarget]
        component_to_factor = {}  # type: Dict[Spec, VennSet[StateTarget]]
        for component in rxncon_sys.components():
            grouped_states = rxncon_sys.states_for_component_grouped(component)
            # component is not part of any state
            if not grouped_states.values():
                component_state_targets.append(ComponentStateTarget(component))
                component_to_factor[component] = ValueSet(ComponentStateTarget(component))
            # component is part of at least one state
            else:
                # mutually exclusive states are combined by OR

                component_to_factor[component] = \
                    Intersection(
                        *(Union(*(ValueSet(StateTarget(x)) for x in group)) for group in grouped_states.values()))

        return component_to_factor, component_state_targets

    def calc_reaction_targets_with_dnf_contingencies(k_plus_strict: bool, k_minus_strict: bool) -> List[ReactionTarget]:
        """Calculates contingency factors for reaction targets.
        Degradation reactions are handled differently then other reactions. An OR contingency will lead to a
        split of the degradation reaction in as many reactions as OR statements. Each OR will be assigned to one
        instance of the reaction."""
        reaction_targets = set()

        for reaction in rxncon_sys.reactions:
            factors = (x.to_venn_set(k_plus_strict=k_plus_strict, k_minus_strict=k_minus_strict, structured=False,
                                     state_wrapper=StateTarget)
                       for x in rxncon_sys.contingencies_for_reaction(reaction))
            cont = Intersection(*factors).to_simplified_set()  # type: VennSet[StateTarget]
            # The reaction is not a degradation reaction or the DNF has just one term.
            if not reaction.degraded_components or len(cont.to_dnf_list()) == 1:
                reaction_targets.add(ReactionTarget(reaction, contingency_factor=cont))
            # The reaction is a degradation reaction
            else:
                # The reaction is split into separated entities according to the number of minterms of the
                # disjunctive normal form (dnf). Each minterm will be assigned to a entity of the degradation reaction.
                for index, factor in enumerate(cont.to_dnf_list()):
                    reaction_targets.add(
                        ReactionTarget(reaction, contingency_variant=index, contingency_factor=factor))

        return list(reaction_targets)

    def update_degs_add_component_states(reaction_targets: List[ReactionTarget],
                                         component_state_targets: List[ComponentStateTarget]) -> List[ReactionTarget]:
        """For degradation reactions, add the stateless components they degrade to the list of targets they degrade."""
        result = deepcopy(reaction_targets)
        for reaction_target in result:
            for degraded_component in reaction_target.degraded_components:
                if ComponentStateTarget(degraded_component) in component_state_targets:
                    reaction_target.degraded_targets.append(ComponentStateTarget(degraded_component))

        return result

    def update_degs_add_contingent_states(reaction_targets: List[ReactionTarget]) -> List[ReactionTarget]:
        """For degradation reactions, add their contingent states to the list of targets they degrade."""
        def degraded_state_targets(component: Spec, soln: Dict[StateTarget, bool]) -> List[StateTarget]:
            # Disregard input states, since they do not influence which states are degraded.
            soln = {k: v for k, v in soln.items() if not k.is_input() and component in k.components}

            # soln evaluates to False if solution is tautology, since when there are no constraints on which
            # states are required to be true/false, soln is an empty dict. Nicely counterintuitive.
            if not soln and ComponentStateTarget(component) in component_state_targets:
                return [ComponentStateTarget(component)]
            elif not soln:
                return [StateTarget(x) for x in rxncon_sys.states_for_component(component)]
            else:
                trues = [target for target, val in soln.items() if val]
                falses = [target for target, val in soln.items() if not val]
                for target in falses:
                    trues += target.complementary_state_targets(rxncon_sys, component)

                return trues

        result = deepcopy(reaction_targets)

        for reaction_target in result:
            solutions = reaction_target.contingency_factor.calc_solutions()

            if reaction_target.degraded_targets and len(solutions) == 1 and not solutions[0]:
                # No contingencies, but targeted degradation. Do not mess with the list of degraded targets.
                continue

            for degraded_component, solution in product(reaction_target.degraded_components, solutions):  # type: ignore
                reaction_target.degraded_targets.extend(degraded_state_targets(degraded_component, solution))

        return result

    def update_degs_add_interaction_state_partner(reaction_targets: List[ReactionTarget]) -> List[ReactionTarget]:
        """Update degradation reactions for interaction states.
        Interaction states are composed out of two components. A degradation reaction degrading an interaction
        state will degrade one of these components. The other component is then released, and thus produced
        by the degrading reaction."""
        result = []

        for reaction_target in reaction_targets:
            appended = False
            degraded_interaction_targets = [x for x in reaction_target.degraded_targets if x.is_interaction]
            for index, interaction_target in enumerate(degraded_interaction_targets):
                empty_partners = [neutral_target for neutral_target in interaction_target.neutral_targets
                                  if not any(component in reaction_target.degraded_components
                                             for component in neutral_target.components)]

                if interaction_target.is_homodimer:
                    assert len(empty_partners) == 0, 'update_degs_add_interaction_state_partner::homodimer error, ' \
                                                     'reaction_target: {}, interaction_target: {}, empty_partners: {}' \
                                                     ''.format(reaction_target, interaction_target, empty_partners)
                    continue

                if len(empty_partners) != 1:
                    raise AssertionError('update_degs_add_interaction_state_partner::empty partners != 1 error\n'
                                         'The full list of degraded targets is {}\n'
                                         'The current reaction target is {}\n'
                                         'The current interaction target is {}\n'
                                         'The current empty partners that have been deduced are: {}\n'
                                         .format(', '.join(str(x) for x in degraded_interaction_targets),
                                                 str(reaction_target), str(interaction_target),
                                                 ', '.join(str(x) for x in empty_partners)))
                new_reaction = deepcopy(reaction_target)
                new_reaction.interaction_variant_index = index if len(degraded_interaction_targets) > 1 else None
                new_reaction.consumed_targets.append(interaction_target)
                new_reaction.produced_targets.append(empty_partners[0])
                result.append(new_reaction)
                appended = True

            if not appended:
                result.append(deepcopy(reaction_target))

        return result

    def update_syns_with_component_states(reaction_targets: List[ReactionTarget],
                                          component_state_targets: List[ComponentStateTarget]) -> List[ReactionTarget]:
        """Update synthesis reaction with component states: stateless components that are synthesised have
        rights too."""
        result = deepcopy(reaction_targets)

        for reaction_target in result:
            for component in reaction_target.synthesised_components:
                if ComponentStateTarget(component) in component_state_targets:
                    reaction_target.synthesised_targets.append(ComponentStateTarget(component))

        return result

    def calc_knockout_targets(knockout_strategy: KnockoutStrategy) -> List[KnockoutTarget]:
        if knockout_strategy == KnockoutStrategy.no_knockout:
            return []
        else:
            return [KnockoutTarget(component) for component in rxncon_sys.components()]

    def calc_overexpression_targets(overexpression_strategy: OverexpressionStrategy) -> List[OverexpressionTarget]:
        if overexpression_strategy == OverexpressionStrategy.no_overexpression:
            return []
        else:
            return [OverexpressionTarget(component) for component in rxncon_sys.components()]

    def calc_reaction_rules() -> None:
        """Calculate the rules of reaction targets: the component factor AND the contingencies."""

        for reaction_target in reaction_targets:
            components = (component_presence_factor[x] for x in reaction_target.components_lhs)
            component_factor = Intersection(*components)  # type: VennSet[StateTarget]
            reaction_rules.append(UpdateRule(reaction_target, Intersection(
                component_factor, reaction_target.contingency_factor).to_simplified_set()))

    def calc_state_rules() -> None:
        """Calculate the rules of the state targets, includes smoothing. For details see our paper."""
        def synthesis_factor(state_target: StateTarget) -> VennSet[Target]:
            fac = EmptySet()  # type: VennSet[Target]
            for rxn in (x for x in reaction_targets if x.synthesises(state_target)):
                fac = Union(fac, ValueSet(rxn))

            for prod_rxn in (x for x in reaction_targets if x.produces(state_target)):
                sources = []
                for source in prod_rxn.consumed_targets:
                    sources.append([source] + [x for x in reaction_targets if x.synthesises(source)])

                for source_combi in product(*sources):
                    # At least one source should be synthesised.
                    if all(isinstance(x, StateTarget) for x in source_combi):
                        continue
                    assert any(isinstance(x, ReactionTarget) and x.synthesised_targets for x in source_combi)

                    fac = Union(fac, Intersection(ValueSet(prod_rxn), *(ValueSet(x) for x in source_combi)))

            return fac

        def component_factor(state_target: StateTarget) -> VennSet[StateTarget]:
            return Intersection(*(component_presence_factor[x] for x in state_target.components))

        def degradation_factor(state_target: StateTarget) -> VennSet[ReactionTarget]:
            return Complement(Union(*(ValueSet(x) for x in reaction_targets if x.degrades(state_target))))

        def pi(state_target: StateTarget, level: int) -> VennSet[Target]:
            res = EmptySet()  # type: VennSet[Target]

            for r in (x for x in reaction_targets if x.produces(state_target)):
                rxn_term = ValueSet(r)  # type: VennSet[Target]
                for s in (x for x in state_targets if r.consumes(x)):
                    if r.degraded_targets:
                        state_term = ValueSet(s)  # type: VennSet[Target]
                    else:
                        state_term = Intersection(ValueSet(s), degradation_factor(s))
                    for l in range(level):
                        state_term = Union(state_term, sigma(s, level - 1))
                    rxn_term = Intersection(rxn_term, state_term)
                res = Union(res, rxn_term)

            return res

        def kappa(state_target: StateTarget, level: int) -> VennSet[Target]:
            res = EmptySet()  # type: VennSet[Target]

            for r in (x for x in reaction_targets if x.consumes(state_target)):
                rxn_term = ValueSet(r)  # type: VennSet[Target]
                for s in (x for x in state_targets if r.consumes(x)):
                    rxn_term = Intersection(rxn_term, ValueSet(s), degradation_factor(s))
                res = Union(res, rxn_term)

            return res

        def sigma(state_target: StateTarget, level: int) -> VennSet[Target]:
            prod_cons_factor = Union(pi(state_target, level),
                                     Intersection(ValueSet(state_target), Complement(kappa(state_target, level))))

            return Union(synthesis_factor(state_target),
                         Intersection(degradation_factor(state_target),
                                      component_factor(state_target),
                                      prod_cons_factor))

        if smoothing_strategy == SmoothingStrategy.no_smoothing:
            level = 0
        elif smoothing_strategy == SmoothingStrategy.smooth_production_sources:
            level = 1
        else:
            raise AssertionError

        for state_target in state_targets:
            state_rules.append(UpdateRule(state_target, sigma(state_target, level).to_simplified_set()))

    def update_state_rules_with_knockouts(knockout_strategy: KnockoutStrategy) -> None:
        if knockout_strategy == KnockoutStrategy.no_knockout:
            return
        elif knockout_strategy in (KnockoutStrategy.knockout_all_states, KnockoutStrategy.knockout_neutral_states):
            for state_rule in state_rules:
                assert isinstance(state_rule.target, StateTarget)

                if knockout_strategy == KnockoutStrategy.knockout_neutral_states and not state_rule.target.is_neutral:
                    continue

                knockout_factor = Complement(
                    Union(*(ValueSet(KnockoutTarget(component)) for component in state_rule.target.components)))
                state_rule.factor = Intersection(knockout_factor, state_rule.factor)

    def update_state_rules_with_overexpressions(overexpression_strategy: OverexpressionStrategy) -> None:
        if overexpression_strategy == OverexpressionStrategy.no_overexpression:
            return
        elif overexpression_strategy in (OverexpressionStrategy.overexpress_all_states,
                                         OverexpressionStrategy.overexpress_neutral_states):
            for state_rule in state_rules:
                assert isinstance(state_rule.target, StateTarget)

                if overexpression_strategy == OverexpressionStrategy.overexpress_neutral_states and not state_rule.target.is_neutral:
                    continue

                overexpression_factor = Intersection(
                    *(ValueSet(OverexpressionTarget(component)) for component in state_rule.target.components))
                state_rule.factor = Union(overexpression_factor, state_rule.factor)

    def calc_knockout_rules() -> None:
        for knockout_target in knockout_targets:
            knockout_rules.append(UpdateRule(knockout_target, ValueSet(knockout_target)))

    def calc_overexpression_rules() -> None:
        for overexpression_target in overexpression_targets:
            overexpression_rules.append(UpdateRule(overexpression_target, ValueSet(overexpression_target)))

    def update_input_output_rules() -> None:
        """If an Input state and an Output reaction share the same name [BLA], they are assumed
        to refer to the same global quantity. Therefore the update rule for the state (which was trivial),
        becomes the update rule for the reaction."""
        to_delete = []  # type: List[UpdateRule]
        for reaction_rule in reaction_rules:
            for state_rule in state_rules:
                if (reaction_rule.target.is_output and state_rule.target.is_input and   # type: ignore
                        str(reaction_rule.target) == str(state_rule.target)):
                    state_rule.factor = reaction_rule.factor
                    to_delete.append(reaction_rule)

        for rule_to_delete in to_delete:
            reaction_targets.remove(rule_to_delete.target)
            reaction_rules.remove(rule_to_delete)

    component_presence_factor, component_state_targets = calc_component_presence_factors()

    state_targets = [StateTarget(x) for x in rxncon_sys.states]  # type: List[StateTarget]
    state_targets += component_state_targets

    reaction_targets = calc_reaction_targets_with_dnf_contingencies(k_plus_strict, k_minus_strict)
    reaction_targets = update_degs_add_component_states(reaction_targets, component_state_targets)
    reaction_targets = update_degs_add_contingent_states(reaction_targets)
    reaction_targets = update_degs_add_interaction_state_partner(reaction_targets)
    reaction_targets = update_syns_with_component_states(reaction_targets, component_state_targets)

    knockout_targets = calc_knockout_targets(knockout_strategy)
    overexpression_targets = calc_overexpression_targets(overexpression_strategy)

    reaction_rules = []  # type: List[UpdateRule]
    state_rules = []  # type: List[UpdateRule]
    knockout_rules = []  # type: List[UpdateRule]
    overexpression_rules = []  # type: List[UpdateRule]

    calc_reaction_rules()
    calc_state_rules()
    update_state_rules_with_knockouts(knockout_strategy)
    update_state_rules_with_overexpressions(overexpression_strategy)
    calc_knockout_rules()
    calc_overexpression_rules()
    update_input_output_rules()

    return BooleanModel(state_targets + reaction_targets + knockout_targets + overexpression_targets,  # type: ignore
                        reaction_rules + state_rules + knockout_rules + overexpression_rules,
                        initial_conditions(reaction_targets, state_targets, knockout_targets, overexpression_targets))
