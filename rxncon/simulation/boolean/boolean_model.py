from copy import deepcopy
from enum import Enum
from itertools import product
from typing import List, Dict, Tuple, Optional
from abc import ABCMeta

from rxncon.core.reaction import Reaction, OutputReaction
from rxncon.core.rxncon_system import RxnConSystem
from rxncon.core.spec import Spec
from rxncon.core.state import State, InteractionState
from rxncon.venntastic.sets import Set as VennSet, ValueSet, Intersection, Union, Complement, UniversalSet


class BooleanModel:
    def __init__(self, targets: List['Target'], update_rules: List['UpdateRule'], initial_conditions: 'BooleanModelConfig') -> None:
        self.update_rules            = sorted(update_rules)
        self.initial_conditions      = initial_conditions
        self._state_targets          = {str(x): x for x in targets if isinstance(x, StateTarget)}
        self._reaction_targets       = {str(x): x for x in targets if isinstance(x, ReactionTarget)}
        self._knockout_targets       = {str(x): x for x in targets if isinstance(x, KnockoutTarget)}
        self._overexpression_targets = {str(x): x for x in targets if isinstance(x, OverexpressionTarget)}
        self._validate_update_rules()
        self._validate_initial_conditions()

    def set_initial_condition(self, target: 'Target', value: bool) -> None:
        self.initial_conditions.set_target(target, value)

    def state_target_by_name(self, name: str) -> 'StateTarget':
        return self._state_targets[name]

    def reaction_target_by_name(self, name: str) -> 'ReactionTarget':
        return self._reaction_targets[name]

    def knockout_target_by_name(self, name: str) -> 'KnockoutTarget':
        return self._knockout_targets[name]

    def overexpression_target_by_name(self, name: str) -> 'OverexpressionTarget':
        return self._overexpression_targets[name]

    def _validate_update_rules(self) -> None:
        """
        Validating the update rules.

        Note:
            All targets in the factor expressions have to be regulated somehow.

        Returns:
            None

        Raises:
            AssertionError: Raise an error if not all targets on the factor side are also targets on the target side.

        """
        all_lhs_targets = []  # type: List[Target]
        all_rhs_targets = []  # type: List[Target]
        for rule in self.update_rules:
            all_lhs_targets.append(rule.target)
            all_rhs_targets += rule.factor_targets

        assert all(x in all_lhs_targets for x in all_rhs_targets)

    def _validate_initial_conditions(self) -> None:
        self.initial_conditions.validate_by_model(self)


class BooleanModelConfig:
    """
    Configuration of the boolean model

    Args:
        target_to_value: Mapping of state and reaction targets to specific boolean values.

    """
    def __init__(self, target_to_value: Dict['Target', bool]) -> None:
        self.target_to_value = target_to_value

    def set_target(self, target: 'Target', value: bool) -> None:
        """
        Assigning a value to a target.

        Args:
            target: StateTarget or ReactionTarget.
            value: boolean value.

        Mutates:
            target_to_value: Mapping of state and reaction targets to specific boolean values.

        Returns:
            None

        """
        self.target_to_value[target] = value

    def validate_by_model(self, model: BooleanModel) -> None:
        """
        Validating the boolean model.

        Args:
            model: boolean model

        Returns:
            None

        Raises:
            AssertionError: A error will be raised if there are not the same targets in the model and the configuration.

        """
        model_targets  = [rule.target for rule in model.update_rules]
        config_targets = self.target_to_value.keys()

        assert set(model_targets) == set(config_targets) and len(model_targets) == len(config_targets)


class Target(metaclass=ABCMeta):  # pylint: disable=too-few-public-methods
    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)


class ReactionTarget(Target):  # pylint: disable=too-many-instance-attributes
    """
    Reaction of the boolean model.

    Args:
        reaction_parent: A elemental reaction of the rxncon system.

    """
    def __init__(self, reaction_parent: Reaction, contingency_variant: Optional[int]=None,
                 interaction_variant: Optional[int]=None, contingency_factor: VennSet['StateTarget']=None) -> None:
        self.reaction_parent     = reaction_parent  # type: Reaction
        self.produced_targets    = [StateTarget(x) for x in reaction_parent.produced_states]     # type: List[StateTarget]
        self.consumed_targets    = [StateTarget(x) for x in reaction_parent.consumed_states]     # type: List[StateTarget]
        self.synthesised_targets = [StateTarget(x) for x in reaction_parent.synthesised_states]  # type: List[StateTarget]
        self.degraded_targets    = [StateTarget(x) for x in reaction_parent.degraded_states]     # type: List[StateTarget]

        self.contingency_variant_index = contingency_variant
        self.interaction_variant_index = interaction_variant

        if contingency_factor is None:
            self.contingency_factor = UniversalSet()      # type: VennSet[StateTarget]
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
        """
        Asking for components of the left hand side.

        Returns:
            List of components.

        """
        return self.reaction_parent.components_lhs

    @property
    def components_rhs(self) -> List[Spec]:
        """
        Asking for components of the right hand side.

        Returns:
            List of components.

        """
        return self.reaction_parent.components_rhs

    @property
    def degraded_components(self) -> List[Spec]:
        """
        Asking for components getting degraded.

        Returns:
            List of components.

        """
        return [component for component in self.components_lhs if component not in self.components_rhs]

    @property
    def synthesised_components(self) -> List[Spec]:
        """
        Asking for components getting synthesised.

        Returns:
            List of components.

        """
        return [component for component in self.components_rhs if component not in self.components_lhs]

    def degrades_component(self, spec: Spec) -> bool:
        """
        Asking if a component is degraded.

        Args:
            spec: Specification

        Returns:
            bool: True if the component is in the list of degraded components, False otherwise.

        """
        assert spec.is_component_spec
        return spec in self.degraded_components

    def synthesises_component(self, spec: Spec) -> bool:
        """
        Asking if a component is synthesised.

        Args:
            spec: Specification

        Returns:
            bool: True if the component is in the list of synthesised components, False otherwise.

        """
        assert spec.is_component_spec
        return spec in self.synthesised_components

    def is_output(self) -> bool:
        """
        Checks if the ReactionTarget is an OUTPUT

        Returns:
            bool

        """
        return isinstance(self.reaction_parent, OutputReaction)


class StateTarget(Target):
    """
    An elemental state of the boolean model.

    Args:
        state_parent: A state of the rxncon system.

    """
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
        """
        Checks if the state target is degraded by the respective reaction target.

        Args:
            reaction_target: Reactions of the boolean model producing, consuming, degrading or synthesising state targets.

        Returns:
            bool: True if the state target is produced by the respective reaction target, False otherwise.

        """
        return reaction_target.produces(self)

    def is_consumed_by(self, reaction_target: ReactionTarget) -> bool:
        """
        Checks if the state target is consumed by the respective reaction target.

        Args:
            reaction_target: Reactions of the boolean model producing, consuming, degrading or synthesising state targets.

        Returns:
            bool: True if the state target is consumed by the respective reaction target, False otherwise.

        """
        return reaction_target.consumes(self)

    def is_synthesised_by(self, reaction_target: ReactionTarget) -> bool:
        """
        Checks if the state target is synthesised by the respective reaction target.

        Args:
            reaction_target: Reactions of the boolean model producing, consuming, degrading or synthesising state targets.

        Returns:
            bool: True if the state target is synthesised by the respective reaction, False otherwise.

        """
        return reaction_target.synthesises(self)

    def is_degraded_by(self, reaction_target: ReactionTarget) -> bool:
        """
        Checks if the state target is degraded by the respective reaction target.

        Args:
            reaction_target: Reactions of the boolean model producing, consuming, degrading or synthesising state targets.

        Returns:
            bool: True if state target is degraded by the respective reaction, False otherwise.

        """
        return reaction_target.degrades(self)

    def is_input(self) -> bool:
        """
        Checks if the StateTarget is an INPUT

        Returns:
            bool

        """
        return self.state_parent.is_global

    @property
    def components(self) -> List[Spec]:
        """
        Asking for the components of the state target.

        Returns:
            List of components.

        """
        return self.state_parent.components

    def shares_component_with(self, other_target: 'StateTarget') -> bool:
        return any(x in other_target.components for x in self.components)

    @property
    def is_neutral(self) -> bool:
        """
        Asking for the neutrality of the state target.

        Returns:
            bool: True if neutral, False otherwise.

        """
        return self.state_parent.is_neutral

    @property
    def is_homodimer(self) -> bool:
        return self.state_parent.is_homodimer

    @property
    def neutral_targets(self) -> List['StateTarget']:
        """
        Asking for the neutral states of state target.

        Returns:
            List of neutral StateTargets.

        """
        return [StateTarget(x) for x in self.state_parent.neutral_states]

    @property
    def is_interaction(self) -> bool:
        return isinstance(self.state_parent, InteractionState)

    def is_mutually_exclusive_with(self, other: 'StateTarget') -> bool:
        """
        Asking if a state is mutually exclusive with another state.

        Args:
            other: StateTarget

        Returns:
            bool: True if state is mutually exclusive, False otherwise.

        """
        return self.state_parent.is_mutually_exclusive_with(other.state_parent)

    def complementary_state_targets(self, rxnconsys: RxnConSystem, component: Spec) -> List['StateTarget']:
        others = rxnconsys.complement_states_for_component(component, self.state_parent)
        return [StateTarget(x) for x in others]


class ComponentStateTarget(StateTarget):
    def __init__(self, component: Spec) -> None:  # pylint: disable=super-init-not-called
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
        """
        Asking for the components of the component state target.

        Returns:
            List of components

        """
        return [self.component]

    @property
    def is_neutral(self) -> bool:
        """
        Asking for the neutrality of the component State target.

        Note:
            ComponentStateTargets are always neutral.

        Returns:
            bool: True

        """
        return True

    @property
    def is_interaction(self) -> bool:
        return False


class KnockoutTarget(ComponentStateTarget):
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
        """
        Asking for the state and reaction targets, which are important for the update of target.

        Returns:
            A list of Targets.

        """
        return self.factor.values


class SmoothingStrategy(Enum):
    """
    Defined smoothing strategies.

    Note:
        To overcome non biological oscillatory behaviour during the simulation we introduced a smoothing strategy.

    """
    no_smoothing              = 'no_smoothing'
    smooth_production_sources = 'smooth_production_sources'


class KnockoutStrategy(Enum):
    no_knockout             = 'no_knockout'
    knockout_neutral_states = 'knockout_neutral_states'
    knockout_all_states     = 'knockout_all_states'


class OverexpressionStrategy(Enum):
    no_overexpression          = 'no_overexpression'
    overexpress_neutral_states = 'overexpress_neutral_states'
    overexpress_all_states     = 'overexpress_all_states'


def boolean_model_from_rxncon(rxncon_sys: RxnConSystem,
                              smoothing_strategy: SmoothingStrategy=SmoothingStrategy.no_smoothing,
                              knockout_strategy: KnockoutStrategy=KnockoutStrategy.no_knockout,
                              overexpression_strategy: OverexpressionStrategy=OverexpressionStrategy.no_overexpression,
                              k_plus_strict: bool=True, k_minus_strict: bool=True) -> BooleanModel:
    """
    Constructs a boolean model from a rxncon system and a smoothing strategy.

    Args:
          rxncon_sys: The rxncon system.
          smoothing_strategy: The smoothing strategy. Defaults to no smoothing.

    Returns:
          The boolean model.

    """
    def initial_conditions(reaction_targets: List[ReactionTarget], state_targets: List[StateTarget],
                           knockout_targets: List[KnockoutTarget], overexpression_targets: List[OverexpressionTarget])\
            -> BooleanModelConfig:
        """
        Calculates default initial conditions of the boolean model.

        Note:
            As default all the neutral state targets are set to True. All other state targets as well as all
            reaction targets are set to False.

        Args:
            reaction_targets: Reactions of the boolean model producing or consuming state targets.
            state_targets: States of the boolean model consumed, produced by reaction targets or regulating reaction targets.

        Returns:
            The initial conditions of the boolean model.

        """
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

        return BooleanModelConfig(conds)

    def calc_component_presence_factors() -> Tuple[Dict[Spec, VennSet[StateTarget]], List[ComponentStateTarget]]:
        """
        Calculates the factors for components.

        Note:
            The form is: (state_a1 | ... | state_an) & (state_b1 | ... | state_bm) & ...

            Non-mutually exclusive states are combined by boolean AND (state_a and state_b).
            Mutually exclusive states are combined by boolean OR (state_a1 to state_an as well as state_b1 to state_bm).

            If a component is not part of any state of the system, the component will hold itself as ValueSet of a ComponentStateTarget.

        Mutates:
            component_to_factor: Mapping of components and VennSets, containing all the states the component is
                involved in.

        Returns:
            None

        """
        component_state_targets = []  # type: List[ComponentStateTarget]
        component_to_factor     = {}  # type: Dict[Spec, VennSet[StateTarget]]
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
                    Intersection(*(Union(*(ValueSet(StateTarget(x)) for x in group)) for group in grouped_states.values()))

        return component_to_factor, component_state_targets

    def calc_reaction_targets_with_dnf_contingencies(k_plus_strict: bool, k_minus_strict: bool) -> List[ReactionTarget]:
        """
        Calculates contingency factors for reaction targets.

        Note: Degradation reactions are handled differently then other reactions. An OR contingency will lead to a
            split of the degradation reaction in as many reactions as OR statements. Each OR will be assigned to one
            instance of the reaction.

        Mutates:
            reaction_target_to_factor: Mapping of target reactions and their corresponding contingencies.

        Returns:
            None
        """

        reaction_targets = []

        for reaction in rxncon_sys.reactions:
            factors = (x.to_venn_set(k_plus_strict=k_plus_strict, k_minus_strict=k_minus_strict, structured=False, state_wrapper=StateTarget)
                       for x in rxncon_sys.contingencies_for_reaction(reaction))
            cont = Intersection(*factors).to_simplified_set()  # type: VennSet[StateTarget]
            # The reaction is not a degradation reaction or the DNF has just one term.
            if not reaction.degraded_components or len(cont.to_dnf_list()) == 1:
                reaction_targets.append(ReactionTarget(reaction, contingency_factor=cont))
            # The reaction is a degradation reaction
            else:
                # The reaction is split into separated entities according to the number of minterms of the
                # disjunctive normal form (dnf). Each minterm will be assigned to a entity of the degradation reaction.
                for index, factor in enumerate(cont.to_dnf_list()):
                    reaction_targets.append(ReactionTarget(reaction, contingency_variant=index, contingency_factor=factor))

        return reaction_targets

    def update_degs_add_component_states(reaction_targets: List[ReactionTarget],
                                         component_state_targets: List[ComponentStateTarget]) -> List[ReactionTarget]:
        result = deepcopy(reaction_targets)
        for reaction_target in result:
            for degraded_component in reaction_target.degraded_components:
                if ComponentStateTarget(degraded_component) in component_state_targets:
                    reaction_target.degraded_targets.append(ComponentStateTarget(degraded_component))

        return result

    def update_degs_add_contingent_states(reaction_targets: List[ReactionTarget]) -> List[ReactionTarget]:
        def degraded_state_targets(component: Spec, soln: Dict[StateTarget, bool]) -> List[StateTarget]:
            # soln evaluates to False if solution is tautology.
            if not soln and ComponentStateTarget(component) in component_state_targets:
                return [ComponentStateTarget(component)]
            elif not soln:
                return [StateTarget(x) for x in rxncon_sys.states_for_component(component)]
            else:
                trues  = [target for target, val in soln.items() if val]
                falses = [target for target, val in soln.items() if not val]
                for target in falses:
                    trues += target.complementary_state_targets(rxncon_sys, component)
                return trues

        result = deepcopy(reaction_targets)

        for reaction_target in result:
            solutions = reaction_target.contingency_factor.calc_solutions()

            for degraded_component, solution in product(reaction_target.degraded_components, solutions):  # type: ignore
                reaction_target.degraded_targets.extend(degraded_state_targets(degraded_component, solution))

        return result

    def update_degs_add_interaction_state_partner(reaction_targets: List[ReactionTarget]) -> List[ReactionTarget]:
        """
        Update degradation reactions for interaction states.

        Note: Interaction states are composed out of two components. A degradation reaction degrading an interaction
            state will degrade one of these components. The other component is assigned as unbound (neutral form).

        Mutates:
            reaction_target_to_factor: Mapping of target reactions and their corresponding contingencies.

        Returns:
            None

        """

        result = []

        for reaction_target in reaction_targets:
            appended = False
            degraded_interaction_targets = [x for x in reaction_target.degraded_targets if x.is_interaction]
            for index, interaction_target in enumerate(degraded_interaction_targets):
                empty_partners = [neutral_target for neutral_target in interaction_target.neutral_targets
                                  if not any(component in reaction_target.degraded_components for component in neutral_target.components)]

                if interaction_target.is_homodimer:
                    assert len(empty_partners) == 0
                    continue

                assert len(empty_partners) == 1
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
        """
        Update synthesis reaction with component states.

        Mutates:
            reaction_target: Reaction of the boolean model producing, consuming, degrading or synthesising state targets.

        Returns:
            None

        """

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
        """
        Calculate the rules of reaction targets.

        Note:
            The factor of a reaction target has the form: components AND contingencies.

        Mutates:
            reaction_rules: Containing the rules of the boolean model

        Returns:
            None

        """

        for reaction_target in reaction_targets:
            components = (component_presence_factor[x] for x in reaction_target.components_lhs)
            component_factor = Intersection(*components)  # type: VennSet[StateTarget]
            reaction_rules.append(UpdateRule(reaction_target, Intersection(
                component_factor, reaction_target.contingency_factor).to_simplified_set()))

    def calc_state_rules() -> None:
        """
        Calculates the rules of state targets.

        Note:
            The factor for a state target has the from:
            synthesis OR (components AND NOT degradation AND ((production AND sources) OR (state AND NOT (consumption AND sources))))

        Returns:
            A list of updating rules for states.

      """

        def reaction_with_sources(reaction_target: ReactionTarget) -> VennSet[Target]:
            """
            Calculates the source states of the respective reaction target

            Args:
                reaction_target: Reaction of the boolean model producing, consuming, degrading or synthesising state targets.

            Returns:
                VennSet of the reaction target and its source states

            """
            sources = Intersection(*(ValueSet(x) for x in reaction_target.consumed_targets))
            return Intersection(ValueSet(reaction_target), sources)

        def indirect_synth_path(state_target: StateTarget) -> VennSet[ReactionTarget]:
            my_brothers = [x for x in state_targets if state_target.shares_component_with(x) and x != state_target]
            return Union(*(ValueSet(rxn) for state in my_brothers for rxn in reaction_targets if rxn.synthesises(state)))

        def synthesis_factor(state_target: StateTarget) -> VennSet[ReactionTarget]:
            return Union(*(ValueSet(x) for x in reaction_targets if x.synthesises(state_target)))

        def component_factor(state_target: StateTarget) -> VennSet[StateTarget]:
            return Intersection(*(component_presence_factor[x] for x in state_target.components))

        def degradation_factor(state_target: StateTarget) -> VennSet[ReactionTarget]:
            return Complement(Union(*(ValueSet(x) for x in reaction_targets if x.degrades(state_target))))

        for state_target in state_targets:
            synt_fac = synthesis_factor(state_target)
            comp_fac = component_factor(state_target)
            degr_fac = degradation_factor(state_target)

            prod_facs = []  # type: List[VennSet]
            cons_facs = []  # type: List[VennSet]

            for reaction_target in (target for target in reaction_targets if target.produces(state_target)):
                if smoothing_strategy == SmoothingStrategy.no_smoothing:
                    prod_facs.append(reaction_with_sources(reaction_target))
                elif smoothing_strategy == SmoothingStrategy.smooth_production_sources:
                    smoothed_prod_facs = []
                    for primary_source in reaction_target.consumed_targets:
                        smooth_source = Union(ValueSet(primary_source),
                                              Intersection(Union(*(reaction_with_sources(rxn) for rxn in reaction_targets if rxn.produces(primary_source))),
                                                           Union(degradation_factor(primary_source), indirect_synth_path(primary_source))))

                        smoothed_prod_facs.append(smooth_source)
                    prod_facs.append(Intersection(ValueSet(reaction_target), *smoothed_prod_facs))
                else:
                    raise AssertionError

            for reaction_target in (target for target in reaction_targets if target.consumes(state_target)):
                cons_facs.append(Complement(reaction_with_sources(reaction_target)))

            tot_prod_fac = Intersection(comp_fac, Union(*prod_facs), Union(degr_fac, indirect_synth_path(state_target)))
            tot_cons_fac = Intersection(comp_fac, ValueSet(state_target), Intersection(*cons_facs), degr_fac)

            state_rules.append(UpdateRule(state_target,
                                          Union(synt_fac,
                                                tot_prod_fac,
                                                tot_cons_fac).to_simplified_set()))

    def update_state_rules_with_knockouts(knockout_strategy: KnockoutStrategy) -> None:
        if knockout_strategy == KnockoutStrategy.no_knockout:
            return
        elif knockout_strategy in (KnockoutStrategy.knockout_all_states, KnockoutStrategy.knockout_neutral_states):
            for state_rule in state_rules:
                assert isinstance(state_rule.target, StateTarget)

                if knockout_strategy == KnockoutStrategy.knockout_neutral_states and not state_rule.target.is_neutral:
                    continue

                knockout_factor = Complement(Union(*(ValueSet(KnockoutTarget(component)) for component in state_rule.target.components)))
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

                overexpression_factor = Intersection(*(ValueSet(OverexpressionTarget(component)) for component in state_rule.target.components))
                state_rule.factor = Union(overexpression_factor, state_rule.factor)

    def calc_knockout_rules() -> None:
        for knockout_target in knockout_targets:
            knockout_rules.append(UpdateRule(knockout_target, ValueSet(knockout_target)))

    def calc_overexpression_rules() -> None:
        for overexpression_target in overexpression_targets:
            overexpression_rules.append(UpdateRule(overexpression_target, ValueSet(overexpression_target)))

    def update_input_output_rules() -> None:
        """
        Updating the input update rule by a output update rule if both target names are matching.

        Mutates:
            reaction_targets: List of reaction targets
            reaction_rules: List of update rules of reaction targets
            state_rules: List of update rules of state targets

        Returns:
            None

        """
        def remove_output_reaction_rules_and_reaction_targets():
            for rule_to_delete in to_delete:
                reaction_targets.remove(rule_to_delete.target)
                reaction_rules.remove(rule_to_delete)

        to_delete = []  # type: List[UpdateRule]
        for reaction_rule in reaction_rules:
            for state_rule in state_rules:
                if reaction_rule.target.is_output and state_rule.target.is_input and \
                                str(reaction_rule.target) == str(state_rule.target):
                    state_rule.factor = reaction_rule.factor
                    to_delete.append(reaction_rule)

        remove_output_reaction_rules_and_reaction_targets()

    component_presence_factor, component_state_targets = calc_component_presence_factors()

    state_targets    = [StateTarget(x) for x in rxncon_sys.states]  # type: List[StateTarget]
    state_targets   += component_state_targets

    reaction_targets = calc_reaction_targets_with_dnf_contingencies(k_plus_strict, k_minus_strict)
    reaction_targets = update_degs_add_component_states(reaction_targets, component_state_targets)
    reaction_targets = update_degs_add_contingent_states(reaction_targets)
    reaction_targets = update_degs_add_interaction_state_partner(reaction_targets)
    reaction_targets = update_syns_with_component_states(reaction_targets, component_state_targets)

    knockout_targets       = calc_knockout_targets(knockout_strategy)
    overexpression_targets = calc_overexpression_targets(overexpression_strategy)

    reaction_rules       = []  # type: List[UpdateRule]
    state_rules          = []  # type: List[UpdateRule]
    knockout_rules       = []  # type: List[UpdateRule]
    overexpression_rules = []  # type: List[UpdateRule]

    calc_reaction_rules()
    calc_state_rules()
    update_state_rules_with_knockouts(knockout_strategy)
    update_state_rules_with_overexpressions(overexpression_strategy)
    calc_knockout_rules()
    calc_overexpression_rules()
    update_input_output_rules()

    return BooleanModel(state_targets + reaction_targets + knockout_targets + overexpression_targets,
                        reaction_rules + state_rules + knockout_rules + overexpression_rules,
                        initial_conditions(reaction_targets, state_targets, knockout_targets, overexpression_targets))
