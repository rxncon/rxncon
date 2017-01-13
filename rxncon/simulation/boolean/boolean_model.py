from typing import List, Dict, Tuple, Iterable
from copy import deepcopy
from enum import Enum

from rxncon.venntastic.sets import Set as VennSet, ValueSet, Intersection, Union, Complement, UniversalSet, EmptySet
from rxncon.core.reaction import Reaction
from rxncon.core.state import State
from rxncon.core.spec import Spec
from rxncon.core.contingency import Contingency, ContingencyType
from rxncon.core.effector import Effector, AndEffector, OrEffector, NotEffector, StateEffector
from rxncon.core.rxncon_system import RxnConSystem


class BooleanModel:
    """
    Definition of the boolean model.

    Args:
        update_rules: Rules for updating the system.
        initial_conditions: Initial conditions of the system.
    """

    def __init__(self, state_targets: List['StateTarget'], reaction_targets: List['ReactionTarget'],
                 update_rules: List['UpdateRule'], initial_conditions: 'BooleanModelConfig') -> None:

        self.update_rules       = update_rules
        self.initial_conditions = initial_conditions
        self._state_targets     = {str(x): x for x in state_targets}
        self._reaction_targets  = {str(x): x for x in reaction_targets}
        self._validate_update_rules()
        self._validate_initial_conditions()

    def set_initial_condition(self, target: 'Target', value: bool) -> None:

        """
        Assigning initial values to the boolean model.

        Args:
            target: StateTarget or ReactionTarget.
            value: boolean value.

        Mutates:
            initial_conditions: Mapping of state and reaction targets to specific boolean values.

        Returns:
            None

        """
        self.initial_conditions.set_target(target, value)

    def state_target_by_name(self, name: str) -> 'StateTarget':
        return self._state_targets[name]

    def reaction_target_by_name(self, name: str) -> 'ReactionTarget':
        return self._reaction_targets[name]

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


class Target:
    """
    Parent class: a ReactionTarget or a StateTarget of the boolean model

    """
    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)


class ReactionTarget(Target):
    """
    Reaction of the boolean model.

    Args:
        reaction_parent: A elemental reaction of the rxncon system.

    """
    def __init__(self, reaction_parent: Reaction) -> None:
        self.reaction_parent     = reaction_parent  # type: Reaction
        self.produced_targets    = [StateTarget(x) for x in reaction_parent.produced_states]  # type: List[StateTarget]
        self.consumed_targets    = [StateTarget(x) for x in reaction_parent.consumed_states]  # type: List[StateTarget]
        self.synthesised_targets = [StateTarget(x) for x in reaction_parent.synthesised_states]  # type: List[StateTarget]
        self.degraded_targets    = [StateTarget(x) for x in reaction_parent.degraded_states]  # type: List[StateTarget]

        self.contingency_variant_index = 0
        self.interaction_variant_index = 0

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
        if self.interaction_variant_index != 0:
            suffix = '#{}/{}'.format(self.contingency_variant_index, self.interaction_variant_index)
        elif self.contingency_variant_index != 0 and self.interaction_variant_index == 0:
            suffix = '#{}'.format(self.contingency_variant_index)

        return str(self.reaction_parent) + suffix

    def __repr__(self) -> str:
        return str(self)

    def produces(self, state_target: 'StateTarget') -> bool:
        """
        Checks if the reaction produces this particular state.

        Args:
            state_target: State of the boolean model consumed, produced, degraded or synthesised by reaction targets or regulating reaction targets.

        Returns:
            bool: True if the state target is in the list of produced targets, False otherwise.

        """
        return state_target in self.produced_targets

    def consumes(self, state_target: 'StateTarget') -> bool:
        """
        Checks if the reaction consumes this particular state.

        Args:
            state_target: State of the boolean model consumed, produced, degraded or synthesised by reaction targets or regulating reaction targets.

        Returns:
            bool: True if the state target is in the list of consumed targets, False otherwise.

        """
        return state_target in self.consumed_targets

    def synthesises(self, state_target: 'StateTarget') -> bool:
        """
        Checks if the reaction synthesises this particular state.

        Args:
            state_target: State of the boolean model consumed, produced, degraded or synthesised by reaction targets or regulating reaction targets.

        Returns:
            bool: True if the state target is in the list of synthesised targets, False otherwise.

        """
        return state_target in self.synthesised_targets

    def degrades(self, state_target: 'StateTarget') -> bool:
        """
        Checks if the reaction degrades this particular state target.

        Args:
            state_target: State of the boolean model consumed, produced, degraded or synthesised by reaction targets or
                regulating reaction targets.

        Returns:
            bool: True if the state target is in the list of degraded targets, False otherwise.

        """
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


class StateTarget(Target):
    """
    An elemental state of the boolean model.

    Args:
        state_parent: A state of the rxncon system.

    """
    def __init__(self, state_parent: State) -> None:
        self._state_parent = state_parent

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return str(self._state_parent)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Target):
            return NotImplemented
        return isinstance(other, StateTarget) and self._state_parent == other._state_parent

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

    @property
    def components(self) -> List[Spec]:
        """
        Asking for the components of the state target.

        Returns:
            List of components.

        """
        return self._state_parent.components

    def shares_component_with(self, other_target: 'StateTarget') -> bool:
        return any(x in other_target.components for x in self.components)

    @property
    def is_neutral(self) -> bool:
        """
        Asking for the neutrality of the state target.

        Returns:
            bool: True if neutral, False otherwise.

        """
        return self._state_parent.is_neutral

    @property
    def neutral_targets(self) -> List['StateTarget']:
        """
        Asking for the neutral states of state target.

        Returns:
            List of neutral StateTargets.

        """
        return [StateTarget(x) for x in self._state_parent.neutral_states]

    def is_mutually_exclusive_with(self, other: 'StateTarget') -> bool:
        """
        Asking if a state is mutually exclusive with another state.

        Args:
            other: StateTarget

        Returns:
            bool: True if state is mutually exclusive, False otherwise.

        """
        return self._state_parent.is_mutually_exclusive_with(other._state_parent)

    def complementary_state_targets(self, rxnconsys: RxnConSystem, component: Spec) -> List['StateTarget']:
        others = rxnconsys.complementary_states_for_component(component, self._state_parent)
        return [StateTarget(x) for x in others]


class ComponentStateTarget(StateTarget):
    """
    ComponentStateTarget corresponding to Components which have no states in the boolean model.

    Args:
        component (Specification): A reaction partner or part of state.

    """
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


class UpdateRule:
    """
    Updating rule of the boolean model.

    Args:
            target: Is a ReactionTarget or StateTarget.
            factor: Is the updating rule for the respective target.

    """
    def __init__(self, target: Target, factor: VennSet) -> None:

        self.target = target
        self.factor = factor

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


def boolean_model_from_rxncon(rxncon_sys: RxnConSystem,
                              smoothing_strategy: SmoothingStrategy=SmoothingStrategy.no_smoothing) -> BooleanModel:
    """
    Constructs a boolean model from a rxncon system and a smoothing strategy.

    Args:
          rxncon_sys: The rxncon system.
          smoothing_strategy: The smoothing strategy. Defaults to no smoothing.

    Returns:
          The boolean model.

    """
    def factor_from_contingency(contingency: Contingency) -> VennSet[StateTarget]:
        """
        Calculates factors from contingency.

        Note:
            Positive quantitative contingencies are handled like required contingencies.
            Negative quantitative contingencies are handled like inhibitions.

        Args:
            contingency: rxncon system contingency (contextual constraint defined on a reaction)

        Returns:
            If there are contingencies a VennSet of StateTargets is returned, the UniversalSet otherwise.

        """
        def parse_effector(eff: Effector) -> VennSet[StateTarget]:
            """
            Reshape VennSet of StateEffectors into VennSet of StateTargets.

            Args:
                eff: VennSet of StateEffectors.

            Returns:
                VennSet of StateTargets.

            Raises:
                AssertionError if effector is not known.

            """
            if isinstance(eff, StateEffector):
                return ValueSet(StateTarget(eff.expr.to_non_structured()))
            elif isinstance(eff, NotEffector):
                return Complement(parse_effector(eff.expr))
            elif isinstance(eff, OrEffector):
                return Union(*(parse_effector(x) for x in eff.exprs))
            elif isinstance(eff, AndEffector):
                return Intersection(*(parse_effector(x) for x in eff.exprs))
            else:
                raise AssertionError

        if contingency.type in [ContingencyType.requirement, ContingencyType.positive]:
            return parse_effector(contingency.effector)
        elif contingency.type in [ContingencyType.inhibition, ContingencyType.negative]:
            return Complement(parse_effector(contingency.effector))
        else:
            return UniversalSet()

    def initial_conditions(reaction_targets: List[ReactionTarget], state_targets: List[StateTarget]) -> BooleanModelConfig:
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

        return BooleanModelConfig(conds)

    def calc_component_factors() -> None:
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

    def calc_contingency_factors() -> None:
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
        for reaction in rxncon_sys.reactions:
            factors = (x.to_venn_set(k_plus_strict=True, k_minus_strict=True, structured=False, state_wrapper=StateTarget)
                       for x in rxncon_sys.contingencies_for_reaction(reaction))
            cont = Intersection(*factors).to_simplified_set()  # type: VennSet[StateTarget]
            # The reaction is not a degradation reaction
            if not reaction.degraded_components:
                reaction_target_to_factor[ReactionTarget(reaction)] = cont
            # The reaction is a degradation reaction
            else:
                # The reaction is split into separated entities according to the number of minterms of the
                # disjunctive normal from (dnf). Each minterm will be assigned to a entity of the degradation reaction.
                for index, factor in enumerate(cont.to_dnf_list()):
                    target = ReactionTarget(reaction)
                    target.contingency_variant_index = index
                    reaction_target_to_factor[target] = factor

    def update_degradations_with_component_states() -> None:
        for reaction_target in reaction_target_to_factor.keys():
            for degraded_component in reaction_target.degraded_components:
                if ComponentStateTarget(degraded_component) in component_state_targets:
                    reaction_target.degraded_targets.append(ComponentStateTarget(degraded_component))

    def update_degradations_with_contingencies() -> None:
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

        for reaction_target, contingency_factor in reaction_target_to_factor.items():
            if contingency_factor.is_equivalent_to(UniversalSet()) and not reaction_target.degraded_components:
                pass

            if contingency_factor.is_equivalent_to(UniversalSet()) and reaction_target.degraded_components:
                for degraded_component in reaction_target.degraded_components:
                    reaction_target.degraded_targets += degraded_state_targets(degraded_component, {})

            # Make sure the contingency concerns the object, not the subject of the degradation.
            contingency_components = []  # type: List[Spec]

            for state in (state for state in contingency_factor.values if state):
                contingency_components.extend(state.components)

            if not any(x in reaction_target.degraded_components for x in contingency_components):
                continue

            solns = contingency_factor.calc_solutions()
            assert len(solns) == 1
            soln = solns[0]
            for degraded_component in reaction_target.degraded_components:
                reaction_target.degraded_targets += degraded_state_targets(degraded_component, soln)

    def update_degradations_for_interaction_states() -> None:
        """
        Update degradation reactions for interaction states.

        Note: Interaction states are composed out of two components. A degradation reaction degrading an interaction
            state will degrade one of these components. The other component is assigned as unbound (neutral form).

        Mutates:
            reaction_target_to_factor: Mapping of target reactions and their corresponding contingencies.

        Returns:
            None

        """

        def state_exclusive_with_contingency(state: StateTarget, contingency_factor: VennSet) -> bool:
            solns = contingency_factor.calc_solutions()
            assert len(solns) == 1
            true_states  = [state for state, val in solns[0].items() if val]
            false_states = [state for state, val in solns[0].items() if not val]

            if state in false_states or any(true_state.is_mutually_exclusive_with(state) for true_state in true_states):
                return True
            else:
                return False

        new_reactions = {}  # type: Dict[ReactionTarget, VennSet]

        for reaction_target, contingency_factor in reaction_target_to_factor.items():
            for num, interaction_state in enumerate(state for state in rxncon_sys.states if len(state.components) > 1
                                                    and any(reaction_target.degrades_component(spec) for spec in state.components)
                                                    and not state_exclusive_with_contingency(StateTarget(state), contingency_factor)):
                neutral_targets = StateTarget(interaction_state).neutral_targets
                new_reaction_target = deepcopy(reaction_target)

                new_reaction_target.consumed_targets.append(StateTarget(interaction_state))
                new_reaction_target.produced_targets += \
                    [x for x in neutral_targets if not any(component in new_reaction_target.degraded_components for component in x.components)]

                new_reaction_target.interaction_variant_index = num + 1
                new_reactions[new_reaction_target] = contingency_factor

        reaction_target_to_factor.update(new_reactions)

    def update_syntheses_with_component_states() -> None:
        """
        Update synthesis reaction with component states.

        Mutates:
            reaction_target: Reaction of the boolean model producing, consuming, degrading or synthesising state targets.

        Returns:
            None

        """
        for reaction_target, _ in reaction_target_to_factor.items():
            for component in reaction_target.synthesised_components:
                if ComponentStateTarget(component) in component_state_targets:
                    reaction_target.synthesised_targets.append(ComponentStateTarget(component))

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

        for reaction_target, contingency_factor in reaction_target_to_factor.items():
            components = (component_to_factor[x] for x in reaction_target.components_lhs)
            component_factor = Intersection(*components)  # type: VennSet[StateTarget]
            reaction_rules.append(UpdateRule(reaction_target, Intersection(component_factor, contingency_factor).to_simplified_set()))

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
            return Intersection(*(component_to_factor[x] for x in state_target.components))

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

    component_to_factor       = {}  # type: Dict[Spec, VennSet]
    reaction_target_to_factor = {}  # type: Dict[ReactionTarget, VennSet]
    component_state_targets   = []  # type: List[ComponentStateTarget]

    calc_component_factors()
    calc_contingency_factors()
    update_degradations_with_component_states()
    update_degradations_with_contingencies()
    update_degradations_for_interaction_states()
    update_syntheses_with_component_states()

    state_targets    = [StateTarget(x) for x in rxncon_sys.states]
    state_targets.extend(component_state_targets)
    reaction_targets = list(reaction_target_to_factor.keys())

    reaction_rules = []  # type: List[UpdateRule]
    state_rules    = []  # type: List[UpdateRule]

    calc_reaction_rules()
    calc_state_rules()

    return BooleanModel(state_targets, reaction_targets, reaction_rules + state_rules,
                        initial_conditions(reaction_targets, state_targets))


### SIMULATION STUFF ###


def boolnet_from_boolean_model(boolean_model: BooleanModel) -> Tuple[str, Dict[str, str], Dict[str, bool]]:
    """
    Translates the boolean model into the BoolNet syntax.

    Note:
        BoolNet is an R package that provides tools for assembling, analyzing and visualizing Boolean networks.

    Args:
        boolean_model: The boolean model.

    Returns:
        1. The first return value is the boolean model in BooleNet syntax.
        2. The second return value is an abbreviation, target mapping.
        3. The third return value is a mapping of the initial condition.

    """
    def str_from_factor(factor: VennSet) -> str:
        """
        Translates a factor into a string.

        Note:
            During this process the names of the targets are replaced by abbreviations. Reaction targets are replaced
            by R{} and states are replaced by S{} where {} is a continuous numerating for state and reaction targets
            respectively.

        Args:
            factor:

        Returns:
            The string of the factor.

        Raises:
            AssertionError: If the factor is not a valid VennSet object an error is raised.

        """
        if isinstance(factor, ValueSet):
            return boolnet_name_from_target(factor.value)
        elif isinstance(factor, Complement):
            return '!({})'.format(str_from_factor(factor.expr))
        elif isinstance(factor, Intersection):
            return '({})'.format(' & '.join(str_from_factor(x) for x in factor.exprs))
        elif isinstance(factor, Union):
            return '({})'.format(' | '.join(str_from_factor(x) for x in factor.exprs))
        elif isinstance(factor, EmptySet):
            return '0'
        else:
            raise AssertionError('Could not parse factor {}'.format(factor))

    def str_from_update_rule(update_rule: UpdateRule) -> str:
        """
        Creates a string from an update rule.

        Args:
            update_rule: A target and its factor.

        Returns:
            The string of the update rule.

        """
        return '{0}, {1}'.format(boolnet_name_from_target(update_rule.target),
                                 str_from_factor(update_rule.factor))

    def boolnet_name_from_target(target: Target) -> str:
        """
        Creates a valid BoolNet name from the respective target.

        Note:
            The target name is replaced by a abbreviation, which is a valid BoolNet name. Reaction targets are replaced
            by R{} and states are replaced by S{} where {} is a continuous numerating for state and reaction targets
            respectively. The replacement is tracked by boolnet_names.

        Args:
            target: A StateTarget or ReactionTarget.

        Mutates:
            boolnet_names: A mapping of target and its abbreviation.

        Returns:
            The string representation of a valid BoolNet name.

        Raises:
            AssertionError: If the target is neither a ReactionTarget nor a StateTarget an error is raised.

        """
        nonlocal reaction_index
        nonlocal state_index

        try:
            return boolnet_names[target]
        except KeyError:
            if isinstance(target, ReactionTarget):
                name = 'R{}'.format(reaction_index)
                boolnet_names[target] = name
                reaction_index += 1
                return name
            elif isinstance(target, StateTarget):
                name = 'S{}'.format(state_index)
                boolnet_names[target] = name
                state_index += 1
                return name
            else:
                raise AssertionError

    boolnet_names  = {}  # type: Dict[Target, str]

    reaction_index = 0
    state_index    = 0

    def sort_key(rule_str: str) -> Tuple[str, int]:
        """
        Function for sorting.

        Args:
            rule_str: A string representing and update rule of the boolean system.

        Returns:
            A tuple of string and integer.

        """
        target = rule_str.split(',')[0].strip()
        return target[0], int(target[1:])

    rule_strs = sorted([str_from_update_rule(x) for x in boolean_model.update_rules], key=sort_key)

    return 'targets, factors\n' + '\n'.join(rule for rule in rule_strs) + '\n', \
           {name: str(target) for target, name in boolnet_names.items()}, \
           {boolnet_names[target]: value for target, value in boolean_model.initial_conditions.target_to_value.items()}
