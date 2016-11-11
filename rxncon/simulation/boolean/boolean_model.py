from typing import List, Dict, Tuple, Union
from copy import deepcopy
from enum import Enum

from rxncon.venntastic.sets import Set as VennSet, ValueSet, Intersection, Union, Complement, UniversalSet
from rxncon.core.reaction import Reaction
from rxncon.core.state import State
from rxncon.core.spec import Spec
from rxncon.core.contingency import Contingency, ContingencyType
from rxncon.core.effector import Effector, AndEffector, OrEffector, NotEffector, StateEffector
from rxncon.core.rxncon_system import RxnConSystem



class BooleanModel:
    def __init__(self, update_rules: List['UpdateRule'], initial_conditions: 'BooleanModelConfig'):
        self.update_rules = update_rules
        self.initial_conditions = initial_conditions
        self._validate_update_rules()
        self._validate_initial_conditions()

    def set_initial_condition(self, target: 'Target', value: bool):
        self.initial_conditions.set_target(target, value)

    def _validate_update_rules(self):
        all_lhs_targets = []
        all_rhs_targets = []
        for rule in self.update_rules:
            all_lhs_targets.append(rule.target)
            all_rhs_targets += rule.factor_targets

        assert all(x in all_lhs_targets for x in all_rhs_targets)

    def _validate_initial_conditions(self):
        self.initial_conditions.validate_by_model(self)


class BooleanModelConfig:
    def __init__(self, target_to_value: Dict['Target', bool]):
        self.target_to_value = target_to_value

    def set_target(self, target: 'Target', value: bool):
        self.target_to_value[target] = value

    def validate_by_model(self, model: BooleanModel):
        model_targets  = [rule.target for rule in model.update_rules]
        config_targets = self.target_to_value.keys()

        assert set(model_targets) == set(config_targets) and len(model_targets) == len(config_targets)


class Target:
    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)


class ReactionTarget(Target):
    def __init__(self, reaction_parent: Reaction):
        self.reaction_parent     = reaction_parent
        self.produced_targets    = [StateTarget(x) for x in reaction_parent.produced_states]
        self.consumed_targets    = [StateTarget(x) for x in reaction_parent.consumed_states]
        self.synthesised_targets = [StateTarget(x) for x in reaction_parent.synthesised_states]
        self.degraded_targets    = [StateTarget(x) for x in reaction_parent.degraded_states]

        self.contingency_variant_index = 0
        self.interaction_variant_index = 0

    def __hash__(self) -> int:
        return hash(str(self))

    def __eq__(self, other: Target):
        return isinstance(other, ReactionTarget) and self.reaction_parent == other.reaction_parent and \
            self.contingency_variant_index == other.contingency_variant_index and \
            self.interaction_variant_index == other.interaction_variant_index

    def __str__(self) -> str:
        suffix = ''
        if self.interaction_variant_index != 0:
            suffix = '#{},{}'.format(self.contingency_variant_index, self.interaction_variant_index)
        elif self.contingency_variant_index != 0 and self.interaction_variant_index == 0:
            suffix = '#{}'.format(self.contingency_variant_index)

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


class StateTarget(Target):
    def __init__(self, state_parent: State):
        self._state_parent = state_parent

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return str(self._state_parent)

    def __eq__(self, other: 'Target') -> bool:
        return isinstance(other, StateTarget) and self._state_parent == other._state_parent

    def is_produced_by(self, reaction_target: ReactionTarget) -> bool:
        return reaction_target.produces(self)

    def is_consumed_by(self, reaction_target: ReactionTarget) -> bool:
        return reaction_target.consumes(self)

    def is_synthesised_by(self, reaction_target: ReactionTarget) -> bool:
        return reaction_target.synthesises(self)

    def is_degraded_by(self, reaction_target: ReactionTarget) -> bool:
        return reaction_target.degrades(self)

    @property
    def components(self) -> List[Spec]:
        return self._state_parent.components

    @property
    def is_neutral(self) -> bool:
        return self._state_parent.is_neutral

    @property
    def neutral_targets(self) -> List['StateTarget']:
        return [StateTarget(x) for x in self._state_parent.neutral_states]


class ComponentStateTarget(StateTarget):
    def __init__(self, component: Spec):
        self.component = component

    def __eq__(self, other: Target):
        return isinstance(other, type(self)) and self.component == other.component

    def __str__(self):
        return str(self.component)

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(str(self))

    @property
    def components(self) -> List[Spec]:
        return [self.component]

    @property
    def is_neutral(self) -> bool:
        return True


class UpdateRule:
    def __init__(self, target: Target, factor: VennSet):
        self.target = target
        self.factor = factor

    def __str__(self):
        return "target: {0}, factors: {1}".format(self.target, self.factor)

    @property
    def factor_targets(self) -> List[Target]:
        return self.factor.values


class SmoothingStrategy(Enum):
    no_smoothing              = 'no_smoothing'
    smooth_production_sources = 'smooth_production_sources'


def boolean_model_from_rxncon(rxncon_sys: RxnConSystem,
                              smoothing_strategy: SmoothingStrategy=SmoothingStrategy.no_smoothing) -> BooleanModel:
    def factor_from_contingency(contingency: Contingency) -> VennSet:
        def parse_effector(eff: Effector) -> VennSet:
            if isinstance(eff, StateEffector):
                return ValueSet(StateTarget(eff.expr.to_non_structured_state()))
            elif isinstance(eff, NotEffector):
                return Complement(parse_effector(eff.expr))
            elif isinstance(eff, OrEffector):
                return Union(parse_effector(eff.left_expr), parse_effector(eff.right_expr))
            elif isinstance(eff, AndEffector):
                return Intersection(parse_effector(eff.left_expr), parse_effector(eff.right_expr))
            else:
                raise AssertionError

        if contingency.type in [ContingencyType.requirement, ContingencyType.positive]:
            return parse_effector(contingency.effector)
        elif contingency.type in [ContingencyType.inhibition, ContingencyType.negative]:
            return Complement(parse_effector(contingency.effector))
        else:
            return UniversalSet()

    def initial_conditions(reaction_targets: List[ReactionTarget], state_targets: List[StateTarget]) -> BooleanModelConfig:
        conds = {}

        for target in reaction_targets:
            conds[target] = False

        for target in state_targets:
            if target.is_neutral:
                conds[target] = True
            else:
                conds[target] = False

        return BooleanModelConfig(conds)

    def calc_component_factors():
        for component in rxncon_sys.components():
            grouped_states = rxncon_sys.states_for_component_grouped(component)
            if not grouped_states.values():
                component_state_targets.append(ComponentStateTarget(component))
                component_to_factor[component] = ValueSet(ComponentStateTarget(component))
            else:
                component_to_factor[component] = \
                    Intersection(*(Union(*(ValueSet(StateTarget(x)) for x in group)) for group in grouped_states.values()))

    def calc_contingency_factors():
        for reaction in rxncon_sys.reactions:
            cont = Intersection(*(factor_from_contingency(x) for x in rxncon_sys.contingencies_for_reaction(reaction))).to_simplified_set()

            if not reaction.degraded_components:
                reaction_target_to_factor[ReactionTarget(reaction)] = cont
            else:
                for index, factor in enumerate(cont.to_dnf_list()):
                    target = ReactionTarget(reaction)
                    target.contingency_variant_index = index
                    reaction_target_to_factor[target] = factor

    def update_degradations_with_contingencies():
        for reaction_target, contingency_factor in reaction_target_to_factor.items():
            # First case: reaction with a non-trivial contingency should degrade only the states appearing
            # in the contingency that are connected to the degraded component.
            if reaction_target.degraded_components and not contingency_factor.is_equivalent_to(UniversalSet()):
                for degraded_component in reaction_target.degraded_components:
                    if ComponentStateTarget(degraded_component) in component_state_targets:
                        reaction_target.degraded_targets.append(ComponentStateTarget(degraded_component))
                    else:
                        reaction_target.degraded_targets += [state_target for state_target in contingency_factor.values
                                                             if degraded_component in state_target.components]

            # Second case: reaction with a trivial contingency should degrade all states for the degraded component.
            elif reaction_target.degraded_components and contingency_factor.is_equivalent_to(UniversalSet()):
                for degraded_component in reaction_target.degraded_components:
                    if ComponentStateTarget(degraded_component) in component_state_targets:
                        reaction_target.degraded_targets.append(ComponentStateTarget(degraded_component))
                    else:
                        reaction_target.degraded_targets += \
                            [StateTarget(x) for x in rxncon_sys.states_for_component(degraded_component)]

    def update_degradations_for_interaction_states():
        new_reactions = {}

        for reaction_target, contingency_factor in reaction_target_to_factor.items():
            for num, interaction_state in enumerate(state for state in rxncon_sys.states if len(state.components) > 1
                                                    and reaction_target.degrades(StateTarget(state))):
                neutral_targets = StateTarget(interaction_state).neutral_targets
                new_reaction_target = deepcopy(reaction_target)
                new_reaction_target.produced_targets += \
                    [x for x in neutral_targets if not any(component in new_reaction_target.degraded_components for component in x.components)]
                new_reaction_target.interaction_variant_index = num + 1

                new_reactions[new_reaction_target] = Intersection(contingency_factor, ValueSet(StateTarget(interaction_state)))

        reaction_target_to_factor.update(new_reactions)

    def update_syntheses_with_component_states():
        for reaction_target, _ in reaction_target_to_factor.items():
            for component in reaction_target.synthesised_components:
                if ComponentStateTarget(component) in component_state_targets:
                    reaction_target.synthesised_targets.append(ComponentStateTarget(component))

    def calc_reaction_rules():
        # Factor for a reaction target is of the form:
        # components AND contingencies
        for reaction_target, contingency_factor in reaction_target_to_factor.items():
            component_factor = Intersection(*(component_to_factor[x] for x in reaction_target.components_lhs))
            reaction_rules.append(UpdateRule(reaction_target, Intersection(component_factor, contingency_factor).to_simplified_set()))

    def calc_state_rules():
        def reaction_with_sources(reaction_target: ReactionTarget) -> VennSet:
            return Intersection(ValueSet(reaction_target),
                                Intersection(*(ValueSet(x) for x in reaction_target.consumed_targets)))

        # Factor for a state target is of the form:
        # synthesis OR (components AND NOT degradation AND ((production AND sources) OR (state AND NOT (consumption AND sources))))
        for state_target in state_targets:
            synt_fac = Union(*(ValueSet(x) for x in reaction_targets if x.synthesises(state_target)))
            comp_fac = Intersection(*(component_to_factor[x] for x in state_target.components))
            degr_fac = Complement(Union(*(ValueSet(x) for x in reaction_targets if x.degrades(state_target))))

            prod_facs = []
            cons_facs = []

            for reaction_target in (target for target in reaction_targets if target.produces(state_target)):
                if smoothing_strategy == SmoothingStrategy.no_smoothing:
                    prod_facs.append(reaction_with_sources(reaction_target))
                elif smoothing_strategy == SmoothingStrategy.smooth_production_sources:
                    smoothed_prod_facs = []
                    for primary_source in reaction_target.consumed_targets:
                        smoothed_prod_facs.append(Union(ValueSet(primary_source),
                                                        *(reaction_with_sources(rxn) for rxn in reaction_targets if rxn.produces(primary_source))))
                    prod_facs.append(Intersection(ValueSet(reaction_target), *smoothed_prod_facs))
                else:
                    raise AssertionError

            for reaction_target in (target for target in reaction_targets if target.consumes(state_target)):
                cons_facs.append(Complement(reaction_with_sources(reaction_target)))

            prod_cons_fac = Union(Union(*prod_facs), Intersection(ValueSet(state_target), Intersection(*cons_facs)))

            state_rules.append(UpdateRule(state_target,
                                          Union(synt_fac,
                                                Intersection(comp_fac, degr_fac, prod_cons_fac)).to_simplified_set()))

    component_to_factor       = {}  # type: Dict[Spec, VennSet]
    reaction_target_to_factor = {}  # type: Dict[ReactionTarget, VennSet]
    component_state_targets   = []  # type: List[ComponentStateTarget]

    calc_component_factors()
    calc_contingency_factors()
    update_degradations_with_contingencies()
    update_degradations_for_interaction_states()
    update_syntheses_with_component_states()

    state_targets    = component_state_targets + [StateTarget(x) for x in rxncon_sys.states]
    reaction_targets = list(reaction_target_to_factor.keys())

    reaction_rules = []  # type: List[UpdateRule]
    state_rules    = []  # type: List[UpdateRule]

    calc_reaction_rules()
    calc_state_rules()

    return BooleanModel(reaction_rules + state_rules, initial_conditions(reaction_targets, state_targets))


### SIMULATION STUFF ###


def boolnet_from_boolean_model(boolean_model: BooleanModel) -> Tuple[str, Dict[str, str], Dict[str, bool]]:
    def str_from_factor(factor: VennSet) -> str:
        if isinstance(factor, ValueSet):
            return boolnet_name_from_target(factor.value)
        elif isinstance(factor, Complement):
            return '!({})'.format(str_from_factor(factor.expr))
        elif isinstance(factor, Intersection):
            return '({})'.format(' & '.join(str_from_factor(x) for x in factor.exprs))
        elif isinstance(factor, Union):
            return '({})'.format(' | '.join(str_from_factor(x) for x in factor.exprs))
        else:
            raise AssertionError

    def str_from_update_rule(update_rule: UpdateRule) -> str:
        return '{0}, {1}'.format(boolnet_name_from_target(update_rule.target),
                                 str_from_factor(update_rule.factor))

    def boolnet_name_from_target(target: Target) -> str:
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
                return AssertionError

    # boolnet_name_from_target closes over these variables.
    boolnet_names  = {}
    reaction_index = 0
    state_index    = 0

    def sort_key(rule_str):
        target = rule_str.split(',')[0].strip()
        return target[0], int(target[1:])

    rule_strs = sorted([str_from_update_rule(x) for x in boolean_model.update_rules], key=sort_key)

    return 'targets, factors\n' + '\n'.join(rule for rule in rule_strs) + '\n', \
           {name: str(target) for target, name in boolnet_names.items()}, \
           {boolnet_names[target]: value for target, value in boolean_model.initial_conditions.target_to_value.items()}
