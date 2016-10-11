from rxncon.venntastic.sets import Set as VennSet, MultiIntersection, MultiUnion, ValueSet, Intersection, Union, Complement, UniversalSet, venn_from_pyeda
from rxncon.core.reaction import Reaction
from rxncon.core.state import State
from rxncon.core.spec import Spec
from rxncon.core.contingency import Contingency, ContingencyType
from rxncon.core.effector import Effector, AndEffector, OrEffector, NotEffector, StateEffector
from rxncon.core.rxncon_system import RxnConSystem

from typing import List, Dict


class BooleanModel:
    def __init__(self, update_rules: List['UpdateRule'], initial_conditions: 'BooleanModelConfig'):
        self.update_rules = update_rules
        self.initial_conditions = initial_conditions
        self._validate_update_rules()
        self._validate_initial_conditions()

    def set_initial_condition(self, target: Target, value: bool):
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
    def __init__(self, target_to_value: Dict[Target, bool]):
        self.target_to_value = target_to_value

    def set_target(self, target: Target, value: bool):
        self.target_to_value[target] = value

    def validate_by_model(self, model: BooleanModel):
        model_targets  = [rule.target for rule in model.update_rules]
        config_targets = self.target_to_value.keys()

        assert set(model_targets) == set(config_targets)


class Target:
    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)


class ReactionTarget(Target):
    def __init__(self, reaction_parent: Reaction):
        # maybe need to copy these objects.
        self.reaction_parent     = reaction_parent
        self.produced_targets    = [StateTarget(x) for x in reaction_parent.produced_states]
        self.consumed_targets    = [StateTarget(x) for x in reaction_parent.consumed_states]
        self.synthesised_targets = [StateTarget(x) for x in reaction_parent.synthesised_states] + [ComponentStateTarget(x) for x in reaction_parent.synthesised_components]
        self.degraded_targets    = [StateTarget(x) for x in reaction_parent.degraded_states] + [ComponentStateTarget(x) for x in reaction_parent.degraded_components]

    def __hash__(self) -> int:
        return hash(str(self))

    def __eq__(self, other: Target):
        #  Possibly more than one ReactionTarget from a single reaction_parent, so also check all its targets.
        return isinstance(other, ReactionTarget) and self.reaction_parent == other.reaction_parent and self.produced_targets == other.produced_targets and \
            self.consumed_targets == other.consumed_targets and self.synthesised_targets == other.synthesised_targets and \
            self.degraded_targets == other.degraded_targets

    def __str__(self) -> str:
        return str(self.reaction_parent)

    def produces(self, state_target: 'StateTarget') -> bool:
        return state_target in self.produced_targets

    def consumes(self, state_target: 'StateTarget') -> bool:
        return state_target in self.consumed_targets

    def synthesises(self, state_target: 'StateTarget') -> bool:
        return state_target in self.synthesised_targets

    def degrades(self, state_target: 'StateTarget') -> bool:
        return state_target in self.degraded_targets

    @property
    def components(self) -> List[Spec]:
        return list(set(self.reaction_parent.components_lhs + self.reaction_parent.components_rhs))


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
        self._validate()

    def __str__(self):
        return "target: {0}, factors: {1}".format(self.target, self.factor)

    @property
    def factor_targets(self) -> List[Target]:
        return self.factor.values

    def _validate(self):
        pass
        # assert any(cls == Target for cls in self.factor.value_type.mro())


def boolean_model_from_rxncon(rxncon_sys: RxnConSystem) -> BooleanModel:
    def naive_component_factor(component: Spec) -> VennSet:
        grouped_states = rxncon_sys.states_for_component_grouped(component)
        if not grouped_states.values():
            return UniversalSet()

        return MultiIntersection(*(MultiUnion(*(ValueSet(StateTarget(x)) for x in group)) for group in grouped_states.values()))

    def component_factor(component: Spec) -> VennSet:
        if component in stateless_targets.keys():
            return ValueSet(stateless_targets[component])
        else:
            return naive_component_factor(component)

    def contingency_factor(contingency: Contingency) -> VennSet:
        def parse_effector(eff: Effector) -> VennSet:
            if isinstance(eff, StateEffector):
                return ValueSet(StateTarget(eff.expr.to_non_struct_state()))
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

    def stateless_component_targets(reaction_targets: List[ReactionTarget]):
        all_components = [component for reaction_target in reaction_targets for component in reaction_target.components]
        stateless_components = [component for component in all_components if naive_component_factor(component) == UniversalSet()]
        return {x: ComponentStateTarget(x) for x in stateless_components}

    reaction_targets  = [ReactionTarget(x) for x in rxncon_sys.reactions]
    stateless_targets = stateless_component_targets(reaction_targets)
    state_targets     = [StateTarget(x.to_non_struct_state()) for x in rxncon_sys.states] + list(stateless_targets.values())

    reaction_rules  = []
    state_rules     = []

    # Factor for a reaction target is of the form:
    # components AND contingencies
    for reaction_target in reaction_targets:
        cont_fac = MultiIntersection(*(contingency_factor(x) for x in rxncon_sys.contingencies_for_reaction(reaction_target.reaction_parent)))
        comp_fac = MultiIntersection(*(component_factor(x) for x in reaction_target.components))
        reaction_rules.append(UpdateRule(reaction_target, Intersection(cont_fac, comp_fac).to_simplified_set()))

    # Factor for a state target is of the form:
    # synthesis OR (components AND NOT degradation AND ((production AND sources) OR (state AND NOT (consumption AND sources))))
    for state_target in state_targets:
        synt_fac = MultiUnion(*(ValueSet(x) for x in reaction_targets if x.synthesises(state_target)))
        comp_fac = MultiIntersection(*(component_factor(x) for x in state_target.components))
        degr_fac = Complement(MultiUnion(*(ValueSet(x) for x in reaction_targets if x.degrades(state_target))))

        prod_facs = []
        cons_facs = []

        for reaction_target in reaction_targets:
            if reaction_target.produces(state_target):
                sources = MultiIntersection(*(ValueSet(x) for x in reaction_target.consumed_targets))
                prod_facs.append(Intersection(ValueSet(reaction_target), sources))

            if reaction_target.consumes(state_target):
                sources = MultiIntersection(*(ValueSet(x) for x in reaction_target.consumed_targets))
                cons_facs.append(Complement(Intersection(ValueSet(reaction_target), sources)))

        prod_cons_fac = Union(MultiUnion(*prod_facs), Intersection(ValueSet(state_target), MultiIntersection(*cons_facs)))

        state_rules.append(UpdateRule(state_target,
                                      Union(synt_fac, MultiIntersection(comp_fac, degr_fac, prod_cons_fac)).to_simplified_set()))

    return BooleanModel(reaction_rules + state_rules, initial_conditions(reaction_targets, state_targets))


### SIMULATION STUFF ###


def boolnet_str_from_boolean_model(boolean_model: BooleanModel) -> str:
    def clean_str(the_str: str) -> str:
        dirty_chars = {
            '[': 'DOM',
            ']': 'DOM',
            '/': '',
            '@': '',
            '{': 'MOD',
            '}': 'MOD',
            ' ': '',
            '(': 'RES',
            ')': 'RES',
            '<': '',
            '>': ''
        }
        for dirty_char, clean_char in dirty_chars.items():
            the_str = the_str.replace(dirty_char, clean_char)

        return the_str

    def str_from_factor(factor: VennSet) -> str:
        if isinstance(factor, ValueSet):
            return clean_str(str(factor.value))
        elif isinstance(factor, Complement):
            return '!({})'.format(str_from_factor(factor.expr))
        elif isinstance(factor, Intersection):
            return '({0} & {1})'.format(str_from_factor(factor.left_expr), str_from_factor(factor.right_expr))
        elif isinstance(factor, Union):
            return '({0} | {1})'.format(str_from_factor(factor.left_expr), str_from_factor(factor.right_expr))
        else:
            raise AssertionError

    def str_from_update_rule(update_rule: UpdateRule) -> str:
        return '{0} , {1}'.format(clean_str(str(update_rule.target)),
                                  str_from_factor(update_rule.factor))

    return 'targets, factors\n' + '\n'.join(str_from_update_rule(x) for x in boolean_model.update_rules) + '\n'


class BooleanModelConfigPath:
    pass


class BooleanSimulator:
    def __init__(self, boolean_model: BooleanModel):
        self.boolean_model = boolean_model

    @property
    def update_rules(self):
        return self.boolean_model.update_rules

    @property
    def initial_conditions(self):
        return self.boolean_model.initial_conditions

    def set_initian_condition(self, target: Target, value: bool):
        self.boolean_model.set_initial_condition(target, value)

    def calc_attractor_path(self) -> BooleanModelConfigPath:
        return BooleanModelConfigPath()



