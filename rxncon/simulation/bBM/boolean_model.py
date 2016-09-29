from rxncon.venntastic.sets import Set as VennSet, MultiIntersection, MultiUnion, ValueSet, Intersection, Union, Complement, UniversalSet
from rxncon.core.reaction import Reaction
from rxncon.core.state import State
from rxncon.core.spec import MolSpec
from rxncon.core.contingency import Contingency, ContingencyType
from rxncon.core.effector import Effector, AndEffector, OrEffector, NotEffector, StateEffector
from rxncon.core.rxncon_system import RxnConSystem
from typecheck import typecheck
from typing import List


class BooleanModel:
    @typecheck
    def __init__(self, update_rules: List['UpdateRule'], initial_conditions: List['InitialCondition']):
        self.update_rules = update_rules
        self.initial_conditions = initial_conditions
        self._validate()

    def _validate(self):
        all_lhs_targets = []
        all_rhs_targets = []
        for rule in self.update_rules:
            all_lhs_targets.append(rule.target)
            all_rhs_targets += rule.factor_targets

        assert all(x in all_lhs_targets for x in all_rhs_targets)


class InitialCondition:
    @typecheck
    def __init__(self, target: 'Target', value: bool):
        self.target = target
        self.value = value

    @typecheck
    def __eq__(self, other: 'InitialCondition') -> bool:
        return self.target == other.target and self.value == self.value

    @typecheck
    def __repr__(self) -> str:
        return str(self)

    @typecheck
    def __str__(self) -> str:
        return "target: {0}, value: {1}".format(self.target, self.value)


class Target:
    @typecheck
    def __hash__(self) -> int:
        return hash(str(self))

    @typecheck
    def __repr__(self) -> str:
        return str(self)


class ReactionTarget(Target):
    @typecheck
    def __init__(self, reaction_parent: Reaction):
        # maybe need to copy these objects.
        self.reaction_parent     = reaction_parent
        self.produced_targets    = [StateTarget(x) for x in reaction_parent.produced_states]
        self.consumed_targets    = [StateTarget(x) for x in reaction_parent.consumed_states]
        self.synthesised_targets = [StateTarget(x) for x in reaction_parent.synthesised_states]
        self.degraded_targets    = [StateTarget(x) for x in reaction_parent.degraded_states]

    @typecheck
    def __eq__(self, other: 'ReactionTarget'):
        #  Possibly more than one ReactionTarget from a single reaction_parent, so also check all its targets.
        return self.reaction_parent == other.reaction_parent and self.produced_targets == other.produced_targets and \
            self.consumed_targets == other.consumed_targets and self.synthesised_targets == other.synthesised_targets and \
            self.degraded_targets == other.degraded_targets

    @typecheck
    def __str__(self) -> str:
        return str(self.reaction_parent)

    @typecheck
    def produces(self, state_target: StateTarget) -> bool:
        return state_target in self.produced_targets

    @typecheck
    def consumes(self, state_target: StateTarget) -> bool:
        return state_target in self.consumed_targets

    @typecheck
    def synthesises(self, state_target: StateTarget) -> bool:
        return state_target in self.synthesised_targets

    @typecheck
    def degrades(self, state_target: StateTarget) -> bool:
        return state_target in self.degraded_targets

    @typecheck
    @property
    def components(self) -> List[MolSpec]:
        return list(set(self.reaction_parent.components_lhs + self.reaction_parent.components_rhs))


class StateTarget(Target):
    @typecheck
    def __init__(self, state_parent: State):
        self._state_parent = state_parent

    @typecheck
    def __eq__(self, other: 'StateTarget') -> bool:
        return self._state_parent == other._state_parent

    @typecheck
    def is_produced_by(self, reaction_target: ReactionTarget) -> bool:
        return reaction_target.produces(self)

    @typecheck
    def is_consumed_by(self, reaction_target: ReactionTarget) -> bool:
        return reaction_target.consumes(self)

    @typecheck
    def is_synthesised_by(self, reaction_target: ReactionTarget) -> bool:
        return reaction_target.synthesises(self)

    @typecheck
    def is_degraded_by(self, reaction_target: ReactionTarget) -> bool:
        return reaction_target.degrades(self)

    @typecheck
    @property
    def components(self) -> List[MolSpec]:
        return self._state_parent.components

    @typecheck
    @property
    def is_neutral(self) -> bool:
        return self._state_parent.is_neutral

class UpdateRule:
    @typecheck
    def __init__(self, target: Target, factor: VennSet):
        self.target = target
        self.factor = factor
        self._validate()

    def __str__(self):
        return "target: {0}, factors: {1}".format(self.target, self.factor)

    @typecheck
    def factor_targets(self) -> List[Target]:
        return self.factor.values

    def _validate(self):
        assert self.factor.value_type == Target


def boolean_model_from_rxncon(rxncon_sys: RxnConSystem) -> BooleanModel:
    def component_factor(component: MolSpec) -> VennSet:
        grouped_states = rxncon_sys.states_for_component_grouped(component)
        factor = UniversalSet()
        for group in grouped_states:
            factor = Intersection(factor, MultiUnion(*(ValueSet(x) for x in group)))

        return factor

    def contingency_factor(contingency: Contingency) -> VennSet:
        def parse_effector(eff: Effector) -> VennSet:
            if isinstance(eff, StateEffector):
                return ValueSet(eff.expr)
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

    def initial_conditions(reaction_targets: List[ReactionTarget], state_targets: List[StateTarget]) -> List[InitialCondition]:
        conds = [InitialCondition(x, False) for x in reaction_targets]
        conds += [InitialCondition(x, True) for x in state_targets if x.is_neutral]
        conds += [InitialCondition(x, False) for x in state_targets if not x.is_neutral]

        return conds

    reaction_targets = [ReactionTarget(x) for x in rxncon_sys.reactions]
    state_targets    = [StateTarget(x) for x in rxncon_sys.states]

    reaction_rules = []
    state_rules    = []

    # Factor for a reaction target is of the form:
    # components AND contingencies
    for reaction_target in reaction_targets:
        cont_fac = MultiIntersection(*(contingency_factor(x) for x in rxncon_sys.contingencies(reaction_target.reaction_parent)))
        comp_fac = MultiIntersection(*(component_factor(x) for x in reaction_target.components))
        reaction_rules.append(UpdateRule(reaction_target, Intersection(cont_fac, comp_fac)))

    # Factor for a state target is of the form:
    # synthesis OR (components AND NOT degradation AND ((production AND sources) OR (state AND NOT (consumption AND sources))))
    for state_target in state_targets:
        synt_fac = MultiUnion(*(ValueSet(x) for x in reaction_targets if x.synthesises(state_target)))
        comp_fac = MultiIntersection(*(component_factor(x) for x in state_target.components))
        degr_fac = Complement(MultiUnion(*(ValueSet(x) for x in reaction_targets if x.degrades(state_target))))

        prod_cons_facs = []
        for reaction_target in reaction_targets:
            if reaction_target.produces(state_target):
                sources = MultiIntersection(*(ValueSet(x) for x in reaction_target.consumed_targets))
                prod_cons_facs.append(Intersection(ValueSet(reaction_target), sources))

            if reaction_target.consumes(state_target):
                sources = MultiIntersection(*(ValueSet(x) for x in reaction_target.consumed_targets))
                prod_cons_facs.append(Intersection(ValueSet(state_target), Complement(Intersection(ValueSet(reaction_target), sources))))

        prod_cons_fac = MultiUnion(*prod_cons_facs)

        state_rules.append(UpdateRule(state_target, Union(synt_fac, MultiIntersection(comp_fac, degr_fac, prod_cons_fac))))

    return BooleanModel(reaction_rules + state_rules, initial_conditions(reaction_rules, state_rules))

