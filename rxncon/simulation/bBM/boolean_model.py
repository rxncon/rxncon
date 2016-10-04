from rxncon.venntastic.sets import Set as VennSet, MultiIntersection, MultiUnion, ValueSet, Intersection, Union, Complement, UniversalSet
from rxncon.core.reaction import Reaction, matching_reaction_def, reaction_from_string
from rxncon.core.state import State, STATE_DEFS, state_from_string
from rxncon.core.spec import MolSpec
from rxncon.core.contingency import Contingency, ContingencyType
from rxncon.core.effector import Effector, AndEffector, OrEffector, NotEffector, StateEffector
from rxncon.core.rxncon_system import RxnConSystem
from typecheck import typecheck
from typing import List

import re

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
        return "target: {0}, value: {1}".format(str(self.target), str(self.value))


class Target:
#    @typecheck
    def __hash__(self) -> int:
        return hash(str(self))

#    @typecheck
    def __repr__(self) -> str:
        return str(self)


class ReactionTarget(Target):
    @typecheck
    def __init__(self, reaction_parent: Reaction):
        self.reaction_parent     = reaction_parent
        self.produced_targets    = [StateTarget(x) for x in reaction_parent.produced_states]
        self.consumed_targets    = [StateTarget(x) for x in reaction_parent.consumed_states]
        self.synthesised_targets = [StateTarget(x) for x in reaction_parent.synthesised_states] + [ComponentStateTarget(x) for x in reaction_parent.synthesised_components]
        self.degraded_targets    = [StateTarget(x) for x in reaction_parent.degraded_states] + [ComponentStateTarget(x) for x in reaction_parent.degraded_components]

    @typecheck
    def __hash__(self) -> int:
        return hash(str(self))

    @typecheck
    def __eq__(self, other: Target):
        #  Possibly more than one ReactionTarget from a single reaction_parent, so also check all its targets.
        return isinstance(other, ReactionTarget) and self.reaction_parent == other.reaction_parent and self.produced_targets == other.produced_targets and \
            self.consumed_targets == other.consumed_targets and self.synthesised_targets == other.synthesised_targets and \
            self.degraded_targets == other.degraded_targets

    @typecheck
    def __str__(self) -> str:
        return "ReactionTarget<{}>".format(str(self.reaction_parent))

    @typecheck
    def produces(self, state_target: 'StateTarget') -> bool:
        return state_target in self.produced_targets

    @typecheck
    def consumes(self, state_target: 'StateTarget') -> bool:
        return state_target in self.consumed_targets

    @typecheck
    def synthesises(self, state_target: 'StateTarget') -> bool:
        return state_target in self.synthesised_targets

    @typecheck
    def degrades(self, state_target: 'StateTarget') -> bool:
        return state_target in self.degraded_targets

    @property
    @typecheck
    def components(self) -> List[MolSpec]:
        return list(set(self.reaction_parent.components_lhs + self.reaction_parent.components_rhs))


class StateTarget(Target):
    @typecheck
    def __init__(self, state_parent: State):
        self._state_parent = state_parent

    @typecheck
    def __hash__(self) -> int:
        return hash(str(self))

    @typecheck
    def __str__(self) -> str:
        return "StateTarget<{}>".format(str(self._state_parent))

    @typecheck
    def __eq__(self, other: 'Target') -> bool:
        return isinstance(other, StateTarget) and self._state_parent == other._state_parent

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

    @property
    @typecheck
    def components(self) -> List[MolSpec]:
        return self._state_parent.components

    @property
    @typecheck
    def is_neutral(self) -> bool:
        return self._state_parent.is_neutral


class ComponentStateTarget(StateTarget):
    def __init__(self, component: MolSpec):
        self.component = component

    @typecheck
    def __eq__(self, other: Target):
        return isinstance(other, type(self)) and self.component == other.component

    def __str__(self):
        return "ComponentTarget<{}>".format(self.component)

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(str(self))

    @property
    @typecheck
    def components(self) -> List[MolSpec]:
        return [self.component]

    @property
    @typecheck
    def is_neutral(self) -> bool:
        return True


class UpdateRule:
    @typecheck
    def __init__(self, target: Target, factor: VennSet):
        self.target = target
        self.factor = factor
        self._validate()

    def __str__(self):
        return "target: {0}, factors: {1}".format(self.target, self.factor)

    @property
    @typecheck
    def factor_targets(self) -> List[Target]:
        return self.factor.values

    def _validate(self):
        assert any(cls == Target for cls in self.factor.value_type.mro())


def boolean_model_from_rxncon(rxncon_sys: RxnConSystem) -> BooleanModel:
    def naive_component_factor(component: MolSpec) -> VennSet:
        grouped_states = rxncon_sys.states_for_component_grouped(component)
        factor = UniversalSet()
        for group in grouped_states.values():
            factor = Intersection(factor, MultiUnion(*(ValueSet(StateTarget(x)) for x in group)))

        return factor

    @typecheck
    def component_factor(component: MolSpec) -> VennSet:
        if component in stateless_targets.keys():
            return ValueSet(stateless_targets[component])
        else:
            return naive_component_factor(component)

    def contingency_factor(contingency: Contingency) -> VennSet:
        def parse_effector(eff: Effector) -> VennSet:
            if isinstance(eff, StateEffector):
                return ValueSet(StateTarget(eff.expr))
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

    def stateless_component_targets(reaction_targets: List[ReactionTarget]):
        all_components = [component for reaction_target in reaction_targets for component in reaction_target.components]
        stateless_components = [component for component in all_components if naive_component_factor(component) == UniversalSet()]
        return {x: ComponentStateTarget(x) for x in stateless_components}

    reaction_targets  = [ReactionTarget(x) for x in rxncon_sys.reactions]
    stateless_targets = stateless_component_targets(reaction_targets)
    state_targets     = [StateTarget(x) for x in rxncon_sys.states] + list(stateless_targets.values())

    reaction_rules  = []
    state_rules     = []

    # Factor for a reaction target is of the form:
    # components AND contingencies
    for reaction_target in reaction_targets:
        cont_fac = MultiIntersection(*(contingency_factor(x) for x in rxncon_sys.contingencies_for_reaction(reaction_target.reaction_parent)))
        comp_fac = MultiIntersection(*(component_factor(x) for x in reaction_target.components))
        reaction_rules.append(UpdateRule(reaction_target,
                                         Intersection(cont_fac, comp_fac).to_full_simplified_form()))

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
                                      Union(synt_fac, MultiIntersection(comp_fac, degr_fac, prod_cons_fac)).to_full_simplified_form()))

    return BooleanModel(reaction_rules + state_rules, initial_conditions(reaction_targets, state_targets))



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

@typecheck
def boolean_rule_from_str(rule_str: str):

    @typecheck
    def is_state_target(target_str: str):
        if next((x for x in STATE_DEFS if x.matches_repr(target_str)), None) is not None:
            return True
        return False

    @typecheck
    def is_reaction_target(target_str: str):
        if matching_reaction_def(target_str):
            return True
        return False

    @typecheck
    def target_str_to_rule_target(target_str):
        if is_state_target(target_str):
            return StateTarget(state_from_string(target_str))
        elif is_reaction_target(target_str):
            return ReactionTarget(reaction_from_string(target_str))
        else:
            raise NotImplementedError

    # syn | C & ( deg & ( prod & ss1 & ss2 ) | ( s & ! deg & ! con )


    def _get_parentheses_pairs(factor_str):
        parentheses_pairs = []
        for idx, char in enumerate(factor_str):
            if char == "<":
                parentheses_pairs.append([idx])
            if char == ">":
                for parentheses_pair in parentheses_pairs[::-1]:
                    if len(parentheses_pair) == 1:
                        parentheses_pair.append(idx)
                        break
        return [tuple(parentheses_pair) for parentheses_pair in parentheses_pairs]

    def venn_conversion(bool_expr: List[str]):
        return [Complement(ValueSet(target_str_to_rule_target(ele))) if ele.startswith('!') else ValueSet(target_str_to_rule_target(ele)) for ele in bool_expr]

    def get_factor(factor_str: str):
        if factor_str.split("&"):
            return MultiIntersection(*venn_conversion(factor_str.split("&")))
        elif factor_str.split("|"):
            return MultiUnion(*venn_conversion(factor_str.split("|")))
        else:
            raise NotImplementedError


    #def factor_str_to_rule_factor(factor_str: str):
        #1) finde die klammer paare
        #3) löse klammern auf die keine weiteren klammern enthalten
        #2) finde benachbarte klammern
    #    factor_str = factor_str.replace(" ", "")
    #    parentheses_pairs = _get_parentheses_pairs(factor_str)
    #    parentheses_set_mapping = {}

        # for parentheses_pair in parentheses_pairs:
        #     sub_factor_str = factor_str[parentheses_pair[0]+1:parentheses_pair[1]]
        #     if not "<" in sub_factor_str or not ">" in sub_factor_str:
        #         parentheses_set_mapping[parentheses_pair] = get_factor(sub_factor_str)
        #
        # def merge_perenthesis():
        #     pass
        # # todo This can be done recursively we have always tuple
        # # merge directly neighboring tuple
        # # if there are no directly neighboring tuple any more check if the start and end are directly continues
        # # if check if the last or first parenthesis is directly neighbored
        # # <<a|b>&<a|c>>
        # # 0        6           16 17         27  29            43 44 45 46
        # # < a--a | <b_ppi+_a | <  < c--a&d--e > | <e_ppi-_a&f--r>   >  > >
        #
        # parentheses_pairs
    #'f--e=<a--a | <b_ppi+_a | <<c--a&d--e> | <e_ppi-_a&f--r> >>>'
    def factor_str_to_rule_factor(factor_str: str):
        if factor_str.startswith('<') or factor_str.startswith('>'):
            return factor_str_to_rule_factor(factor_str[1:])
        elif re.match('[\w\]\[\}\{()-_]',factor_str):
            rule_target_match_AND_inner_OR = re.match('([\w\[\]\{\}()-_]+[^>])[&]([\w\[\]\{\}()-_]+)>[|]', factor_str)
            rule_target_match_AND_inner_AND = re.match('([\w\[\]\{\}()-_]+[^>])[&]([\w\[\]\{\}()-_]+)>[&]', factor_str)
            rule_target_match_OR_inner_OR = re.match('([\w\[\]\{\}()-_]+[^>])[|]([\w\[\]\{\}()-_]+)>[|]', factor_str)
            rule_target_match_OR_inner_AND = re.match('([\w\[\]\{\}()-_]+[^>])[|]([\w\[\]\{\}()-_]+)>[&]', factor_str)
            rule_target_match_AND = re.match('([\w\]\[\}\{()-_]+[^>])[&]', factor_str)
            rule_target_match_OR  = re.match('([\w\]\[\}\{()-_]+[^>])[|]', factor_str)
            rule_target_match = re.match('[\w\]\[\}\{()-_]+[^>|&]', factor_str)
            if rule_target_match_AND_inner_OR:
                return Union(Intersection(ValueSet(target_str_to_rule_target(rule_target_match_AND_inner_OR.group(1))), ValueSet(target_str_to_rule_target(rule_target_match_AND_inner_OR.group(2)))),
                             factor_str_to_rule_factor(factor_str[rule_target_match_AND_inner_OR.end():]))
            elif rule_target_match_AND_inner_AND:
                return Union(Intersection(ValueSet(target_str_to_rule_target(rule_target_match_AND_inner_AND.group(1))),
                                          ValueSet(target_str_to_rule_target(rule_target_match_AND_inner_AND.group(2)))),
                             factor_str_to_rule_factor(factor_str[rule_target_match_AND_inner_AND.end():]))
            elif rule_target_match_OR_inner_OR:
                return Union(Intersection(ValueSet(target_str_to_rule_target(rule_target_match_OR_inner_OR.group(1))),
                                          ValueSet(target_str_to_rule_target(rule_target_match_OR_inner_OR.group(2)))),
                             factor_str_to_rule_factor(factor_str[rule_target_match_OR_inner_OR.end():]))
            elif rule_target_match_OR_inner_AND:
                return Union(Intersection(ValueSet(target_str_to_rule_target(rule_target_match_OR_inner_AND.group(1))),
                                          ValueSet(target_str_to_rule_target(rule_target_match_OR_inner_AND.group(2)))),
                             factor_str_to_rule_factor(factor_str[rule_target_match_OR_inner_AND.end():]))
            elif rule_target_match_AND:
                return Intersection(ValueSet(target_str_to_rule_target(rule_target_match_AND.group(1))),factor_str_to_rule_factor(factor_str[rule_target_match_AND.end():]))
            elif rule_target_match_OR:
                return Union(ValueSet(target_str_to_rule_target(rule_target_match_OR.group(1))), factor_str_to_rule_factor(factor_str[rule_target_match_OR.end():]))
            elif rule_target_match:
                return ValueSet(target_str_to_rule_target(rule_target_match.group()))

    rule_str = rule_str.replace(" ", "")
    target_str, factor_str = rule_str.split('=')
    assert factor_str.count("<") == factor_str.count('>')
    target = target_str_to_rule_target(target_str)
    factor = factor_str_to_rule_factor(factor_str)

    return UpdateRule(target, factor)
