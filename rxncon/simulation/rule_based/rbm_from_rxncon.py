from typing import Dict, List, Set, Iterable, Tuple
import itertools as itt

from rxncon.core.rxncon_system import RxnConSystem
from rxncon.core.specification import Specification
from rxncon.core.reaction import Reaction
from rxncon.core.state import State
from rxncon.core.effector import Effector, AndEffector, NotEffector, OrEffector, StateEffector
from rxncon.semantics.molecule_definition import MoleculeDefinition
from rxncon.simulation.rule_based.rule_based_model import RuleBasedModel, Rule
from rxncon.semantics.molecule_definition_from_rxncon import mol_defs_from_rxncon_sys
from rxncon.simulation.rule_based.molecule_from_rxncon import mol_instance_set_from_state_set, mol_instance_set_pair_from_reaction
from rxncon.venntastic.sets import Set as VennSet, PropertySet, Complement, Union, Intersection, nested_expression_from_list_and_binary_op


def rbm_from_rxncon_sys(rxconsys: RxnConSystem) -> RuleBasedModel:
    mol_defs = mol_defs_from_rxncon_sys(rxconsys)

    rules = set()

    for reaction in rxconsys.reactions:
        rules = rules.union(_rules_from_reaction(mol_defs, reaction))

    return RuleBasedModel(set(mol_defs.values()), rules, set(), set())


def _rules_from_reaction(mol_defs: Dict[Specification, MoleculeDefinition], reaction: Reaction) -> Set[Rule]:
    pass


def _state_set_from_effector(effector: Effector) -> VennSet:
    if isinstance(effector, StateEffector):
        return PropertySet(effector.expr)
    if isinstance(effector, NotEffector):
        return Complement(_state_set_from_effector(effector.expr))
    elif isinstance(effector, AndEffector):
        return Intersection(_state_set_from_effector(effector.left_expr), _state_set_from_effector(effector.right_expr))
    elif isinstance(effector, OrEffector):
        return Union(_state_set_from_effector(effector.left_expr), _state_set_from_effector(effector.right_expr))
    else:
        raise AssertionError


def _quant_contingency_configs_from_reaction(rxnconsys: RxnConSystem, reaction: Reaction) -> Set['_QuantitativeContingencyConfiguration']:
    states = set()
    for contingency in rxnconsys.quantitative_contingencies_for_reaction(reaction):
        states = states.union(set(contingency.effector.states))

    combis = _true_false_combinations(states)

    return {_QuantitativeContingencyConfiguration(combi[0], combi[1]) for combi in combis}


def _true_false_combinations(xs: Set) -> List[Tuple[Set, Set]]:
    result = []
    for n in range(len(xs)):
        true_combis = itt.combinations(xs, n)
        for trues in true_combis:
            result.append((set(trues), {f for f in xs if f not in trues}))

    return result


class _QuantitativeContingencyConfiguration:
    def __init__(self, present_states: Set[State], absent_states: Set[State]):
        self.present_states = present_states
        self.absent_states = absent_states

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return '!{0}'.join(sorted(self.present_states)) + 'x{0}'.join(sorted(self.absent_states))

    def to_state_set(self) -> VennSet:
        return Intersection(
            nested_expression_from_list_and_binary_op([PropertySet(x) for x in self.present_states], Intersection),
            nested_expression_from_list_and_binary_op([Complement(PropertySet(x)) for x in self.absent_states], Intersection)
        )













