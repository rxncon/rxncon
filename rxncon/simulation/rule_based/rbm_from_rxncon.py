from typing import Dict, List, Set

from rxncon.core.rxncon_system import RxnConSystem
from rxncon.core.specification import Specification
from rxncon.core.reaction import Reaction
from rxncon.core.effector import Effector, AndEffector, NotEffector, OrEffector, StateEffector
from rxncon.semantics.molecule_definition import MoleculeDefinition
from rxncon.simulation.rule_based.rule_based_model import RuleBasedModel, Rule
from rxncon.semantics.molecule_definition_from_rxncon import mol_defs_from_rxncon_sys
from rxncon.simulation.rule_based.molecule_from_rxncon import mol_instance_set_from_state_set, mol_instance_set_pair_from_reaction
from rxncon.venntastic.sets import Set, PropertySet, Complement, Union, Intersection

def rbm_from_rxncon_sys(rxconsys: RxnConSystem) -> RuleBasedModel:
    mol_defs = mol_defs_from_rxncon_sys(rxconsys)

    rules = set()

    for reaction in rxconsys.reactions:
        rules = rules.union(_rules_from_reaction(mol_defs, reaction))

    return RuleBasedModel(set(mol_defs.values()), rules, set(), set())


def _rules_from_reaction(mol_defs: Dict[Specification, MoleculeDefinition], reaction: Reaction) -> Set[Rule]:
    pass


def _state_set_from_effector(effector: Effector) -> Set:
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

