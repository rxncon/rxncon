from typing import Dict, List

from rxncon.core.rxncon_system import RxnConSystem
from rxncon.core.specification import Specification
from rxncon.core.reaction import Reaction
from rxncon.semantics.molecule_definition import MoleculeDefinition
from rxncon.simulation.rule_based.rule_based_model import RuleBasedModel, Rule
from rxncon.semantics.molecule_definition_from_rxncon import mol_defs_from_rxncon_sys
from rxncon.simulation.rule_based.molecule_from_rxncon import mol_instance_set_from_state_set, mol_instance_set_pair_from_reaction

def rbm_from_rxncon_sys(rxconsys: RxnConSystem) -> RuleBasedModel:
    mol_defs = mol_defs_from_rxncon_sys(rxconsys)

    rules = []

    for reaction in rxconsys.reactions:
        rules += _rules_from_reaction(mol_defs, reaction)

    return RuleBasedModel(sorted(mol_defs.values()), rules, [], [])


def _rules_from_reaction(mol_defs: Dict[Specification, MoleculeDefinition], reaction: Reaction) -> List[Rule]:
    pass
