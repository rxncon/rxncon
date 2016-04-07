from rxncon.core.rxncon_system import RxnConSystem
from rxncon.simulation.rule_based.rule_based_model import RuleBasedModel, Rule
from rxncon.semantics.molecule_definition_from_rxncon import MoleculeDefinitionSupervisor


def rbm_from_rxncon(rxconsys: RxnConSystem) -> RuleBasedModel:
    mol_defs = MoleculeDefinitionSupervisor(rxconsys).molecule_definitions