from typing import List
from collections import defaultdict
import typecheck as tc

import rxncon.core.rxncon_system as rxs
import rxncon.core.reaction as rxn
import rxncon.simulation.rule_based.rule_based_model as rbm
import rxncon.core.state as sta
import rxncon.core.contingency as con


@tc.typecheck
def rule_based_model_from_rxncon(rxnconsys: rxs.RxnConSystem) -> rbm.RuleBasedModel:
    mol_defs = molecule_defs_from_rxncon(rxnconsys)

    rules = []
    for reaction in rxnconsys.reactions:
        lhs_base_reactants = base_reactants(mol_defs,
                                            rxnconsys.explicit_contingencies_for_reaction(reaction) + rxnconsys.implicit_contingencies_for_reaction(reaction),
                                            [])

        rhs_base_reactants = base_reactants(mol_defs,
                                            rxnconsys.explicit_contingencies_for_reaction(reaction),
                                            rxnconsys.implicit_contingencies_for_reaction(reaction))

        base_rules = [] # rbm.Rule(lhs_base_reactants, rhs_base_reactants, arrow_type_from_reaction(reaction), base_rates_from_reaction(reaction))
        for base_rule in base_rules:
            rules.extend(derived_rules_from_base_rule_and_contingencies(base_rule, rxnconsys.explicit_contingencies_for_reaction(reaction)))

    non_overlapping_rules = non_overlapping_rules_from_rules(rules)

    return rbm.RuleBasedModel(mol_defs, non_overlapping_rules, None, None)


@tc.typecheck
def molecule_defs_from_rxncon(rxnconsys: rxs.RxnConSystem) -> List[rbm.MoleculeDefinition]:
    unmodified = 'u'
    assoc_prefix = 'Assoc'
    mod_prefix = 'Mod'

    names = set()
    name_to_assoc_defs = defaultdict(set)
    name_to_mod_defs = defaultdict(set)
    name_to_valid_compartments = defaultdict(set)

    for reaction in rxnconsys.reactions:
        for state in [x for x in (reaction.source, reaction.product) if x]:
            if isinstance(state, sta.CovalentModificationState):
                names.add(state.substrate.name)
                modifiers = [unmodified, state.modifier.value]
                domain = state.substrate.domain if state.substrate.domain else '{0}{1}'.format(mod_prefix, reaction.subject.name)

                name_to_mod_defs[state.substrate.name].add(rbm.ModificationDefinition(domain, modifiers))

            elif isinstance(state, sta.InterProteinInteractionState) or isinstance(state, sta.IntraProteinInteractionState):
                names.add(state.first_component.name)
                names.add(state.second_component.name)

                first_domain = state.first_component.domain if state.first_component.domain \
                    else '{0}{1}'.format(assoc_prefix, state.second_component.name)

                second_domain = state.second_component.domain if state.second_component.domain \
                    else '{0}{1}'.format(assoc_prefix, state.first_component.name)

                name_to_assoc_defs[state.first_component.name].add(rbm.AssociationDefinition(first_domain))
                name_to_assoc_defs[state.second_component.name].add(rbm.AssociationDefinition(second_domain))

            elif isinstance(state, sta.SynthesisDegradationState):
                names.add(state.component.name)

            elif isinstance(state, sta.TranslocationState):
                names.add(state.substrate.name)
                name_to_valid_compartments[state.substrate.name].add(state.compartment)

    mol_defs = []

    for name in names:
        assoc_defs = list(name_to_assoc_defs[name]) if name in name_to_assoc_defs else None
        mod_defs = list(name_to_mod_defs[name]) if name in name_to_mod_defs else None
        loc_def = rbm.LocalizationDefinition(list(name_to_valid_compartments[name])) if name in name_to_valid_compartments else None

        mol_defs.append(rbm.MoleculeDefinition(name, mod_defs, assoc_defs, loc_def))

    return mol_defs


def base_reactants(mol_defs: List[rbm.MoleculeDefinition], contingencies: List[con.Contingency],
                   negated_contingencies: List[con.Contingency]):
    return []


def arrow_type_from_reaction(reaction: rxn.Reaction) -> rbm.Arrow:
    pass


def base_rates_from_reaction(reaction: rxn.Reaction) -> List[rbm.Parameter]:
    pass


def derived_rules_from_base_rule_and_contingencies(base_rule: rbm.Rule, contingencies: List[con.Contingency]) -> List[rbm.Rule]:
    pass


def non_overlapping_rules_from_rules(rules: List[rbm.Rule]) -> List[rbm.Rule]:
    pass





