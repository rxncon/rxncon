from typing import List
from collections import defaultdict
import typecheck as tc

import rxncon.core.rxncon_system as rxs
import rxncon.simulation.rule_based.rule_based_model as rbm
import rxncon.core.state as sta


@tc.typecheck
def rule_based_model_from_rxncon(rxnconsys: rxs.RxnConSystem) -> rbm.RuleBasedModel:
    mol_defs = molecule_defs_from_rxncon(rxnconsys)


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
        assoc_defs = name_to_assoc_defs[name] if name in name_to_assoc_defs else None
        mod_defs = name_to_mod_defs if name in name_to_mod_defs else None
        loc_def = rbm.LocalizationDefinition(list(name_to_valid_compartments[name])) if name in name_to_valid_compartments else None

        mol_defs.append(rbm.MoleculeDefinition(name, mod_defs, assoc_defs, loc_def))

    return mol_defs













