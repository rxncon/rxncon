from collections import defaultdict
import typing as tg

import rxncon.core.rxncon_system as rxs
import rxncon.simulation.rule_based.rule_based_model as rbm
import rxncon.core.state as sta
import rxncon.core.reaction as rxn
import rxncon.core.component as com


def molecule_defs_from_rxncon(rxnconsys: rxs.RxnConSystem) -> tg.Dict[str, rbm.MoleculeDefinition]:
    names = set()

    name_to_ass_defs = defaultdict(list)
    name_to_locs = defaultdict(list)
    name_to_mod_defs = defaultdict(list)

    for reaction in rxnconsys.reactions:
        for state in [x for x in [reaction.source, reaction.product] if x]:
            if isinstance(state, sta.CovalentModificationState):
                names.add(state.substrate.name)
                name_to_mod_defs[state.substrate.name].append(_mod_def_from_state_and_reaction(state, reaction))

            elif isinstance(state, sta.InterProteinInteractionState) or isinstance(state, sta.IntraProteinInteractionState):
                names.add(state.first_component.name)
                names.add(state.second_component.name)
                ass_defs = _state_to_ass_defs(state)
                name_to_ass_defs[state.first_component.name].append(ass_defs[0])
                name_to_ass_defs[state.second_component.name].append(ass_defs[1])

            elif isinstance(state, sta.TranslocationState):
                names.add(state.substrate.name)
                name_to_locs[state.substrate.name].append(state.compartment)

            else:
                raise NotImplementedError

    mol_defs = {}

    for name in names:
        mol_defs[name] = rbm.MoleculeDefinition(name,
                                                name_to_mod_defs[name],
                                                name_to_ass_defs[name],
                                                rbm.LocalizationDefinition(name_to_locs[name]))

    return mol_defs


def _mod_def_from_state_and_reaction(state: sta.CovalentModificationState, reaction: rxn.Reaction):
    if not state.substrate.domain and not state.substrate.subdomain and not state.substrate.residue:
        domain_name = 'Mod{}'.format(reaction.subject.name)

    else:
        domain_name = state.substrate.domain

    mod_def = rbm.ModificationDefinition(domain_name, [sta.StateModifier.unmodified.value, state.modifier.value])
    mod_def.matching_state = state

    return mod_def


def _state_to_ass_defs(state: tg.Union[sta.InterProteinInteractionState, sta.IntraProteinInteractionState]):
    first_dom = state.first_component.domain if state.first_component.domain else 'Ass{}'.format(state.second_component.name)
    second_dom = state.second_component.domain if state.second_component.domain else 'Ass{}'.format(state.first_component.name)

    first_def = rbm.AssociationDefinition(first_dom)
    first_def.matching_state = state

    second_def = rbm.AssociationDefinition(second_dom)
    second_def.matching_state = state

    return first_def, second_def




