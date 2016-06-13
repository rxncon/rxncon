import typing as tg
from collections import defaultdict
import copy

import rxncon.core.rxncon_system as rxs
import rxncon.semantics.molecule
import rxncon.core.state as sta
import rxncon.core.reaction as rxn
import rxncon.core.specification as spe


def mol_defs_from_rxncon_sys(rxnconsys: rxs.RxnConSystem) -> tg.Dict[spe.Specification, rxncon.semantics.molecule.MoleculeDefinition]:
    return _MoleculeDefinitionSupervisor(rxnconsys).molecule_definitions


def mol_modifier_from_state_modifier(state_mod: sta.StateModifier) -> rxncon.semantics.molecule.Modifier:
    if state_mod == sta.StateModifier.phosphor:
        return rxncon.semantics.molecule.Modifier.phosphorylated
    elif state_mod == sta.StateModifier.ubiquitin:
        return rxncon.semantics.molecule.Modifier.ubiquitinated
    elif state_mod == sta.StateModifier.truncated:
        return rxncon.semantics.molecule.Modifier.truncated
    elif state_mod == sta.StateModifier.guanosintriphosphat:
        return rxncon.semantics.molecule.Modifier.guanosintriphosphat
    else:
        raise NotImplementedError


def mod_domain_spec_from_state_and_reaction(state: sta.CovalentModificationState, reaction: rxn.Reaction):
    # Need to copy, since otherwise we will mutate the specs appearing in the original state/reaction.
    spec = copy.copy(state.substrate)

    if not spec.residue and reaction.influence == rxn.Influence.transfer:
        if spec == reaction.subject:
            spec.residue = _kinase_residue_name(reaction.object)
        elif spec == reaction.object:
            spec.residue = _kinase_residue_name(reaction.subject)
        else:
            raise NotImplementedError
    elif not spec.residue and reaction.influence in [rxn.Influence.positive, rxn.Influence.negative]:
        spec.residue = _kinase_residue_name(reaction.subject)

    return spec


def ass_domain_specs_from_state(state: tg.Union[sta.InteractionState, sta.SelfInteractionState]):
    # Need to copy, since otherwise we will mutate the specs appearing in the original state/reaction.
    first_spec = copy.copy(state.first_component)
    second_spec = copy.copy(state.second_component)

    if not first_spec.domain:
        first_spec.domain = _assoc_domain_from_partner_spec(state.second_component)

    if not second_spec.domain:
        second_spec.domain = _assoc_domain_from_partner_spec(state.first_component)

    return first_spec, second_spec


class _MoleculeDefinitionSupervisor:
    def __init__(self, rxnconsys: rxs.RxnConSystem):
        self.rxnconsys = rxnconsys
        self.molecule_definitions = {}
        self._generate_molecule_definitions()

    def mol_def_for_name(self, name: str) -> rxncon.semantics.molecule.MoleculeDefinition:
        return self.molecule_definitions[name]

    def _generate_molecule_definitions(self):
        specs = set()

        spec_to_assoc_defs = defaultdict(set)
        spec_to_mod_defs = defaultdict(set)
        spec_to_locs = defaultdict(set)

        for reaction in self.rxnconsys.reactions:
            for state in [x for x in reaction.sources + reaction.products if x]:
                if isinstance(state, sta.CovalentModificationState):
                    specs.add(state.substrate.to_component_specification())
                    specs.add(reaction.subject.to_component_specification())

                    mod_def = _mod_def_from_state_and_reaction(state, reaction)
                    _update_defs(spec_to_mod_defs[state.substrate.to_component_specification()], mod_def)
                elif isinstance(state, sta.InteractionState) or isinstance(state, sta.SelfInteractionState):
                    specs.add(state.first_component.to_component_specification())
                    specs.add(state.second_component.to_component_specification())

                    assoc_defs = _assoc_defs_from_state(state)
                    _update_defs(spec_to_assoc_defs[state.first_component.to_component_specification()], assoc_defs[0])
                    _update_defs(spec_to_assoc_defs[state.second_component.to_component_specification()], assoc_defs[1])
                elif isinstance(state, sta.TranslocationState):
                    specs.add(state.substrate.to_component_specification())
                    specs.add(reaction.subject.to_component_specification())
                    spec_to_locs[state.substrate.to_component_specification()].add(state.compartment)
                else:
                    raise NotImplementedError

        for spec in specs:
            self.molecule_definitions[spec] = rxncon.semantics.molecule.MoleculeDefinition(spec,
                                                                                           spec_to_mod_defs[spec],
                                                                                           spec_to_assoc_defs[spec],
                                                                                           rxncon.semantics.molecule.LocalizationPropertyDefinition(spec_to_locs[spec]))


def _mod_def_from_state_and_reaction(state: sta.CovalentModificationState, reaction: rxn.Reaction):
    spec = mod_domain_spec_from_state_and_reaction(state, reaction)
    mod_def = rxncon.semantics.molecule.ModificationPropertyDefinition(spec,
                                                                       {rxncon.semantics.molecule.Modifier.unmodified, mol_modifier_from_state_modifier(state.modifier)})

    return mod_def


def _assoc_defs_from_state(state: tg.Union[sta.InteractionState, sta.SelfInteractionState]):
    first_spec, second_spec = ass_domain_specs_from_state(state)

    first_def = rxncon.semantics.molecule.AssociationPropertyDefinition(first_spec, {second_spec})
    second_def = rxncon.semantics.molecule.AssociationPropertyDefinition(second_spec, {first_spec})

    return first_def, second_def


def _update_defs(defs: tg.Set[rxncon.semantics.molecule.PropertyDefinition], new_def: rxncon.semantics.molecule.PropertyDefinition):
    if isinstance(new_def, rxncon.semantics.molecule.AssociationPropertyDefinition):
        matches = lambda x, y: x == y
        update = lambda old_def, new_def: old_def.valid_partners.update(new_def.valid_partners)
    else:
        matches = lambda x, y: x.is_equivalent_to(y)
        update = lambda old_def, new_def: old_def.valid_modifiers.update(new_def.valid_modifiers)

    # if update is used the function will return a None which is stored in found_updateable_def resulting in an non-empty list
    found_updatable_def = [update(the_def, new_def) for the_def in defs if matches(the_def.spec, new_def.spec)]

    if not found_updatable_def:
        defs.add(new_def)


def _kinase_residue_name(spec: spe.Specification) -> str:
    return '{0}site'.format(spec.name)


def _assoc_domain_from_partner_spec(spec: spe.Specification) -> str:
    return '{0}assoc'.format(spec.name)
