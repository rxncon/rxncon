import typing as tg
from collections import defaultdict

from rxncon.core import rxncon_system as rxs, state, state, state, state, specification, state, reaction, state, state, \
    state, state, state, state, state
from rxncon.semantics import molecule_definition as mol
from rxncon.semantics.molecule_instance_from_rxncon import _mod_spec_domain_from_state, _assoc_spec_domain_from_state


class MoleculeDefinitionSupervisor:
    def __init__(self, rxnconsys: rxs.RxnConSystem):
        self.rxnconsys = rxnconsys
        self.molecule_definitions = {}
        self._generate_molecule_definitions()
        self.molecules = self.molecule_definitions.keys()

    def molecule_definition_for_name(self, name: str) -> mol.MoleculeDefinition:
        return self.molecule_definitions[name]

    def _generate_molecule_definitions(self):
        names = set()

        name_to_assoc_defs = defaultdict(set)
        name_to_mod_defs = defaultdict(set)
        name_to_locs = defaultdict(set)

        for reaction in self.rxnconsys.reactions:
            for state in [x for x in [reaction.source, reaction.product] if x]:
                if isinstance(state, sta.CovalentModificationState):
                    names.add(state.substrate.name)
                    names.add(reaction.subject.name)

                    mod_def = _mod_def_from_state_and_reaction(state, reaction)
                    _update_defs(name_to_mod_defs[state.substrate.name], mod_def)

                elif isinstance(state, sta.InterProteinInteractionState) or isinstance(state, sta.IntraProteinInteractionState):
                    names.add(state.first_component.name)
                    names.add(state.second_component.name)

                    assoc_defs = _assoc_defs_from_state(state)
                    _update_defs(name_to_assoc_defs[state.first_component.name], assoc_defs[0])
                    _update_defs(name_to_assoc_defs[state.second_component.name], assoc_defs[1])

                elif isinstance(state, sta.TranslocationState):
                    names.add(state.substrate.name)
                    names.add(reaction.subject.name)
                    name_to_locs[state.substrate.name].add(state.compartment)

                else:
                    raise NotImplementedError

        for name in names:
            universal_specification = spe.Specification(name, None, None, None)
            self.molecule_definitions[name] = mol.MoleculeDefinition(universal_specification,
                                                                     name_to_mod_defs[name],
                                                                     name_to_assoc_defs[name],
                                                                     mol.LocalizationDefinition(name_to_locs[name]))


def _mod_def_from_state_and_reaction(state: sta.CovalentModificationState, reaction: rxn.Reaction):

    state = _mod_spec_domain_from_state(state, reaction)
    mod_def = mol.ModificationDefinition(state.substrate,
                                         {mol.Modifier.unmodified, _molecule_modifier_from_state_modifier(state.modifier)})

    return mod_def


def _assoc_defs_from_state(state: tg.Union[sta.InterProteinInteractionState, sta.IntraProteinInteractionState]):
    first_spec, second_spec = _assoc_spec_domain_from_state(state)

    first_def = mol.AssociationDefinition(first_spec, {second_spec})
    second_def = mol.AssociationDefinition(second_spec, {first_spec})

    return first_def, second_def


def _molecule_modifier_from_state_modifier(state_mod: sta.StateModifier) -> mol.Modifier:
    if state_mod == sta.StateModifier.phosphor:
        return mol.Modifier.phosphorylated

    elif state_mod == sta.StateModifier.ubiquitin:
        return mol.Modifier.ubiquitinated

    elif state_mod == sta.StateModifier.truncated:
        return mol.Modifier.truncated

    elif state_mod == sta.StateModifier.unmodified:
        raise NotImplementedError

    else:
        raise NotImplementedError


def _update_defs(defs: tg.Set[mol.Definition], new_def: mol.Definition):
    if isinstance(new_def, mol.AssociationDefinition):
        matches = lambda x, y: x == y
        update = lambda old_def, new_def: old_def.valid_partners.update(new_def.valid_partners)
    else:
        matches = lambda x, y: x.is_equivalent_to(y)
        update = lambda old_def, new_def: old_def.valid_modifiers.update(new_def.valid_modifiers)

    # if update is used the function will return a None which is stored in found_updateable_def resulting in an non-empty list
    found_updatable_def = [update(the_def, new_def) for the_def in defs if matches(the_def.spec, new_def.spec)]

    if not found_updatable_def:
        defs.add(new_def)