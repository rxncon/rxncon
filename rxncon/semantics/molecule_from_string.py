from rxncon.semantics.molecule_definition import MoleculeDefinition, \
    AssociationPropertyDefinition, LocalizationPropertyDefinition, ModificationPropertyDefinition, \
    Modifier, Compartment

from rxncon.syntax.rxncon_from_string import component_from_string


def mol_def_from_string(mol_def_string: str) -> MoleculeDefinition:
    name_spec = component_from_string(mol_def_string.split('#')[0])
    def_strings = mol_def_string.split('#')[1].split(',')

    property_defs = [_property_def_from_string(def_string) for def_string in def_strings if def_string]

    ass_defs = {x for x in property_defs if isinstance(x, AssociationPropertyDefinition)}
    mod_defs = {x for x in property_defs if isinstance(x, ModificationPropertyDefinition)}
    loc_defs = {x for x in property_defs if isinstance(x, LocalizationPropertyDefinition)}

    assert len(loc_defs) <= 1
    if len(loc_defs) == 1:
        loc_def = list(loc_defs)[0]
    else:
        loc_def = None

    return MoleculeDefinition(name_spec, mod_defs, ass_defs, loc_def)


def _property_def_from_string(def_string: str):
    identifier = def_string[0:3]
    if identifier == 'ass':
        return _ass_property_def_from_string(def_string)
    elif identifier == 'mod':
        return _mod_property_def_from_string(def_string)
    elif identifier == 'loc':
        return _loc_property_def_from_string(def_string)
    else:
        raise NotImplementedError


def _ass_property_def_from_string(def_string):
    assert def_string[0:4] == 'ass/'
    def_string = def_string[4:]
    ass_domain = component_from_string(def_string.split(':')[0])

    partner_domains = set()

    partner_domain_strings = def_string.split(':')[1].split('~')
    for partner_domain_string in partner_domain_strings:
        partner_domains.add(component_from_string(partner_domain_string))

    return AssociationPropertyDefinition(ass_domain, partner_domains)


def _mod_property_def_from_string(def_string):
    assert def_string[0:4] == 'mod/'
    def_string = def_string[4:]
    mod_domain = component_from_string(def_string.split(':')[0])

    modifiers = set()

    modifier_strings = def_string.split(':')[1].split('~')
    for modifier_string in modifier_strings:
        modifiers.add(Modifier(modifier_string))

    return ModificationPropertyDefinition(mod_domain, modifiers)


def _loc_property_def_from_string(def_string):
    assert def_string[0:4] == 'loc/'
    def_string = def_string[4:]

    compartments = set()

    compartment_strings = def_string.split('~')
    for compartment_string in compartment_strings:
        compartments.add(Compartment(compartment_string))

    return LocalizationPropertyDefinition(compartments)
