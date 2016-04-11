from rxncon.semantics.molecule_definition import MoleculeDefinition, \
    AssociationPropertyDefinition, LocalizationPropertyDefinition, ModificationPropertyDefinition, \
    Modifier, Compartment, OccupationStatus
from rxncon.semantics.molecule_instance import MoleculeInstance, \
    AssociationPropertyInstance, LocalizationPropertyInstance, ModificationPropertyInstance, Binding

from rxncon.syntax.rxncon_from_string import specification_from_string
from rxncon.simulation.rule_based.rule_based_model import Binding, Rule, MoleculeReactant, ComplexReactant, Arrow, \
    Parameter


def mol_def_from_string(mol_def_string: str) -> MoleculeDefinition:
    name_spec = specification_from_string(mol_def_string.split('#')[0])
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


def mol_instance_from_string(mol_def, mol_ins_string: str) -> MoleculeInstance:
    if isinstance(mol_def, str):
        mol_def = mol_def_from_string(mol_def)

    assert isinstance(mol_def, MoleculeDefinition)

    property_instance_strings = mol_ins_string.split('#')[1].split(',')
    property_instances = [_property_ins_from_string(mol_def, x) for x in property_instance_strings if x]

    mod_props = {x for x in property_instances if isinstance(x, ModificationPropertyInstance)}
    ass_props = {x for x in property_instances if isinstance(x, AssociationPropertyInstance)}
    loc_props = {x for x in property_instances if isinstance(x, LocalizationPropertyInstance)}

    assert len(loc_props) <= 1
    if loc_props:
        loc_prop = list(loc_props)[0]
    else:
        loc_prop = None

    return MoleculeInstance(mol_def, mod_props, ass_props, loc_prop)


def mol_instances_and_bindings_from_string(mol_defs, mol_instances_string: str):
    mol_instances = []
    bindings = []

    mol_def_dict = {}

    for mol_def in mol_defs:
        if isinstance(mol_def, str):
            mol_def = mol_def_from_string(mol_def)
        mol_def_dict[mol_def.spec] = mol_def

    mol_ins_strings = mol_instances_string.split('.')

    for mol_ins_string in mol_ins_strings:
        mol_def = mol_def_dict[specification_from_string(mol_ins_string.split('#')[0])]
        mol_ins = mol_instance_from_string(mol_def, mol_ins_string)
        mol_instances.append(mol_ins)

    finished_bindings_nums = []

    for i, left_mol_ins in enumerate(mol_instances):
        for left_ass_property in left_mol_ins.association_properties:
            if hasattr(left_ass_property, 'binding_num') and left_ass_property.binding_num not in finished_bindings_nums:
                left_partner = (i, left_ass_property)

                for j, right_mol_ins in enumerate(mol_instances):
                    for right_ass_property in right_mol_ins.association_properties:
                        if hasattr(right_ass_property, 'binding_num') and left_ass_property.binding_num == right_ass_property.binding_num:
                            right_partner = (j, right_ass_property)

                assert right_partner
                bindings.append(Binding(left_partner, right_partner))

                finished_bindings_nums.append(left_ass_property.binding_num)

    return mol_instances, bindings


def rule_from_string(mol_defs, rule_string):
    mol_def_dict = {}

    for mol_def in mol_defs:
        if isinstance(mol_def, str):
            mol_def = mol_def_from_string(mol_def)
        mol_def_dict[mol_def.spec] = mol_def

    if '<->' in rule_string:
        split = rule_string.split('<->')
        arrow = Arrow('<->')
        parameters = [Parameter('k_f', None), Parameter('k_r', None)]
    elif '->' in rule_string:
        split = rule_string.split('->')
        arrow = Arrow('->')
        parameters = [Parameter('k', None)]
    else:
        raise AssertionError

    lhs_string, rhs_string = split[0].strip(), split[1].strip()

    lhs_terms = [x.strip() for x in lhs_string.split('+')]
    rhs_terms = [x.strip() for x in rhs_string.split('+')]

    lhs_reactants = []
    for term in lhs_terms:
        if '.' in term:
            instances, bindings = mol_instances_and_bindings_from_string(mol_defs, term)
            lhs_reactants.append(ComplexReactant(instances, bindings))
        else:
            instance = mol_instance_from_string(mol_def_dict[specification_from_string(term.split('#')[0])], term)
            lhs_reactants.append(MoleculeReactant(instance))

    rhs_reactants = []
    for term in rhs_terms:
        if '.' in term:
            instances, bindings = mol_instances_and_bindings_from_string(mol_defs, term)
            rhs_reactants.append(ComplexReactant(instances, bindings))
        else:
            instance = mol_instance_from_string(mol_def_dict[specification_from_string(term.split('#')[0])], term)
            rhs_reactants.append(MoleculeReactant(instance))

    return Rule(lhs_reactants, rhs_reactants, arrow, parameters)


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
    ass_domain = specification_from_string(def_string.split(':')[0])

    partner_domains = set()

    partner_domain_strings = def_string.split(':')[1].split('~')
    for partner_domain_string in partner_domain_strings:
        partner_domains.add(specification_from_string(partner_domain_string))

    return AssociationPropertyDefinition(ass_domain, partner_domains)


def _mod_property_def_from_string(def_string):
    assert def_string[0:4] == 'mod/'
    def_string = def_string[4:]
    mod_domain = specification_from_string(def_string.split(':')[0])

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


def _property_ins_from_string(mol_def, prop_string):
    identifier = prop_string[0:3]
    if identifier == 'ass':
        if len(prop_string.split('~')) == 2:
            binding_num = prop_string.split('~')[1]
            prop_string = prop_string.split('~')[0]
            ass_prop_instance = _ass_property_ins_from_string(mol_def, prop_string)
            ass_prop_instance.binding_num = binding_num
            return ass_prop_instance
        else:
            return _ass_property_ins_from_string(mol_def, prop_string)

    elif identifier == 'mod':
        return _mod_property_ins_from_string(mol_def, prop_string)
    elif identifier == 'loc':
        return _loc_property_ins_from_string(mol_def, prop_string)
    else:
        raise NotImplementedError


def _ass_property_ins_from_string(mol_def, prop_string):
    assert prop_string[0:4] == 'ass/'
    prop_string = prop_string[4:]

    assert len(prop_string.split(':')) == 2

    ass_domain = specification_from_string(prop_string.split(':')[0])
    ass_defs = [x for x in mol_def.association_defs if x.spec == ass_domain]
    assert len(ass_defs) == 1
    ass_def = ass_defs[0]

    partner_domain_string = prop_string.split(':')[1]
    if partner_domain_string:
        return AssociationPropertyInstance(ass_def, OccupationStatus.occupied_known_partner, specification_from_string(partner_domain_string))
    else:
        return AssociationPropertyInstance(ass_def, OccupationStatus.not_occupied, None)


def _mod_property_ins_from_string(mol_def, prop_string):
    assert prop_string[0:4] == 'mod/'
    prop_string = prop_string[4:]

    mod_domain = specification_from_string(prop_string.split(':')[0])
    mod_defs = [x for x in mol_def.modification_defs if x.spec == mod_domain]
    assert len(mod_defs) == 1
    mod_def = mod_defs[0]

    mod_string = prop_string.split(':')[1]
    return ModificationPropertyInstance(mod_def, Modifier(mod_string))


def _loc_property_ins_from_string(mol_def, prop_string):
    assert prop_string[0:4] == 'mod/'
    prop_string = prop_string[4:]

    return ModificationPropertyInstance(mol_def.localization_def, Modifier(prop_string))







