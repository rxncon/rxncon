from rxncon.semantics.molecule_definition import MoleculeDefinition, AssociationPropertyDefinition, LocalizationPropertyDefinition, ModificationPropertyDefinition, Modifier, Compartment, OccupationStatus
from rxncon.semantics.molecule_instance import MoleculeInstance, AssociationPropertyInstance, LocalizationPropertyInstance, ModificationPropertyInstance, Binding

from rxncon.syntax.rxncon_from_string import specification_from_string
from rxncon.simulation.rule_based.rule_based_model import Rule, MoleculeReactant, ComplexReactant, Arrow, Parameter

from typecheck import typecheck
import typing as tg

ID_ASS = 'ass/'
ID_MOD = 'mod/'
ID_LOC = 'loc/'

TOK_MOL_SEP  = '#'
TOK_PROP_SEP = ','
TOK_PROP_VAL = ':'
TOK_RATE_SEP = '@'
TOK_IRREV    = '->'
TOK_REV      = '<->'
TOK_BIND     = '.'


@typecheck
def mol_def_from_string(mol_def_string: str) -> MoleculeDefinition:
    def _property_def_from_string(def_string: str):
        def _ass_property_def_from_string(def_string):
            assert def_string[0:4] == ID_ASS
            def_string = def_string[4:]
            ass_domain = specification_from_string(def_string.split(TOK_PROP_VAL)[0])

            partner_domains = set()

            partner_domain_strings = def_string.split(TOK_PROP_VAL)[1].split('~')
            for partner_domain_string in partner_domain_strings:
                partner_domains.add(specification_from_string(partner_domain_string))

            return AssociationPropertyDefinition(ass_domain, partner_domains)

        def _mod_property_def_from_string(def_string):
            assert def_string[0:4] == ID_MOD
            def_string = def_string[4:]
            mod_domain = specification_from_string(def_string.split(TOK_PROP_VAL)[0])

            modifiers = set()

            modifier_strings = def_string.split(TOK_PROP_VAL)[1].split('~')
            for modifier_string in modifier_strings:
                modifiers.add(Modifier(modifier_string))

            return ModificationPropertyDefinition(mod_domain, modifiers)

        def _loc_property_def_from_string(def_string):
            assert def_string[0:4] == ID_LOC
            def_string = def_string[4:]

            compartments = set()

            compartment_strings = def_string.split('~')
            for compartment_string in compartment_strings:
                compartments.add(Compartment(compartment_string))

            return LocalizationPropertyDefinition(compartments)

        identifier = def_string[0:4]
        if identifier == ID_ASS:
            return _ass_property_def_from_string(def_string)
        elif identifier == ID_MOD:
            return _mod_property_def_from_string(def_string)
        elif identifier == ID_LOC:
            return _loc_property_def_from_string(def_string)
        else:
            raise NotImplementedError

    name_spec = specification_from_string(mol_def_string.split(TOK_MOL_SEP)[0])

    property_defs = [_property_def_from_string(def_string) for def_string in
                     mol_def_string.split(TOK_MOL_SEP)[1].split(TOK_PROP_SEP) if def_string]

    ass_defs = {x for x in property_defs if isinstance(x, AssociationPropertyDefinition)}
    mod_defs = {x for x in property_defs if isinstance(x, ModificationPropertyDefinition)}
    loc_defs = {x for x in property_defs if isinstance(x, LocalizationPropertyDefinition)}

    assert len(loc_defs) <= 1, 'Number of LocalizationPropertyDefinition in {0} exceeds one.'.format(mol_def_string)
    loc_def = list(loc_defs)[0] if len(loc_defs) == 1 else None

    return MoleculeDefinition(name_spec, mod_defs, ass_defs, loc_def)


@typecheck
def mol_instance_from_string(mol_def: tg.Union[str, MoleculeDefinition], mol_instance_string: str) -> MoleculeInstance:
    def _property_ins_from_string(mol_def, prop_string):
        def _ass_property_ins_from_string(mol_def, prop_string):
            assert prop_string[0:4] == ID_ASS
            prop_string = prop_string[4:]

            assert len(prop_string.split(TOK_PROP_VAL)) == 2

            ass_domain = specification_from_string(prop_string.split(TOK_PROP_VAL)[0])
            ass_defs = [x for x in mol_def.association_defs if x.spec == ass_domain]
            assert len(ass_defs) == 1
            ass_def = ass_defs[0]

            partner_domain_string = prop_string.split(TOK_PROP_VAL)[1]
            if partner_domain_string:
                return AssociationPropertyInstance(ass_def, OccupationStatus.occupied_known_partner,
                                                   specification_from_string(partner_domain_string))
            else:
                return AssociationPropertyInstance(ass_def, OccupationStatus.not_occupied, None)

        def _mod_property_ins_from_string(mol_def, prop_string):
            assert prop_string[0:4] == ID_MOD
            prop_string = prop_string[4:]

            mod_domain = specification_from_string(prop_string.split(TOK_PROP_VAL)[0])
            mod_defs = [x for x in mol_def.modification_defs if x.spec == mod_domain]
            assert len(mod_defs) == 1
            mod_def = mod_defs[0]

            mod_string = prop_string.split(TOK_PROP_VAL)[1]
            return ModificationPropertyInstance(mod_def, Modifier(mod_string))

        def _loc_property_ins_from_string(mol_def, prop_string):
            assert prop_string[0:4] == ID_LOC
            prop_string = prop_string[4:]

            return ModificationPropertyInstance(mol_def.localization_def, Modifier(prop_string))

        identifier = prop_string[0:4]
        if identifier == ID_ASS:
            if len(prop_string.split('~')) == 2:
                binding_num = prop_string.split('~')[1]
                prop_string = prop_string.split('~')[0]
                ass_prop_instance = _ass_property_ins_from_string(mol_def, prop_string)
                ass_prop_instance.binding_num = binding_num
                return ass_prop_instance
            else:
                return _ass_property_ins_from_string(mol_def, prop_string)
        elif identifier == ID_MOD:
            return _mod_property_ins_from_string(mol_def, prop_string)
        elif identifier == ID_LOC:
            return _loc_property_ins_from_string(mol_def, prop_string)
        else:
            raise NotImplementedError('Unknown property type identifier {0}'.format(identifier))

    if isinstance(mol_def, str):
        mol_def = mol_def_from_string(mol_def)

    assert isinstance(mol_def, MoleculeDefinition)

    property_instances = [_property_ins_from_string(mol_def, prop_string) for prop_string in
                          mol_instance_string.split(TOK_MOL_SEP)[1].split(TOK_PROP_SEP) if prop_string]

    mod_props = {x for x in property_instances if isinstance(x, ModificationPropertyInstance)}
    ass_props = {x for x in property_instances if isinstance(x, AssociationPropertyInstance)}
    loc_props = {x for x in property_instances if isinstance(x, LocalizationPropertyInstance)}

    assert len(loc_props) <= 1, 'Number of LocalizationPropertyInstance in {0} exceeds one.'.format(mol_instance_string)
    loc_prop = list(loc_props)[0] if len(loc_props) == 1 else None

    return MoleculeInstance(mol_def, mod_props, ass_props, loc_prop)


@typecheck
def mol_instances_and_bindings_from_string(mol_defs: tg.Union[tg.Iterable[str], tg.Iterable[MoleculeDefinition]], mol_instances_string: str) \
        -> tg.Tuple[tg.Set[MoleculeInstance], tg.Set[Binding]]:
    mol_instances = set()
    bindings = set()

    mol_def_dict = {}

    for mol_def in mol_defs:
        if isinstance(mol_def, str):
            mol_def = mol_def_from_string(mol_def)
        assert isinstance(mol_def, MoleculeDefinition)
        mol_def_dict[mol_def.spec] = mol_def

    mol_ins_strings = mol_instances_string.split(TOK_BIND)

    for mol_ins_string in mol_ins_strings:
        mol_def = mol_def_dict[specification_from_string(mol_ins_string.split(TOK_MOL_SEP)[0])]
        mol_ins = mol_instance_from_string(mol_def, mol_ins_string)
        mol_instances.add(mol_ins)

    for mol_ins in mol_instances:
        for ass_prop in mol_ins.association_properties:
            if ass_prop.partner:
                bindings.add(Binding(ass_prop.property_def.spec, ass_prop.partner))

    return mol_instances, bindings


@typecheck
def rule_from_string(mol_defs, rule_string: str) -> Rule:
    def params_from_string(param_string):
        return {Parameter(name.strip(), None) for name in param_string.split(TOK_PROP_SEP)}

    mol_def_dict = {}

    for mol_def in mol_defs:
        if isinstance(mol_def, str):
            mol_def = mol_def_from_string(mol_def)
        mol_def_dict[mol_def.spec] = mol_def

    assert TOK_RATE_SEP in rule_string, 'Rate constants missing from rule string.'

    if TOK_REV in rule_string:
        split, param_string = rule_string.split(TOK_RATE_SEP)
        split = split.split(TOK_REV)
        arrow = Arrow(TOK_REV)
        parameters = params_from_string(param_string)
    elif TOK_IRREV in rule_string:
        split, param_string = rule_string.split(TOK_RATE_SEP)
        split = split.split(TOK_IRREV)
        arrow = Arrow(TOK_IRREV)
        parameters = params_from_string(param_string)
    else:
        raise AssertionError

    lhs_string, rhs_string = split[0].strip(), split[1].strip()

    lhs_terms = [x.strip() for x in lhs_string.split('+')]
    rhs_terms = [x.strip() for x in rhs_string.split('+')]

    lhs_reactants = set()
    for term in lhs_terms:
        if TOK_BIND in term:
            instances, bindings = mol_instances_and_bindings_from_string(mol_defs, term)
            lhs_reactants.add(ComplexReactant(instances, bindings))
        else:
            instance = mol_instance_from_string(mol_def_dict[specification_from_string(term.split(TOK_MOL_SEP)[0])], term)
            lhs_reactants.add(MoleculeReactant(instance))

    rhs_reactants = set()
    for term in rhs_terms:
        if TOK_BIND in term:
            instances, bindings = mol_instances_and_bindings_from_string(mol_defs, term)
            rhs_reactants.add(ComplexReactant(instances, bindings))
        else:
            instance = mol_instance_from_string(mol_def_dict[specification_from_string(term.split(TOK_MOL_SEP)[0])], term)
            rhs_reactants.add(MoleculeReactant(instance))

    return Rule(lhs_reactants, rhs_reactants, arrow, parameters)

