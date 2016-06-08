import typing as tp
import typing as tg
from enum import unique

import typecheck as tc

from rxncon.core import specification as spe
from rxncon.core.specification import Specification
from rxncon.util.utils import OrderedEnum


@unique
class Modifier(OrderedEnum):
    unmodified          = 'u'
    phosphorylated      = 'p'
    ubiquitinated       = 'ub'
    truncated           = 'truncated'
    guanosintriphosphat = 'gtp'


@unique
class OccupationStatus(OrderedEnum):
    not_specified = 0
    not_occupied = 1
    occupied_known_partner = 2
    occupied_unknown_partner = 3


@unique
class Compartment(OrderedEnum):
    cell = 'cell'
    cytosole = 'cytosole'
    nucleus = 'nucleus'


class MutualExclusivityError(Exception):
    pass


class MoleculeDefinition:
    @tc.typecheck
    def __init__(self, spec: spe.Specification,
                 modification_defs: tp.Set['ModificationPropertyDefinition'],
                 association_defs: tp.Set['AssociationPropertyDefinition'],
                 localization_def: tp.Optional['LocalizationPropertyDefinition']):
        self.spec = spec
        self.modification_defs = modification_defs
        self.association_defs = association_defs

        if localization_def is None:
            localization_def = LocalizationPropertyDefinition(set())

        self.localization_def = localization_def

        self._validate()

    def _validate(self):
        def definitions_validation(definitions: tp.Union[tp.Set[ModificationPropertyDefinition],
                                                       tp.Set[AssociationPropertyDefinition]]):
            definitions = list(definitions)
            i = 0
            while i < len(definitions):
                for definition in definitions[i+1:]:
                    if definitions[i].spec == definition.spec:
                        raise AssertionError('Same specification in different context: {0} and {1}'.format(definitions[i], definition))
                i += 1

        definitions_validation(self.association_defs)
        definitions_validation(self.modification_defs)


    def __hash__(self) -> int:
        return hash(str(self.spec))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'MoleculeDefinition: {0} [Mod: {1}, Ass: {2}, Loc: {3}]'\
            .format(self.spec,
                    '/'.join(str(x) for x in sorted(self.modification_defs)),
                    '/'.join(str(x) for x in sorted(self.association_defs)),
                    str(self.localization_def))

    @tc.typecheck
    def __eq__(self, other: 'MoleculeDefinition') -> bool:
        return isinstance(other, MoleculeDefinition) and self.spec == other.spec and self.localization_def == other.localization_def and \
            other.modification_defs == self.modification_defs and other.association_defs == self.association_defs

    def __lt__(self, other: 'MoleculeDefinition'):
        return self.spec < other.spec


class PropertyDefinition:
    pass


class ModificationPropertyDefinition(PropertyDefinition):
    @tc.typecheck
    def __init__(self, spec: spe.Specification, valid_modifiers: tp.Set['Modifier']):
        self.spec = spec
        self.valid_modifiers = valid_modifiers

    def __hash__(self) -> int:
        return hash('*mod-def* {0}'.format(self.spec.name))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'ModificationPropertyDefinition: Domain = {0}, Modifiers = {1}'\
            .format(self.spec, ', '.join(mod.value for mod in sorted(self.valid_modifiers, key=str)))

    @tc.typecheck
    def __eq__(self, other: PropertyDefinition):
        return isinstance(other, ModificationPropertyDefinition) and self.spec == other.spec and \
            self.valid_modifiers == other.valid_modifiers

    def __lt__(self, other: 'ModificationPropertyDefinition'):
        return self.spec < other.spec


class AssociationPropertyDefinition(PropertyDefinition):
    @tc.typecheck
    def __init__(self, spec: spe.Specification, valid_partners: tp.Set[spe.Specification]):
        self.spec = spec
        self.valid_partners = valid_partners

    def __hash__(self) -> int:
        return hash('*ass-def* {}'.format(self.spec.name))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'AssociationDefinition: Domain = {0}, valid_partners = {1}'\
            .format(self.spec, ', '.join(str(x) for x in sorted(self.valid_partners)))

    @tc.typecheck
    def __eq__(self, other: PropertyDefinition) -> bool:
        return isinstance(other, AssociationPropertyDefinition) and self.spec == other.spec and \
            self.valid_partners == other.valid_partners

    def __lt__(self, other: 'AssociationPropertyDefinition'):
        return  self.spec < other.spec


class LocalizationPropertyDefinition(PropertyDefinition):
    @tc.typecheck
    def __init__(self, valid_compartments: tp.Set[Compartment]):
        self.valid_compartments = valid_compartments

    def __hash__(self) -> int:
        return hash('*loc-def* with num of compartments {}'.format(len(self.valid_compartments)))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'LocalizationDefinition: {0}'.format(', '.join(str(x) for x in sorted(self.valid_compartments)))


    @tc.typecheck
    def __eq__(self, other: PropertyDefinition):
        return isinstance(other, LocalizationPropertyDefinition) and self.valid_compartments == other.valid_compartments

    def __lt__(self, other: 'LocalizationPropertyDefinition') -> bool:
        return sorted(list(self.valid_compartments)) < sorted(list(other.valid_compartments))


class MoleculeInstance:
    @tc.typecheck
    def __init__(self,
                 mol_def: MoleculeDefinition,
                 modification_properties: tg.Set['ModificationPropertyInstance'],
                 association_properties: tg.Set['AssociationPropertyInstance'],
                 localization_property: tg.Optional['LocalizationPropertyInstance']):
        self.mol_def = mol_def
        self.modification_properties = modification_properties
        self.association_properties = association_properties
        self.localization_property = localization_property

        self._validate()

    @tc.typecheck
    def __eq__(self, other: 'MoleculeInstance'):
        return isinstance(other, MoleculeInstance) and \
            self.mol_def == other.mol_def and \
            self.localization_property == other.localization_property and \
            self.modification_properties == other.modification_properties and \
            self.association_properties == other.association_properties

    @tc.typecheck
    def __lt__(self, other: 'MoleculeInstance') -> bool:
        if self.mol_def.spec < other.mol_def.spec:
            return True
        elif self.mol_def.spec == other.mol_def.spec:
            if sorted(self.modification_properties) < sorted(other.modification_properties):
                return True
            if sorted(self.association_properties) < sorted(other.association_properties):
                return True
            if self.localization_property is not None and other.localization_property is not None \
                    and self.localization_property < other.localization_property:
                return True

        return False

    def __hash__(self) -> int:
        return hash(str(self.mol_def.spec))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        props = sorted(self.modification_properties) + sorted(self.association_properties)
        if self.localization_property:
            props += [self.localization_property]

        return '{0}#{1}'\
            .format(self.mol_def.spec, ','.join(str(x) for x in props))

    def add_property(self, property_instance: PropertyInstance):
        if isinstance(property_instance, ModificationPropertyInstance):
            occupied_residues = [x.property_def.spec for x in self.modification_properties]
            if property_instance.property_def.spec in occupied_residues:
                raise MutualExclusivityError
            self.modification_properties.add(property_instance)
        elif isinstance(property_instance, AssociationPropertyInstance):
            occupied_domains = [x.property_def.spec for x in self.association_properties]
            if property_instance.property_def.spec in occupied_domains:
                raise MutualExclusivityError
            self.association_properties.add(property_instance)
        elif isinstance(property_instance, LocalizationPropertyInstance):
            if self.localization_property:
                raise MutualExclusivityError
            self.localization_property = property_instance
        else:
            raise NotImplementedError

        self._validate()

    @property
    def bindings(self) -> tg.Set['Binding']:
        return {Binding(ass_prop.property_def.spec, ass_prop.partner) for ass_prop in self.association_properties
                if ass_prop.occupation_status == OccupationStatus.occupied_known_partner}

    @property
    def properties(self):
        if self.localization_property:
            return self.association_properties.union(self.modification_properties).union({self.localization_property})
        else:
            return self.association_properties.union(self.modification_properties)

    def _validate(self):
        # Assert each modification domain and each association domain is only present once.
        assert len([mod_prop.property_def for mod_prop in self.modification_properties]) == \
            len(set([mod_prop.property_def for mod_prop in self.modification_properties]))

        assert len([ass_prop.property_def for ass_prop in self.association_properties]) == \
            len(set([ass_prop.property_def for ass_prop in self.association_properties]))


class PropertyInstance:
    pass


class ModificationPropertyInstance(PropertyInstance):
    @tc.typecheck
    def __init__(self, property_def: ModificationPropertyDefinition, modifier: 'Modifier'):
        self.property_def = property_def
        self.modifier = modifier
        self._validate()

    @tc.typecheck
    def __eq__(self, other: PropertyInstance) -> bool:
        return isinstance(other, ModificationPropertyInstance) and self.property_def == other.property_def and \
            self.modifier == other.modifier

    def __lt__(self, other: 'ModificationPropertyInstance'):
        if self.property_def.spec < other.property_def.spec:
            return True

        return self.modifier.value < other.modifier.value

    def __hash__(self) -> bool:
        return hash(str(self.property_def.spec))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'mod/{0}:{1}'.format(self.property_def.spec, self.modifier)

    def _validate(self):
        if self.modifier not in self.property_def.valid_modifiers:
            raise ValueError('Modifier {0} does not appear in Set of valid modifiers for domain {1}.'
                             .format(self.modifier, self.property_def.spec))

    @property
    def complementary_instances(self):
        return [ModificationPropertyInstance(self.property_def, modifier) for modifier
                in self.property_def.valid_modifiers if modifier != self.modifier]


class AssociationPropertyInstance(PropertyInstance):
    @tc.typecheck
    def __init__(self, property_def: AssociationPropertyDefinition, occupation_status: 'OccupationStatus',
                 partner: tg.Optional[spe.Specification]):
        self.property_def = property_def
        self.occupation_status = occupation_status
        self.partner = partner

        self._validate()

    @tc.typecheck
    def __eq__(self, other: PropertyInstance) -> bool:
        if not isinstance(other, AssociationPropertyInstance):
            return False

        if self.partner and other.partner:
            return self.property_def == other.property_def and self.occupation_status == other.occupation_status and self.partner == other.partner
        elif not self.partner and not other.partner:
            return self.property_def == other.property_def and self.occupation_status == other.occupation_status
        else:
            return False

    def __lt__(self, other: 'AssociationPropertyInstance'):
        if self.property_def.spec < other.property_def.spec:
            return True
        return False

    def __hash__(self) -> int:
        return hash(str(self.property_def.spec))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'ass/{0}:{1}'.format(self.property_def.spec, self.partner)

    @property
    def complementary_instances(self):
        if self.occupation_status == OccupationStatus.occupied_known_partner:
            unoccupied = AssociationPropertyInstance(self.property_def, OccupationStatus.not_occupied, None)
            other_partners = [AssociationPropertyInstance(self.property_def, OccupationStatus.occupied_known_partner, x)
                              for x in self.property_def.valid_partners if x != self.partner]

            return other_partners + [unoccupied]
        elif self.occupation_status == OccupationStatus.not_occupied:
            return [AssociationPropertyInstance(self.property_def, OccupationStatus.occupied_known_partner, x)
                    for x in self.property_def.valid_partners]
        else:
            raise NotImplementedError

    def _validate(self):
        if self.partner:
            assert self.partner in self.property_def.valid_partners, \
                'For association domain {0}, binding partner {1} does not appear in list of valid partners {2}.' \
                .format(str(self.property_def.spec), str(self.partner), ','.join(str(x) for x in sorted(self.property_def.valid_partners)))

        if self.occupation_status == OccupationStatus.not_occupied:
            assert not self.partner


class LocalizationPropertyInstance(PropertyInstance):
    @tc.typecheck
    def __init__(self, property_def: LocalizationPropertyDefinition, compartment: Compartment):
        self.property_def = property_def
        self.compartment = compartment
        self._validate()

    @tc.typecheck
    def __eq__(self, other: PropertyInstance) -> bool:
        return isinstance(other, LocalizationPropertyInstance) and self.property_def == other.localization_def and \
            self.compartment == other.compartment

    def __lt__(self, other: 'LocalizationPropertyInstance') -> bool:
        return self.compartment < other.compartment

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'loc/{0}'.format(self.compartment)

    def _validate(self):
        if self.compartment not in self.property_def.valid_compartments:
            raise ValueError('Compartment {0} does not appear in Set of valid compartments {1}.'
                             .format(self.compartment, ', '.join(str(x) for x in self.property_def.valid_compartments)))

    @property
    def complementary_instances(self):
        return [LocalizationPropertyInstance(self.property_def, compartment) for compartment
                in self.property_def.valid_compartments if compartment != self.compartment]


class Binding:
    @tc.typecheck
    def __init__(self, left_spec: Specification, right_spec: Specification):
        left_spec, right_spec = sorted([left_spec, right_spec])
        self.left_partner = left_spec
        self.right_partner = right_spec
        self._validate()

    @tc.typecheck
    def __eq__(self, other: 'Binding'):
        return self.left_partner == other.left_partner and self.right_partner == other.right_partner

    def __hash__(self):
        return hash(str(self))

    def __lt__(self, other: 'Binding'):
        if self.left_partner != other.left_partner:
            return self.left_partner < other.left_partner
        else:
            return self.right_partner < other.right_partner

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'Binding: L = {0}, R = {1}'.format(self.left_partner, self.right_partner)

    def _validate(self):
        if self.left_partner == self.right_partner:
            raise ValueError('Left - right binding specs are required to be distinct.')


def create_partner_mol_instance(mol_defs: tg.Dict[Specification, MoleculeDefinition], mol_instance: MoleculeInstance) -> MoleculeInstance:
    assert len(mol_instance.association_properties) == 1, \
        'MoleculeInstance passed to create_partner_mol_instance should contain exactly one AssociationPropertyInstance.'

    the_assoc_prop = create_partner_ass_prop_instance(mol_defs, list(mol_instance.association_properties)[0])
    the_mol_def = mol_defs[list(mol_instance.association_properties)[0].partner.to_component_specification()]

    return MoleculeInstance(the_mol_def, set(), {the_assoc_prop}, None)


def create_partner_ass_prop_instance(mol_defs, ass_prop_instance):
    the_other_assoc_spec = ass_prop_instance.property_def.spec

    the_mol_def = mol_defs[ass_prop_instance.partner.to_component_specification()]
    the_assoc_spec = ass_prop_instance.partner
    the_assoc_def = [ass_def for ass_def in the_mol_def.association_defs if ass_def.spec == the_assoc_spec]
    assert len(the_assoc_def) == 1
    the_assoc_def = the_assoc_def[0]

    return AssociationPropertyInstance(the_assoc_def, OccupationStatus.occupied_known_partner, the_other_assoc_spec)