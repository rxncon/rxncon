import typing as tg

import typecheck as tc

import rxncon.core.specification as spe
from rxncon.core.specification import Specification
from rxncon.semantics.molecule_definition import MoleculeDefinition, ModificationPropertyDefinition, AssociationPropertyDefinition, \
    LocalizationPropertyDefinition, Modifier, OccupationStatus, Compartment


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

    @property
    def bindings(self) -> tg.Set['Binding']:
        return {Binding(ass_prop.association_def.spec, ass_prop.partner) for ass_prop in self.association_properties
                if ass_prop.occupation_status == OccupationStatus.occupied_known_partner}

    def _validate(self):
        # Assert each modification domain and each association domain is only present once.
        assert len([mod_prop.modification_def for mod_prop in self.modification_properties]) == \
            len(set([mod_prop.modification_def for mod_prop in self.modification_properties]))

        assert len([ass_prop.association_def for ass_prop in self.association_properties]) == \
            len(set([ass_prop.association_def for ass_prop in self.association_properties]))


class PropertyInstance:
    pass


class ModificationPropertyInstance(PropertyInstance):
    @tc.typecheck
    def __init__(self, modification_def: ModificationPropertyDefinition, modifier: 'Modifier'):
        self.modification_def = modification_def
        self.modifier = modifier
        self._validate()

    @tc.typecheck
    def __eq__(self, other: PropertyInstance) -> bool:
        return isinstance(other, ModificationPropertyInstance) and self.modification_def == other.modification_def and \
            self.modifier == other.modifier

    def __lt__(self, other: 'ModificationPropertyInstance'):
        if self.modification_def.spec < other.modification_def.spec:
            return True

        return self.modifier.value < other.modifier.value

    def __hash__(self) -> bool:
        return hash(str(self.modification_def.spec))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'mod/{0}:{1}'.format(self.modification_def.spec, self.modifier)

    def _validate(self):
        if self.modifier not in self.modification_def.valid_modifiers:
            raise ValueError('Modifier {0} does not appear in Set of valid modifiers for domain {1}.'
                             .format(self.modifier, self.modification_def.spec))

    def complementary_instances(self):
        return [ModificationPropertyInstance(self.modification_def, modifier) for modifier
                in self.modification_def.valid_modifiers if modifier != self.modifier]


class AssociationPropertyInstance(PropertyInstance):
    @tc.typecheck
    def __init__(self, association_def: AssociationPropertyDefinition, occupation_status: 'OccupationStatus',
                 partner: tg.Optional[spe.Specification]):
        self.association_def = association_def
        self.occupation_status = occupation_status
        self.partner = partner

        self._validate()

    @tc.typecheck
    def __eq__(self, other: PropertyInstance) -> bool:
        if not isinstance(other, AssociationPropertyInstance):
            return False

        if self.partner and other.partner:
            return self.association_def == other.association_def and self.occupation_status == other.occupation_status and self.partner == other.partner
        elif not self.partner and not other.partner:
            return self.association_def == other.association_def and self.occupation_status == other.occupation_status
        else:
            return False

    def __lt__(self, other: 'AssociationPropertyInstance'):
        if self.association_def.spec < other.association_def.spec:
            return True
        return False

    def __hash__(self) -> int:
        return hash(str(self.association_def.spec))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'ass/{0}:{1}'.format(self.association_def.spec, self.partner)

    def complementary_instances(self):
        if self.occupation_status == OccupationStatus.occupied_known_partner:
            unoccupied = AssociationPropertyInstance(self.association_def, OccupationStatus.not_occupied, None)
            other_partners = [AssociationPropertyInstance(self.association_def, OccupationStatus.occupied_known_partner, x)
                              for x in self.association_def.valid_partners if x != self.partner]

            return other_partners + [unoccupied]
        elif self.occupation_status == OccupationStatus.not_occupied:
            return [AssociationPropertyInstance(self.association_def, OccupationStatus.occupied_known_partner, x)
                    for x in self.association_def.valid_partners]
        else:
            raise NotImplementedError

    def _validate(self):
        if self.partner:
            assert self.partner in self.association_def.valid_partners, \
                'For association domain {0}, binding partner {1} does not appear in list of valid partners {2}.' \
                .format(str(self.association_def.spec), str(self.partner), ','.join(str(x) for x in sorted(self.association_def.valid_partners)))

        if self.occupation_status == OccupationStatus.not_occupied:
            assert not self.partner


class LocalizationPropertyInstance(PropertyInstance):
    @tc.typecheck
    def __init__(self, localization_def: LocalizationPropertyDefinition, compartment: Compartment):
        self.localization_def = localization_def
        self.compartment = compartment
        self._validate()

    @tc.typecheck
    def __eq__(self, other: PropertyInstance) -> bool:
        return isinstance(other, LocalizationPropertyInstance) and self.localization_def == other.localization_def and \
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
        if self.compartment not in self.localization_def.valid_compartments:
            raise ValueError('Compartment {0} does not appear in Set of valid compartments {1}.'
                             .format(self.compartment, ', '.join(str(x) for x in self.localization_def.valid_compartments)))

    def complementary_instances(self):
        return [LocalizationPropertyInstance(self.localization_def, compartment) for compartment
                in self.localization_def.valid_compartments if compartment != self.compartment]


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

    the_other_assoc_spec = list(mol_instance.association_properties)[0].association_def.spec

    the_mol_def = mol_defs[list(mol_instance.association_properties)[0].partner.to_component_specification()]
    the_assoc_spec = list(mol_instance.association_properties)[0].partner
    the_assoc_def = [ass_def for ass_def in the_mol_def.association_defs if ass_def.spec == the_assoc_spec]
    assert len(the_assoc_def) == 1
    the_assoc_def = the_assoc_def[0]

    the_assoc_prop = AssociationPropertyInstance(the_assoc_def, OccupationStatus.occupied_known_partner, the_other_assoc_spec)

    return MoleculeInstance(the_mol_def, set(), {the_assoc_prop}, None)
