from typing import Set, Optional, Tuple

import typecheck as tc

from rxncon.semantics.molecule_definition import MoleculeDefinition, ModificationDefinition, AssociationDefinition, \
    LocalizationDefinition, Modifier, OccupationStatus, Compartment

import rxncon.core.specification as spe

class Instance:
    pass


class MoleculeInstance(Instance):
    @tc.typecheck
    def __init__(self,
                 molecule_def: MoleculeDefinition,
                 modification_instances: Set['ModificationInstance'],
                 association_instances: Set['AssociationInstance'],
                 localization_instance: Optional['LocalizationInstance']):
        self.molecule_def = molecule_def
        self.modification_instances = modification_instances
        self.association_instances = association_instances
        self.localization_instance = localization_instance

    @tc.typecheck
    def __eq__(self, other: Instance):
        return isinstance(other, MoleculeInstance) and \
            self.molecule_def == other.molecule_def and \
            self.localization_instance == other.localization_instance and \
            self.modification_instances == other.modification_instances and \
            self.association_instances == other.association_instances

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'MoleculeSpecification: {0}, mod_instances = {1}. assoc_instances = {2}. loc_instances = {3}'\
            .format(self.molecule_def.spec, ', '.join([str(x) for x in self.modification_instances]),
                    ', '.join(str(x) for x in self.association_instances), str(self.localization_instance))


class ModificationInstance(Instance):
    @tc.typecheck
    def __init__(self, modification_def: ModificationDefinition, modifier: 'Modifier'):
        self.modification_def = modification_def
        self.modifier = modifier
        self._validate()

    @tc.typecheck
    def __eq__(self, other: Instance) -> bool:
        return isinstance(other, ModificationInstance) and self.modification_def == other.modification_def and \
            self.modifier == other.modifier

    def __hash__(self) -> bool:
        return hash(str(self))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'ModificationSpecification: Domain = {0}, Modifier = {1}'\
            .format(self.modification_def.spec, self.modifier)

    def _validate(self):
        if self.modifier not in self.modification_def.valid_modifiers:
            raise ValueError('Modifier {0} does not appear in Set of valid modifiers for domain {1}.'
                             .format(self.modifier, self.modification_def.spec))

    def complementary_instances(self):
        return [ModificationInstance(self.modification_def, modifier) for modifier
                in self.modification_def.valid_modifiers if modifier != self.modifier]


class AssociationInstance(Instance):
    @tc.typecheck
    def __init__(self, association_def: AssociationDefinition, occupation_status: 'OccupationStatus',
                 partner: Optional[spe.Specification]):
        self.association_def = association_def
        self.occupation_status = occupation_status
        self.partner = partner

        self._validate()

    @tc.typecheck
    def __eq__(self, other: Instance) -> bool:
        return isinstance(other, AssociationInstance) and self.association_def == other.association_def and \
            self.occupation_status == other.occupation_status and self.partner == other.partner

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'AssociationInstance: Domain = {0}, occupation_status = {1}'\
            .format(self.association_def.spec, self.occupation_status)

    def complementary_instances(self):
        # todo: is this correct?
        if self.occupation_status == OccupationStatus.occupied_known_partner:
            return [AssociationInstance(self.association_def, OccupationStatus.not_occupied, self.partner)]
        elif self.occupation_status == OccupationStatus.not_occupied:
            return [AssociationInstance(self.association_def, OccupationStatus.occupied_known_partner, self.partner)]
        else:
            raise NotImplementedError

    def _validate(self):
        # For associations the molecule/domain/subdomain spec should match exactly.
        assert self.partner in self.association_def.valid_partners


class LocalizationInstance(Instance):
    @tc.typecheck
    def __init__(self, localization_def: LocalizationDefinition, compartment: Compartment):
        self.localization_def = localization_def
        self.compartment = compartment
        self._validate()

    @tc.typecheck
    def __eq__(self, other: Instance) -> bool:
        return isinstance(other, LocalizationInstance) and self.localization_def == other.localization_def and \
            self.compartment == other.compartment

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'LocalizationSpecification: {0}'.format(self.compartment)

    def _validate(self):
        if self.compartment not in self.localization_def.valid_compartments:
            raise ValueError('Compartment {0} does not appear in Set of valid compartments {1}.'
                             .format(self.compartment, ', '.join(str(x) for x in self.localization_def.valid_compartments)))

    def complementary_instances(self):
        return [LocalizationInstance(self.localization_def, compartment) for compartment
                in self.localization_def.valid_compartments if compartment != self.compartment]


class Binding:
    @tc.typecheck
    def __init__(self, left_partner: Tuple[int, AssociationInstance], right_partner: Tuple[int, AssociationInstance]):
        self.left_partner = left_partner
        self.right_partner = right_partner
        self._validate()

    @tc.typecheck
    def __eq__(self, other: 'Binding'):
        return self.left_partner == other.left_partner and self.right_partner == other.right_partner

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'Binding: L_molecule_index = {0}, L_domain = {1}, R_molecule_index = {2}, R_domain = {3}'\
            .format(self.left_partner[0], self.left_partner[1].association_def.spec,
                    self.right_partner[0], self.right_partner[1].association_def.spec)

    def _validate(self):
        if not self.left_partner[1].occupation_status or not self.right_partner[1].occupation_status:
            raise ValueError('Binding requires both partners to have occupied association domains.')

        if self.left_partner[0] == self.right_partner[0]:
            raise ValueError('Binding-molecule-indices are required to be distinct for each binding.')