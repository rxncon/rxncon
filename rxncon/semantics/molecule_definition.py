from enum import unique, Enum
from typing import Optional, Set

import typecheck as tc

import rxncon.core.specification as spe


@unique
class Modifier(Enum):
    unmodified     = 'u'
    phosphorylated = 'p'
    ubiquitinated  = 'ub'
    truncated       = 'truncated'

@unique
class OccupationStatus(Enum):
    not_specified = 0
    not_occupied = 1
    occupied_known_partner = 2
    occupied_unknown_partner = 3

@unique
class Compartment(Enum):
    cell = 'cell'
    cytosole = 'cytosole'
    nucleus = 'nucleus'


class Definition:
    pass


class MoleculeDefinition(Definition):
    @tc.typecheck
    def __init__(self, spec: spe.Specification,
                 modification_defs: Set['ModificationDefinition'],
                 association_defs: Set['AssociationDefinition'],
                 localization_def: Optional['LocalizationDefinition']):
        self.spec = spec
        self.modification_defs = modification_defs
        self.association_defs = association_defs

        if localization_def is None:
            localization_def = LocalizationDefinition(set())

        self.localization_def = localization_def

    @tc.typecheck
    def __eq__(self, other: Definition) -> bool:
        return isinstance(other, MoleculeDefinition) and self.spec == other.spec and self.localization_def == other.localization_def and \
            other.modification_defs == self.modification_defs and other.association_defs == self.association_defs

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'MoleculeDefinition: {0}'.format(self.spec)


class ModificationDefinition(Definition):
    @tc.typecheck
    def __init__(self, spec: spe.Specification, valid_modifiers: Set['Modifier']):
        self.spec = spec
        self.valid_modifiers = valid_modifiers

    @tc.typecheck
    def __eq__(self, other: Definition):
        return isinstance(other, ModificationDefinition) and self.spec == other.spec and \
            self.valid_modifiers == other.valid_modifiers

    def __hash__(self) -> int:
        return hash('*mod-def* {0}'.format(self.spec.name))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return 'ModificationDefinition: Domain = {0}, Modifiers = {1}'\
            .format(self.spec, ', '.join(mod.value for mod in self.valid_modifiers))


class AssociationDefinition(Definition):
    @tc.typecheck
    def __init__(self, spec: spe.Specification, valid_partners: Set[spe.Specification]):
        self.spec = spec
        self.valid_partners = valid_partners

    @tc.typecheck
    def __eq__(self, other: Definition) -> bool:
        return isinstance(other, AssociationDefinition) and self.spec == other.spec and \
            self.valid_partners == other.valid_partners

    def __hash__(self) -> int:
        return hash('*ass-def* {}'.format(self.spec.name))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'AssociationDefinition: Domain = {0}, valid_partners = {1}'\
            .format(self.spec, ', '.join(str(x) for x in self.valid_partners))


class LocalizationDefinition(Definition):
    @tc.typecheck
    def __init__(self, valid_compartments: Set[Compartment]):
        self.valid_compartments = valid_compartments

    @tc.typecheck
    def __eq__(self, other: Definition):
        return isinstance(other, LocalizationDefinition) and self.valid_compartments == other.valid_compartments

    def __hash__(self) -> int:
        return hash('*loc-def* with num of compartments {}'.format(len(self.valid_compartments)))

    def __str__(self) -> str:
        return 'LocalizationDefinition: {0}'.format(', '.join(str(x) for x in self.valid_compartments))
