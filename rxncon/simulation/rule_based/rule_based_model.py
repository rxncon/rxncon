from typing import Optional, List, Tuple
from enum import Enum, unique
import typecheck as tc


class RuleBasedModel:
    @tc.typecheck
    def __init__(self, molecule_defs: List['MoleculeDefinition'], rules: List['Rule'],
                 parameters: List['Parameter'], initial_conditions: List['InitialCondition']):
        self.molecule_defs = molecule_defs
        self.rules = rules
        self.parameters = parameters
        self.initial_conditions = initial_conditions

        self._validate()

    def _validate(self):
        for initial_condition in self.initial_conditions:
            if initial_condition.molecule_specification.molecule_def not in self.molecule_defs:
                raise ValueError('Initial condition {0} refers to unknown molecule def {1}.'
                                 .format(initial_condition, initial_condition.molecule_specification.molecule_def))

        for rule in self.rules:
            for molecule in rule.molecules:
                if molecule not in self.molecule_defs:
                    raise ValueError('Rule {0} contains molecule def {1}, which is absent in the model'
                                     .format(rule, molecule))


class MoleculeDefinition:
    @tc.typecheck
    def __init__(self, name: str, modification_defs: Optional[List['ModificationDefinition']],
                 association_defs: Optional[List['AssociationDefinition']],
                 localization_def: Optional['LocalizationDefinition']):
        self.name = name
        self.modification_defs = modification_defs
        self.association_defs = association_defs
        self.localization_def = localization_def

    @tc.typecheck
    def __eq__(self, other: 'MoleculeDefinition') -> bool:
        return self.name == other.name and self.localization_def == other.localization_def and \
            other.modification_defs == self.modification_defs and other.association_defs == self.association_defs

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return 'MoleculeDefinition: {0}'.format(self.name)


class MoleculeSpecification:
    @tc.typecheck
    def __init__(self, molecule_def: MoleculeDefinition,
                 modification_specs: List['ModificationSpecification'],
                 association_specs: List['AssociationSpecification'],
                 localization_spec: Optional['LocalizationSpecification']):
        self.molecule_def = molecule_def
        self.modification_specs = modification_specs
        self.association_specs = association_specs
        self.localization_spec = localization_spec

    @tc.typecheck
    def __eq__(self, other: 'MoleculeSpecification'):
        assert isinstance(other, MoleculeSpecification)
        return self.molecule_def == other.molecule_def and self.localization_spec == other.localization_spec and \
            other.modification_specs == self.modification_specs and other.association_specs == other.association_specs

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return 'MoleculeSpecification: {0}, mod_specs = {1}. ass_specs = {2}. loc_spec = {3}'\
            .format(self.molecule_def.name, ', '.join([str(x) for x in self.modification_specs]),
                    ', '.join(str(x) for x in self.association_specs), str(self.localization_spec))


class ModificationDefinition:
    @tc.typecheck
    def __init__(self, domain: str, valid_modifiers: List[str]):
        self.domain = domain
        self.valid_modifiers = valid_modifiers
        self._validate()

    @tc.typecheck
    def __eq__(self, other: 'ModificationDefinition'):
        return self.domain == other.domain and self.valid_modifiers == other.valid_modifiers

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return 'ModificationDefinition: Domain = {0}, Modifiers = {1}'\
            .format(self.domain, ', '.join(self.valid_modifiers))

    def _validate(self):
        if len(self.valid_modifiers) > len(set(self.valid_modifiers)):
            modifiers = ', '.join(self.valid_modifiers)
            raise ValueError('Valid modifier list {0} for domain {1} contains non-unique elements.'
                             .format(modifiers, self.domain))


class ModificationSpecification:
    @tc.typecheck
    def __init__(self, modification_def: ModificationDefinition, modifier: str):
        self.modification_def = modification_def
        self.modifier = modifier
        self._validate()

    @tc.typecheck
    def __eq__(self, other: 'ModificationSpecification') -> bool:
        return self.modification_def == other.modification_def and self.modifier == other.modifier

    def __hash__(self) -> bool:
        return hash(str(self))

    def __str__(self) -> str:
        return 'ModificationSpecification: Domain = {0}, Modifier = {1}'\
            .format(self.modification_def.domain, self.modifier)

    def _validate(self):
        if self.modifier not in self.modification_def.valid_modifiers:
            raise ValueError('Modifier {0} does not appear in list of valid modifiers for domain {1}.'
                             .format(self.modifier, self.modification_def.domain))


class AssociationDefinition:
    @tc.typecheck
    def __init__(self, domain: str):
        self.domain = domain

    @tc.typecheck
    def __eq__(self, other: 'AssociationDefinition') -> bool:
        return self.domain == other.domain

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return 'AssociationDefinition: Domain = {0}'.format(self.domain)


class AssociationSpecification:
    @tc.typecheck
    def __init__(self, association_def: AssociationDefinition, occupation_status: 'OccupationStatus'):
        self.association_def = association_def
        self.occupation_status = occupation_status

    @tc.typecheck
    def __eq__(self, other: 'AssociationSpecification') -> bool:
        return self.association_def == other.association_def and self.occupation_status == other.occupation_status

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return 'AssociationSpecification: Domain = {0}, occupation_status = {1}'\
            .format(self.association_def.domain, self.occupation_status)


@unique
class OccupationStatus(Enum):
    not_specified = 0
    not_unoccupied = 1
    occupied_known_partner = 2
    occupied_unknown_partner = 3


class LocalizationDefinition:
    @tc.typecheck
    def __init__(self, valid_compartments: List[str]):
        self.valid_compartments = valid_compartments

    @tc.typecheck
    def __eq__(self, other: 'LocalizationDefinition'):
        return self.valid_compartments == other.valid_compartments

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return 'LocalizationDefinition: {0}'.format(', '.join(self.valid_compartments))


class LocalizationSpecification:
    @tc.typecheck
    def __init__(self, localization_def: LocalizationDefinition, compartment: str):
        self.localization_def = localization_def
        self.compartment = compartment
        self._validate()

    @tc.typecheck
    def __eq__(self, other: 'LocalizationSpecification') -> bool:
        return self.localization_def == other.localization_def and self.compartment == other.compartment

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return 'LocalizationSpecification: {0}'.format(self.compartment)

    def _validate(self):
        if self.compartment not in self.localization_def.valid_compartments:
            raise ValueError('Compartment {0} does not appear in list of valid compartments {1}.'
                             .format(self.compartment, ', '.join(self.localization_def.valid_compartments)))


class Rule:
    @tc.typecheck
    def __init__(self, left_hand_side: List['Reactant'], right_hand_side: List['Reactant'], arrow_type: 'Arrow', rates: List['Parameter']):
        self.left_hand_side = left_hand_side
        self.right_hand_side = right_hand_side
        self.arrow_type = arrow_type
        self.rates = rates
        self._validate()

    @tc.typecheck
    def __eq__(self, other: 'Rule'):
        return self.left_hand_side == other.left_hand_side and self.right_hand_side == other.right_hand_side and \
            self.arrow_type == other.arrow_type and self.rates == other.rates

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return 'Rule: {0} {1} {2}, {3}'.format('+'.join(str(x) for x in self.left_hand_side), self.arrow_type,
                                               '+'.join(str(x) for x in self.right_hand_side), ', '.join(str(x) for x in self.rates))

    @property
    def molecules(self):
        molecules = []
        for side in [self.left_hand_side, self.right_hand_side]:
            if isinstance(side, MoleculeReactant) and side.molecule_specification.molecule_definition not in molecules:
                molecules.append(side.molecule_specification.molecule_def)

            elif isinstance(side, ComplexReactant):
                [molecules.append(x) for x in side.molecules if x not in molecules]

        return molecules

    def _validate(self):
        if self.arrow_type == Arrow.irreversible and len(self.rates) != 1:
            raise ValueError('Rule {0} is irreversible and thus requires exactly one rate constant, {1} given'
                             .format(str(self), len(self.rates)))

        if self.arrow_type == Arrow.reversible and len(self.rates) != 2:
            raise ValueError('Rule {0} is reversible and thus requires exactly two rate constants, {1} given'
                             .format(str(self), len(self.rates)))


class Reactant:
    pass


class MoleculeReactant(Reactant):
    @tc.typecheck
    def __init__(self, molecule_specification: MoleculeSpecification):
        self.molecule_specification = molecule_specification

    @tc.typecheck
    def __eq__(self, other: 'Reactant') -> bool:
        return isinstance(other, MoleculeReactant) and self.molecule_specification == other.molecule_specification

    def __str__(self) -> str:
        return 'MoleculeReactant: [{0}]'.format(self.molecule_specification)


class ComplexReactant(Reactant):
    @tc.typecheck
    def __init__(self, molecules: List[MoleculeSpecification], bindings: List['Binding']):
        self.molecules = molecules
        self.bindings = bindings
        self._validate()

    @tc.typecheck
    def __eq__(self, other: 'Reactant'):
        return isinstance(other, ComplexReactant) and self.molecules == other.molecules and self.bindings == other.bindings

    def __str__(self) -> str:
        return 'ComplexReactant: Molecules = [{0}], Bindings = [{1}]'\
            .format(', '.join(str(x) for x in self.molecules), ', '.join(str(x) for x in self.bindings))

    def _validate(self):
        unique_localizations = {molecule.localization_spec for molecule in self.molecules}
        if len(unique_localizations) > 1:
            raise ValueError('Molecules making up a ComplexReactant cannot be in different localizations: {0}.'
                             .format(', '.join(str(x) for x in unique_localizations)))


class Binding:
    @tc.typecheck
    def __init__(self, left_partner: Tuple[int, AssociationSpecification], right_partner: Tuple[int, AssociationSpecification]):
        self.left_partner = left_partner
        self.right_partner = right_partner
        self._validate()

    @tc.typecheck
    def __eq__(self, other: 'Binding'):
        return self.left_partner == other.left_partner and self.right_partner == other.right_partner

    def __str__(self) -> str:
        return 'Binding: L_molecule_index = {0}, L_domain = {1}, R_molecule_index = {2}, R_domain = {3}'\
            .format(self.left_partner[0], self.left_partner[1].association_def.domain,
                    self.right_partner[0], self.right_partner[1].association_def.domain)

    def _validate(self):
        if not self.left_partner[1].occupation_status or not self.right_partner[1].occupation_status:
            raise ValueError('Binding requires both partners to have occupied association domains.')

        if self.left_partner[0] == self.right_partner[0]:
            raise ValueError('Binding-molecule-indices are required to be distinct for each binding.')


class Arrow(Enum):
    irreversible = '->'
    reversible   = '<->'


class Parameter:
    @tc.typecheck
    def __init__(self, name: str, value: str):
        self.name = name
        self.value = value

    @tc.typecheck
    def __eq__(self, other: 'Parameter') -> bool:
        assert isinstance(other, Parameter)
        return self.name == other.name and self.value == other.value

    def __str__(self) -> bool:
        return 'Parameter: {0} = {1}'.format(self.name, self.value)


class InitialCondition:
    def __init__(self, molecule_specification: MoleculeSpecification, value):
        self.molecule_specification = molecule_specification
        self.value = value

    def __eq__(self, other):
        assert isinstance(other, InitialCondition)
        return self.molecule_specification == other.molecule_specification and self.value == other.value

    def __str__(self):
        return 'InitialCondition: {0} = {1}'.format(self.molecule_specification, self.value)
