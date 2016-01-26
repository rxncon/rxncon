from typing import Optional, List, Tuple
from enum import Enum


class RuleBasedModel:
    def __init__(self, molecule_definitions: List['MoleculeDefinition'], rules: List['Rule'],
                 parameters: List['Parameter'], initial_conditions: List['InitialCondition']):
        self.molecule_definitions = molecule_definitions
        self.rules = rules
        self.parameters = parameters
        self.initial_conditions = initial_conditions

        self._validate()

    def _validate(self):
        for initial_condition in self.initial_conditions:
            if initial_condition.molecule_specification.molecule_definition not in self.molecule_definitions:
                raise ValueError('Initial condition {0} refers to unknown molecule definition {1}.'
                                 .format(initial_condition, initial_condition.molecule_specification.molecule_definition))


class MoleculeDefinition:
    def __init__(self, name: str, modification_definitions: Optional[List['ModificationDefinition']],
                 association_definitions: Optional[List['AssociationDefinition']],
                 localization_definition: Optional['LocalizationDefinition']):
        self.name = name

        self.modification_definitions = modification_definitions
        self.association_definitions = association_definitions

        self.localization_definition = localization_definition


class MoleculeSpecification:
    def __init__(self, molecule_definition: MoleculeDefinition,
                 modification_specifications: List['ModificationSpecification'],
                 association_specifications: List['AssociationSpecification'],
                 localization_specification: Optional['LocalizationSpecification']):
        self.molecule_definition = molecule_definition

        self.modification_specifications = modification_specifications
        self.association_specifications = association_specifications

        self.localization_specification = localization_specification


class ModificationDefinition:
    def __init__(self, domain_name: str, valid_modifiers: List[str]):
        self.domain_name = domain_name
        self.valid_modifiers = valid_modifiers
        self._validate()

    def _validate(self):
        if len(self.valid_modifiers) > len(set(self.valid_modifiers)):
            modifiers = ', '.join(self.valid_modifiers)
            raise ValueError('Modifier list {0} for domain {1} contains non-unique elements.'
                             .format(modifiers, self.domain_name))


class ModificationSpecification:
    def __init__(self, modification_definition: ModificationDefinition, value: str):
        self.modification_definition = modification_definition
        self.value = value
        self._validate()

    def _validate(self):
        if self.value not in self.modification_definition.valid_modifiers:
            raise ValueError('Modifier {0} does not appear in list of valid modifiers for domain {1}.'
                             .format(self.value, self.modification_definition.domain_name))


class AssociationDefinition:
    def __init__(self, domain_name: str):
        self.domain_name = domain_name


class AssociationSpecification:
    def __init__(self, association_definition: AssociationDefinition, is_occupied: bool):
        self.association_definition = association_definition
        self.is_occupied = is_occupied


class LocalizationDefinition:
    def __init__(self, compartments: List[str]):
        self.compartments = compartments


class LocalizationSpecification:
    def __init__(self, localization_definition: LocalizationDefinition, current_compartment: str):
        self.localization_definition = localization_definition
        self.current_compartment = current_compartment
        self._validate()

    def _validate(self):
        if self.current_compartment not in self.localization_definition.compartments:
            raise ValueError('Compartment {0} does not appear in list of valid compartments {1}.'
                             .format(self.current_compartment, ', '.join(self.localization_definition.compartments)))


class Rule:
    def __init__(self, left_hand_side: List['Reactant'], right_hand_side: List['Reactant'], arrow_type: 'Arrow'):
        self.left_hand_side = left_hand_side
        self.right_hand_side = right_hand_side
        self.arrow_type = arrow_type
        self._validate()

    def _validate(self):
        pass


class Reactant:
    pass


class MoleculeReactant(Reactant):
    def __init__(self, molecule_specification: MoleculeSpecification):
        self.molecule_specification = molecule_specification


class ComplexReactant(Reactant):
    def __init__(self, molecules: List[MoleculeSpecification], bindings: List['Binding']):
        self.molecules = molecules
        self.bindings = bindings
        self._validate()

    def _validate(self):
        localizations = [molecule.localization_specification for molecule in self.molecules]
        if len(set(localizations)) > 1:
            raise ValueError('Molecules making up a ComplexReactant cannot be in different localizations: {0}.'
                             .format(', '.join(localizations)))


class Binding:
    def __init__(self, left_partner: Tuple[int, AssociationSpecification], right_partner: Tuple[int, AssociationSpecification]):
        self.left_partner = left_partner
        self.right_partner = right_partner
        self._validate()

    def __str__(self):
        return 'Binding: L_index = {0}, L_domain = {1}, R_index = {2}, R_domain = {3}'\
            .format(self.left_partner[0], self.left_partner[1].association_definition.domain_name,
                    self.right_partner[0], self.right_partner[1].association_definition.domain_name)

    def _validate(self):
        if not self.left_partner[1].is_occupied or not self.right_partner[1].is_occupied:
            raise ValueError('Binding requires both partners to have occupied association domains.')

        if self.left_partner[0] == self.right_partner[0]:
            raise ValueError('Binding-molecule-indices are required to be distinct for each binding.')


class Arrow(Enum):
    irreversible = '->'
    reversible   = '<->'


class Parameter:
    def __init__(self, name: str, value):
        self.name = name
        self.value = value

    def __str__(self):
        return 'Parameter: {0} = {1}'.format(self.name, self.value)


class InitialCondition:
    def __init__(self, molecule_specification: MoleculeSpecification, value):
        self.molecule_specification = molecule_specification
        self.value = value

    def __str__(self):
        return 'InitialCondition: {0} = {1}'.format(self.molecule_specification, self.value)
