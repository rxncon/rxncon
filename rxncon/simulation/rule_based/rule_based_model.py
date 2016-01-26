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

        for rule in self.rules:
            for molecule in rule.molecules:
                if molecule not in self.molecule_definitions:
                    raise ValueError('Rule {0} contains molecule definition {1}, which is absent in the model'
                                     .format(rule, molecule))


class MoleculeDefinition:
    def __init__(self, name: str, modification_definitions: Optional[List['ModificationDefinition']],
                 association_definitions: Optional[List['AssociationDefinition']],
                 localization_definition: Optional['LocalizationDefinition']):
        self.name = name
        self.modification_definitions = modification_definitions
        self.association_definitions = association_definitions
        self.localization_definition = localization_definition

    def __str__(self):
        return 'MoleculeDefinition: {0}'.format(self.name)


class MoleculeSpecification:
    def __init__(self, molecule_definition: MoleculeDefinition,
                 modification_specifications: List['ModificationSpecification'],
                 association_specifications: List['AssociationSpecification'],
                 localization_specification: Optional['LocalizationSpecification']):
        self.molecule_definition = molecule_definition
        self.modification_specifications = modification_specifications
        self.association_specifications = association_specifications
        self.localization_specification = localization_specification

    def __str__(self):
        return 'MoleculeSpecification: {0}, mod_specs = {1}. ass_specs = {2}. loc_spec = {3}'\
            .format(self.molecule_definition.name, ', '.join(self.modification_specifications),
                    ', '.join(self.association_specifications), self.localization_specification)


class ModificationDefinition:
    def __init__(self, domain: str, valid_modifiers: List[str]):
        self.domain = domain
        self.valid_modifiers = valid_modifiers
        self._validate()

    def __str__(self):
        return 'ModificationDefinition: Domain = {0}, Modifiers = {1}'\
            .format(self.domain, ', '.join(self.valid_modifiers))

    def _validate(self):
        if len(self.valid_modifiers) > len(set(self.valid_modifiers)):
            modifiers = ', '.join(self.valid_modifiers)
            raise ValueError('Modifier list {0} for domain {1} contains non-unique elements.'
                             .format(modifiers, self.domain))


class ModificationSpecification:
    def __init__(self, modification_definition: ModificationDefinition, modifier: str):
        self.modification_definition = modification_definition
        self.modifier = modifier
        self._validate()

    def __str__(self):
        return 'ModificationSpecification: Domain = {0}, Modifier = {1}'\
            .format(self.modification_definition.domain, self.modifier)

    def _validate(self):
        if self.modifier not in self.modification_definition.valid_modifiers:
            raise ValueError('Modifier {0} does not appear in list of valid modifiers for domain {1}.'
                             .format(self.modifier, self.modification_definition.domain))


class AssociationDefinition:
    def __init__(self, domain: str):
        self.domain = domain

    def __str__(self):
        return 'AssociationDefinition: Domain = {0}'.format(self.domain)


class AssociationSpecification:
    def __init__(self, association_definition: AssociationDefinition, is_occupied: bool):
        self.association_definition = association_definition
        self.is_occupied = is_occupied

    def __str__(self):
        return 'AssociationSpecification: Domain = {0}, occupied = {1}'\
            .format(self.association_definition.domain, self.is_occupied)


class LocalizationDefinition:
    def __init__(self, compartments: List[str]):
        self.compartments = compartments

    def __str__(self):
        return 'LocalizationDefinition: {0}'.format(', '.join(self.compartments))


class LocalizationSpecification:
    def __init__(self, localization_definition: LocalizationDefinition, compartment: str):
        self.localization_definition = localization_definition
        self.compartment = compartment
        self._validate()

    def __str__(self):
        return 'LocalizationSpecification: {0}'.format(self.compartment)

    def _validate(self):
        if self.compartment not in self.localization_definition.compartments:
            raise ValueError('Compartment {0} does not appear in list of valid compartments {1}.'
                             .format(self.compartment, ', '.join(self.localization_definition.compartments)))


class Rule:
    def __init__(self, left_hand_side: List['Reactant'], right_hand_side: List['Reactant'], arrow_type: 'Arrow'):
        self.left_hand_side = left_hand_side
        self.right_hand_side = right_hand_side
        self.arrow_type = arrow_type
        self._validate()

    @property
    def molecules(self):
        molecules = []
        for side in [self.left_hand_side, self.right_hand_side]:
            if isinstance(side, MoleculeReactant):
                molecules.append(side.molecule_specification.molecule_definition)

            elif isinstance(side, ComplexReactant):
                molecules += side.molecules

        return molecules

    def _validate(self):
        pass


class Reactant:
    pass


class MoleculeReactant(Reactant):
    def __init__(self, molecule_specification: MoleculeSpecification):
        self.molecule_specification = molecule_specification

    def __str__(self):
        return 'MoleculeReactant: [{0}]'.format(self.molecule_specification)


class ComplexReactant(Reactant):
    def __init__(self, molecules: List[MoleculeSpecification], bindings: List['Binding']):
        self.molecules = molecules
        self.bindings = bindings
        self._validate()

    def __str__(self):
        return 'ComplexReactant: Molecules = [{0}], Bindings = [{1}]'\
            .format(', '.join(self.molecules), ', '.join(self.bindings))

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
            .format(self.left_partner[0], self.left_partner[1].association_definition.domain,
                    self.right_partner[0], self.right_partner[1].association_definition.domain)

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
