from typing import Optional, List, Tuple



class MoleculeDefinition:
    def __init__(self, name: str, modification_definitions: Optional[List['ModificationDefinition']],
                 association_definitions: Optional[List['AssociationDefinition']],
                 localization_definitions: Optional[List['LocalizationDefinition']]):
        self.name = name
        self.modification_definitions = modification_definitions
        self.association_definitions = association_definitions
        self.localization_definitions = localization_definitions


class MoleculeSpecification:
    def __init__(self, molecule_definition: MoleculeDefinition, modification_specifications: List['ModificationSpecification'],
                 association_specifications: List['AssociationSpecification'],
                 localization_specifications: List['LocalizationSpecification']):
        self.molecule_definition = molecule_definition
        self.modification_specifications = modification_specifications
        self.association_specifications = association_specifications
        self.localization_specifications = localization_specifications


class ModificationDefinition:
    def __init__(self, domain_name: str, valid_modifiers: List[str]):
        self.domain_name = domain_name
        self.valid_modifiers = valid_modifiers


class ModificationSpecification:
    def __init__(self, modification_definition: ModificationDefinition, value: str):
        self.modification_definition = modification_definition
        self.value = value


class AssociationDefinition:
    def __init__(self, domain_name: str):
        self.domain_name = domain_name


class AssociationSpecification:
    def __init__(self, association_definition: AssociationDefinition, is_occupied: bool):
        self.association_definition = association_definition
        self.is_occupied = is_occupied


class LocalizationDefinition:
    def __init__(self, compartment: str):
        self.compartment = compartment


class LocalizationSpecification:
    def __init__(self, localization_definition: LocalizationDefinition, is_localized: bool):
        self.localization_definition = localization_definition
        self.is_localized = is_localized


class Rule:
    def __init__(self, left_hand_side: List['Reactant'], right_hand_side: List['Reactant'], arrow_type: 'ArrowType'):
        self.left_hand_side = left_hand_side
        self.right_hand_side = right_hand_side
        self.arrow_type = arrow_type


class Reactant:
    pass


class MoleculeReactant(Reactant):
    def __init__(self, molecule_specification: MoleculeSpecification):
        self.molecule_specification = molecule_specification


class ComplexReactant(Reactant):
    def __init__(self, complex_parts: List[MoleculeSpecification], complex_bindings: List['Binding']):
        self.complex_parts = complex_parts
        self.complex_bindings = complex_bindings


class Binding:
    def __init__(self, left_partner: Tuple[int, AssociationSpecification], right_partner: Tuple[int, AssociationSpecification]):
        self.left_partner = left_partner
        self.right_partner = right_partner



