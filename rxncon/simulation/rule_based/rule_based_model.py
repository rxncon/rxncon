from typing import Optional, List, Tuple
from enum import Enum


class RuleBasedModel:
    def __init__(self, molecule_definitions: List['MoleculeDefinition'], rules: List['Rule'],
                 parameters: List['Parameter'], initial_conditions: List['InitialCondition']):
        self.molecule_definitions = molecule_definitions
        self.rules = rules
        self.parameters = parameters
        self.initial_conditions = initial_conditions

        self.validate()

    def validate(self):
        [x.validate() for x in self.rules]
        [x.validate() for x in self.initial_conditions]


class MoleculeDefinition:
    def __init__(self, name: str, modification_definitions: Optional[List['ModificationDefinition']],
                 association_definitions: Optional[List['AssociationDefinition']],
                 localization_definitions: Optional[List['LocalizationDefinition']]):
        self.name = name
        self.modification_definitions = modification_definitions
        self.association_definitions = association_definitions
        self.localization_definitions = localization_definitions


class MoleculeSpecification:
    def __init__(self, molecule_definition: MoleculeDefinition,
                 modification_specifications: List['ModificationSpecification'],
                 association_specifications: List['AssociationSpecification'],
                 localization_specifications: List['LocalizationSpecification']):
        self.molecule_definition = molecule_definition
        self.modification_specifications = modification_specifications
        self.association_specifications = association_specifications
        self.localization_specifications = localization_specifications
        self.validate()

    def validate(self):

        def validate_unique_specification_property():
            for specification_name in ["modification", "association"]:
                specification_property_list = [getattr(specification, specification_name +"_definition").domain_name for specification in getattr(self , specification_name + "_specifications")]
                assert all(specification_property_list.count(domain) == 1 for domain in specification_property_list)
        for modification_specification in self.modification_specifications:
            if modification_specification.modification_definition in self.molecule_definition.modification_definitions:
                pass
        assert all(modification_specification.modification_definition in self.molecule_definition.modification_definitions
                   for modification_specification in self.modification_specifications)

        assert all(association_specification.association_definition in self.molecule_definition.association_definitions
                   for association_specification in self.association_specifications)

        assert all(localization_specification.localization_definition in self.molecule_definition.localization_definitions
                   for localization_specification in self.localization_specifications)

        validate_unique_specification_property()

        # a molecule can only be localized at one place at a time
        assert len(self.localization_specifications) == 1
        # if a localization is defined it should be localized there in case of specification
        assert self.localization_specifications[0].is_localized


class ModificationDefinition:
    def __init__(self, domain_name: str, valid_modifiers: List[str]):
        self.domain_name = domain_name
        self.valid_modifiers = valid_modifiers


class ModificationSpecification:
    def __init__(self, modification_definition: ModificationDefinition, value: str):
        self.modification_definition = modification_definition
        self.value = value
        self.validate()

    def validate(self):
        assert self.value in self.modification_definition.valid_modifiers


class AssociationDefinition:
    def __init__(self, domain_name: str):
        self.domain_name = domain_name


class AssociationSpecification:
    def __init__(self, association_definition: AssociationDefinition, is_occupied: bool):
        self.association_definition = association_definition
        self.is_occupied = is_occupied

    def validate(self):
        assert isinstance(self.association_definition, AssociationDefinition)
        assert isinstance(self.is_occupied, bool)


class LocalizationDefinition:
    def __init__(self, compartment: str):
        self.compartment = compartment


class LocalizationSpecification:
    def __init__(self, localization_definition: LocalizationDefinition, is_localized: bool):
        self.localization_definition = localization_definition
        self.is_localized = is_localized
        self.validate()

    def validate(self):
        assert isinstance(self.is_localized, bool)


class Rule:
    def __init__(self, left_hand_side: List['Reactant'], right_hand_side: List['Reactant'], arrow_type: 'Arrow'):
        self.left_hand_side = left_hand_side
        self.right_hand_side = right_hand_side
        self.arrow_type = arrow_type
        self.validate()

    def validate(self):
        def localization_validation(hand_side: List[Reactant]):
            reactant_localization = set()
            for reactant in hand_side:
                if isinstance(reactant, MoleculeReactant):
                    if reactant.molecule_specification.localization_specifications[0].is_localized:
                        reactant_localization.add(reactant.molecule_specification.localization_specifications[0].localization_definition.compartment)
                elif isinstance(reactant, ComplexReactant):
                    cp_localization = complex_part_localization(reactant)
                    reactant_localization = reactant_localization | cp_localization

            assert len(reactant_localization) == 1

        assert [left_hand_side_reactant.validate() for left_hand_side_reactant in self.left_hand_side]
        assert [right_hand_side_reactant.validate() for right_hand_side_reactant in self.right_hand_side]

        localization_validation(self.left_hand_side)
        localization_validation(self.right_hand_side)


class Reactant:
    pass


class MoleculeReactant(Reactant):
    def __init__(self, molecule_specification: MoleculeSpecification):
        self.molecule_specification = molecule_specification

    def validate(self):
        self.molecule_specification.validate()


class ComplexReactant(Reactant):
    def __init__(self, complex_parts: List[MoleculeSpecification], complex_bindings: List['Binding']):
        self.complex_parts = complex_parts
        self.complex_bindings = complex_bindings
        self.validate()

    def validate(self):
        [part.validate() for part in self.complex_parts]
        [bind.validate() for bind in self.complex_bindings]
        # each binding pair should appear only once
        assert all(self.complex_bindings.count(bind) == 1 for bind in self.complex_bindings)

        cp_localization = complex_part_localization(self)
        assert len(cp_localization) == 1


class Binding:
    def __init__(self, left_partner: Tuple[int, AssociationSpecification], right_partner: Tuple[int, AssociationSpecification]):
        self.left_partner = left_partner
        self.right_partner = right_partner

    def validate(self):
        # @todo validate the binding indices as well?
        self.left_partner[1].validate()
        self.right_partner[1].validate()
        # the assoc domain should be occupied if we bound something to it
        assert self.left_partner[1].is_occupied
        assert self.right_partner[1].is_occupied


class Arrow(Enum):
    irreversible = '->'
    reversible   = '<->'


class Parameter:
    def __init__(self, name: str, value):
        self.name = name
        self.value = value


class InitialCondition:
    def __init__(self, molecule_specification: MoleculeSpecification, value):
        self.molecule_specification = molecule_specification
        self.value = value

    def validate(self):
        self.molecule_specification.validate()


def complex_part_localization(complex_molecule) -> set:
    complex_part_localization_set = set()
    for part in complex_molecule.complex_parts:
        complex_part_localization_set.add(part.localization_specifications[0].localization_definition.compartment)
    return complex_part_localization_set
