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
        assert all(mod_spec.modification_definition in self.molecule_definition.modification_definitions
                   for mod_spec in self.modification_specifications)
        # a modification can only appear once
        assert all(self.modification_specifications.count(mod_spec) == 1 for mod_spec in self.modification_specifications)

        mod_dom = [mod_spec.modification_definition.domain_name for mod_spec in self.modification_specifications]
        # each mod domain should appear only once in molecule
        assert all(mod_dom.count(dom) == 1 for dom in mod_dom)

        assert all(ass_spec.association_definition in self.molecule_definition.association_definitions
                   for ass_spec in self.association_specifications)

        assoc_dom = [assoc_spec.association_definition.domain_name for assoc_spec in self.association_specifications]
        # each mod domain should appear only once in molecule
        assert all(assoc_dom.count(dom) == 1 for dom in assoc_dom)

        assert all(loc_spec.localization_definition in self.molecule_definition.localization_definitions
                   for loc_spec in self.localization_specifications)

        loc_dom = [loc_spec.localization_definition.compartment for loc_spec in self.localization_specifications]
        # each mod domain should appear only once in molecule
        assert all(loc_dom.count(dom) == 1 for dom in loc_dom)

        # a molecule can only be localised at one place at a time
        localisation = [loc_spec.localization_definition.compartment for loc_spec in self.localization_specifications if loc_spec.is_localized]
        assert len(localisation) == 1

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
        # the assoc domain should be occupied if we bound something to it
        # called in Binding()
        assert self.is_occupied


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

    def validate(self):
        self.left_hand_side.validate()
        self.right_hand_side.validate()


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
        localisation = [loc.localization_definition.compartment for loc in self.complex_parts[0].localization_specifications if loc.is_localized]
        # a molecule can be localised only at one place at a time
        assert len(localisation) == 1
        assert all(loc.localization_definition.compartment == localisation[0] for part in self.complex_parts for loc in part.localization_specifications if loc.is_localized)


class Binding:
    def __init__(self, left_partner: Tuple[int, AssociationSpecification], right_partner: Tuple[int, AssociationSpecification]):
        self.left_partner = left_partner
        self.right_partner = right_partner

    def validate(self):
        # @todo validate the binding indices as well?
        self.left_partner[1].validate()
        self.right_partner[1].validate()


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


