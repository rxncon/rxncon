from typing import List

import rxncon.simulation.rule_based.rule_based_model as rbm


class BNGLSystem:
    def __init__(self, rule_based_model: rbm.RuleBasedModel):
        self.rule_based_model = rule_based_model

    def to_string(self):
        pass

def string_from_molecule_reactant(molecule_reactant: rbm.MoleculeReactant) -> str:
    return string_from_molecule_specification(molecule_reactant.molecule_specification)

def string_from_complex_reactant(complex_reactant: rbm.ComplexReactant) -> str:

    complex = [string_from_molecule_specification(complex_part) for complex_part in complex_reactant.complex_parts]
    return ".".join(complex)

def string_from_molecule_definition(molecule_definition: rbm.MoleculeDefinition) -> str:
    if not molecule_definition.modification_definitions and not molecule_definition.association_definitions and not\
            molecule_definition.localization_definitions:
        return molecule_definition.name

    # todo how to handle loc domains it should be loc~Cyto~Nuc in the definition case but loc~Cytosol in the specific case
    domain_strs = [string_from_localization_definitions(molecule_definition.localization_definitions), \
                   ','.join([string_from_modification_definition(x) for x in molecule_definition.modification_definitions]), \
                   ','.join([string_from_association_definition(x) for x in molecule_definition.association_definitions])]
    return "{0}({1})".format(molecule_definition.name,",".join(domain_strs))


def string_from_modification_definition(modification_definition: rbm.ModificationDefinition) -> str:
    return '{0}~{1}'.format(modification_definition.domain_name, '~'.join(modification_definition.valid_modifiers))


def string_from_association_definition(association_definition: rbm.AssociationDefinition) -> str:
    return association_definition.domain_name


def string_from_localization_definitions(localization_defintions: List[rbm.LocalizationDefinition]) -> str:
    locs = [localization_defintion.compartment for localization_defintion in localization_defintions]
    return 'loc~{0}'.format("~".join(locs))


def string_from_localization_definition(localization_definition: rbm.LocalizationDefinition) -> str:
    return 'loc~{0}'.format(localization_definition.compartment)


def string_from_molecule_specification(molecule_specification: rbm.MoleculeSpecification) -> str:
    if not molecule_specification.modification_specifications and not molecule_specification.association_specifications and not\
            molecule_specification.localization_specifications:
        return molecule_specification.molecule_definition.name

    domain_strs = [','.join(string_from_localization_specification(x) for x in molecule_specification.localization_specifications), \
        ','.join(string_from_modification_specification(x) for x in molecule_specification.modification_specifications), \
        ','.join(string_from_association_specification(x) for x in molecule_specification.association_specifications)]

    return "{0}({1})".format(molecule_specification.molecule_definition.name,",".join(domain_strs))


def string_from_modification_specification(modification_specification: rbm.ModificationSpecification) -> str:
    return "{0}~{1}".format(modification_specification.modification_definition.domain_name, modification_specification.value)


def string_from_association_specification(association_specification: rbm.AssociationSpecification) -> str:
    return association_specification.association_definition.domain_name


def string_from_localization_specification(localization_specification: rbm.LocalizationSpecification) -> str:
    return 'loc~' + localization_specification.localization_definition.compartment





