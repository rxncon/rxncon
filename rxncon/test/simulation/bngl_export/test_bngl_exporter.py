import pytest

import rxncon.simulation.rule_based.rule_based_model as rbm
import rxncon.simulation.rule_based.bngl_export as be


def test_string_from_modification_definition():
    modification_definition = rbm.ModificationDefinition("ModDomain", ["U", "P"])
    mod_def_str = be.string_from_modification_definition(modification_definition)
    assert mod_def_str == "ModDomain~U~P"


def test_string_from_modification_specification():
    modification_definition = rbm.ModificationDefinition("ModDomain", ["U", "P"])
    modification_specification = rbm.ModificationSpecification(modification_definition, "P")

    mod_spec_str = be.string_from_modification_specification(modification_specification)
    assert mod_spec_str == "ModDomain~P"


def test_string_from_association_definition():
    association_definition = rbm.AssociationDefinition("AssociationDomain")
    assoc_def_str = be.string_from_association_definition(association_definition)
    assert assoc_def_str == "AssociationDomain"


def test_string_from_association_specification():
    association_definition = rbm.AssociationDefinition("AssociationDomain")
    association_specification_occupied = rbm.AssociationSpecification(association_definition, True)
    association_specification_not_occupied = rbm.AssociationSpecification(association_definition, False)

    assoc_spec_str_occ = be.string_from_association_specification(association_specification_occupied)
    assert assoc_spec_str_occ == "AssociationDomain"

    assoc_spec_str_not_occ = be.string_from_association_specification(association_specification_not_occupied)
    assert assoc_spec_str_not_occ == "AssociationDomain"


def test_string_from_localization_definition():
    localization_definition = rbm.LocalizationDefinition("Cell")
    loc_str = be.string_from_localization_definition(localization_definition)

    assert loc_str == "loc~Cell"


def test_string_from_localization_specification():
    localization_definition = rbm.LocalizationDefinition("Cell")
    localization_specification = rbm.LocalizationSpecification(localization_definition, True)
    loc_str = be.string_from_localization_specification(localization_specification)

    assert loc_str == "loc~Cell"


def test_string_from_molecule_definition():
    modification_definitions = [rbm.ModificationDefinition("ModDomain1", ["U", "P"]),
                                rbm.ModificationDefinition("ModDomain2", ["U", "GTP"])]
    association_definitions = [rbm.AssociationDefinition("AssociationDomain1"),
                               rbm.AssociationDefinition("AssociationDomain2")]
    localization_definitions = [rbm.LocalizationDefinition("Cytosole"),
                                rbm.LocalizationDefinition("Nucleus")]

    molecule_definition = rbm.MoleculeDefinition("A", modification_definitions, association_definitions,
                                                 localization_definitions)

    mol_def_str = be.string_from_molecule_definition(molecule_definition)

    assert mol_def_str == "A(loc~Cytosole~Nucleus,ModDomain1~U~P,ModDomain2~U~GTP,AssociationDomain1,AssociationDomain2)"


def test_string_from_molecule_specification():
    modification_definitions = [rbm.ModificationDefinition("ModDomain1", ["U", "P"]),
                                rbm.ModificationDefinition("ModDomain2", ["U", "GTP"])]
    modification_specifications = [rbm.ModificationSpecification(modification_definitions[0], "P")]
    association_definitions = [rbm.AssociationDefinition("AssociationDomain1"),
                               rbm.AssociationDefinition("AssociationDomain2")]
    association_specifications = [rbm.AssociationSpecification(association_definitions[0], True)]
    localization_definitions = [rbm.LocalizationDefinition("Cytosole"),
                                rbm.LocalizationDefinition("Nucleus")]
    localization_specifications = [rbm.LocalizationSpecification(localization_definitions[0], True)]

    molecule_definition = rbm.MoleculeDefinition("A", modification_definitions, association_definitions,
                                                 localization_definitions)

    molecule_specification = rbm.MoleculeSpecification(molecule_definition, modification_specifications,
                                                       association_specifications, localization_specifications)
    molecule_specification_str = be.string_from_molecule_specification(molecule_specification)

    assert molecule_specification_str == "A(loc~Cytosole,ModDomain1~P,AssociationDomain1)"
