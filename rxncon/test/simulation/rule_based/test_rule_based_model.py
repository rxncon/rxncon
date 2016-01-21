#mol_definition = MoleculeDefinition(…)
#rule1 = Rule(….)
#etc. etc.
#model = RuleBasedModel(mol_def, rules …)
import pytest
import rxncon.simulation.rule_based.rule_based_model as rbm


def test_modification_definition():
    modification_definition = rbm.ModificationDefinition("ModDomain", ["U","P"])
    assert modification_definition.domain_name == "ModDomain"
    assert modification_definition.valid_modifiers == ["U","P"]


def test_modification_specification():
    modification_definition = rbm.ModificationDefinition("ModDomain", ["U","P"])
    modification_specification = rbm.ModificationSpecification(modification_definition, "P")
    assert modification_specification.value == "P"

    with pytest.raises(AssertionError):
        rbm.ModificationSpecification(modification_definition, "G")


def test_association_definition():
    association_definition = rbm.AssociationDefinition("AssociationDomain")
    assert association_definition.domain_name == "AssociationDomain"


def test_association_specification():
    association_definition = rbm.AssociationDefinition("AssociationDomain")
    association_specification_occupied = rbm.AssociationSpecification(association_definition,True)
    association_specification_not_occupied = rbm.AssociationSpecification(association_definition, False)

    assert association_specification_occupied.association_definition.domain_name == "AssociationDomain"
    assert association_specification_occupied.is_occupied

    assert association_specification_not_occupied.association_definition.domain_name == "AssociationDomain"
    assert not association_specification_not_occupied.is_occupied


def test_localization_definition():
    localization_definition = rbm.LocalizationDefinition("Cell")
    assert localization_definition.compartment == "Cell"

def test_locatlization_specification():
    localization_definition = rbm.LocalizationDefinition("Cell")
    localization_specification = rbm.LocalizationSpecification(localization_definition, True)
    assert localization_specification.localization_definition.compartment == "Cell"

    with pytest.raises(AssertionError):
        rbm.LocalizationSpecification(localization_definition, "Cell")

def test_molecule_definition():
    modification_definitions = [rbm.ModificationDefinition("ModDomain1", ["U","P"]),
                                rbm.ModificationDefinition("ModDomain2", ["U","GTP"])]
    association_definitions = [rbm.AssociationDefinition("AssociationDomain1"),
                               rbm.AssociationDefinition("AssociationDomain2")]
    localization_definitions = [rbm.LocalizationDefinition("Cytosole"),
                               rbm.LocalizationDefinition("Nucleus")]

    molecule_definition = rbm.MoleculeDefinition("A",modification_definitions,association_definitions,localization_definitions)

    assert molecule_definition.name == "A"
    assert molecule_definition.modification_definitions == modification_definitions
    assert molecule_definition.association_definitions == association_definitions
    assert molecule_definition.localization_definitions == localization_definitions


def test_molecule_specification():
    modification_definitions = [rbm.ModificationDefinition("ModDomain1", ["U","P"]),
                                rbm.ModificationDefinition("ModDomain2", ["U","GTP"])]
    modification_specifications = [rbm.ModificationSpecification(modification_definitions[0], "P")]
    association_definitions = [rbm.AssociationDefinition("AssociationDomain1"),
                               rbm.AssociationDefinition("AssociationDomain2")]
    association_specifications = [rbm.AssociationSpecification(association_definitions[0],True)]
    localization_definitions = [rbm.LocalizationDefinition("Cytosole"),
                               rbm.LocalizationDefinition("Nucleus")]
    localization_specifications = [rbm.LocalizationSpecification(localization_definitions[0], True)]

    molecule_definition = rbm.MoleculeDefinition("A",modification_definitions,association_definitions,localization_definitions)

    molecule_specification = rbm.MoleculeSpecification(molecule_definition,modification_specifications,
                                                       association_specifications,localization_specifications)

    assert molecule_specification.molecule_definition == molecule_definition

    with pytest.raises(AssertionError):
        mod_def = rbm.ModificationDefinition("ModDomain1", ["U","ATP"])
        modification_specifications = [rbm.ModificationSpecification(mod_def, "ATP")]
        molecule_specification = rbm.MoleculeSpecification(molecule_definition,modification_specifications,
                                                       association_specifications,localization_specifications)

    with pytest.raises(AssertionError):
        association_def = rbm.AssociationDefinition("AssociationDomain3")
        association_specifications = [rbm.AssociationSpecification(association_def,True)]
        molecule_specification = rbm.MoleculeSpecification(molecule_definition,modification_specifications,
                                                       association_specifications,localization_specifications)
    with pytest.raises(AssertionError):
        localization_def = rbm.LocalizationDefinition("Cell")
        localization_specifications = [rbm.LocalizationSpecification(localization_def, True)]
        molecule_specification = rbm.MoleculeSpecification(molecule_definition,modification_specifications,
                                                       association_specifications,localization_specifications)
