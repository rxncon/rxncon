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
    localization_specification = rbm.LocalizationSpecification(localization_definition, "Cell")
    assert localization_specification.localization_definition.compartment == "Cell"



    #mol_defintion = rbm.MoleculeDefinition("A",)