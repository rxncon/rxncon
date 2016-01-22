#mol_definition = MoleculeDefinition(…)
#rule1 = Rule(….)
#etc. etc.
#model = RuleBasedModel(mol_def, rules …)
import pytest
import rxncon.simulation.rule_based.rule_based_model as rbm


def test_modification_definition():
    modification_definition = rbm.ModificationDefinition("ModDomain", ["U", "P"])
    assert modification_definition.domain_name == "ModDomain"
    assert modification_definition.valid_modifiers == ["U", "P"]


def test_modification_specification():
    modification_definition = rbm.ModificationDefinition("ModDomain", ["U", "P"])
    modification_specification = rbm.ModificationSpecification(modification_definition, "P")
    assert modification_specification.value == "P"

    with pytest.raises(AssertionError):
        rbm.ModificationSpecification(modification_definition, "G")


def test_association_definition():
    association_definition = rbm.AssociationDefinition("AssociationDomain")
    assert association_definition.domain_name == "AssociationDomain"


def test_association_specification():
    association_definition = rbm.AssociationDefinition("AssociationDomain")
    association_specification_occupied = rbm.AssociationSpecification(association_definition, True)
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
    modification_definitions = [rbm.ModificationDefinition("ModDomain1", ["U", "P"]),
                                rbm.ModificationDefinition("ModDomain2", ["U", "GTP"])]
    association_definitions = [rbm.AssociationDefinition("AssociationDomain1"),
                               rbm.AssociationDefinition("AssociationDomain2")]
    localization_definitions = [rbm.LocalizationDefinition("Cytosole"),
                                rbm.LocalizationDefinition("Nucleus")]

    molecule_definition = rbm.MoleculeDefinition("A", modification_definitions, association_definitions,
                                                 localization_definitions)

    assert molecule_definition.name == "A"
    assert molecule_definition.modification_definitions == modification_definitions
    assert molecule_definition.association_definitions == association_definitions
    assert molecule_definition.localization_definitions == localization_definitions


def test_molecule_specification():
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

    assert molecule_specification.molecule_definition == molecule_definition

    with pytest.raises(AssertionError):
        modification_specifications_err = [rbm.ModificationSpecification(modification_definitions[0], "P"),
                                       rbm.ModificationSpecification(modification_definitions[0], "P")]
        molecule_specification = rbm.MoleculeSpecification(molecule_definition, modification_specifications_err,
                                                           association_specifications, localization_specifications)

    with pytest.raises(AssertionError):
        association_specifications_err = [rbm.AssociationSpecification(association_definitions[0], True),
                                          rbm.AssociationSpecification(association_definitions[0], True)]
        molecule_specification = rbm.MoleculeSpecification(molecule_definition, modification_specifications,
                                                           association_specifications_err, localization_specifications)

    with pytest.raises(AssertionError):
        # loc dom twice
        localization_specifications_err_dom = [rbm.LocalizationSpecification(localization_definitions[0], False),
                                               rbm.LocalizationSpecification(localization_definitions[0], False)]
        molecule_specification = rbm.MoleculeSpecification(molecule_definition, modification_specifications,
                                                           association_specifications, localization_specifications_err_dom)

    with pytest.raises(AssertionError):
        # in different compartments at the same time
        localization_specifications_err_loc = [rbm.LocalizationSpecification(localization_definitions[0], True),
                                       rbm.LocalizationSpecification(localization_definitions[1], True)]
        molecule_specification = rbm.MoleculeSpecification(molecule_definition, modification_specifications,
                                                           association_specifications, localization_specifications_err_loc)



def test_molecule_reactant():
    modification_definitions = [rbm.ModificationDefinition("ModDomain1", ["U", "P"])]
    modification_specifications = [rbm.ModificationSpecification(modification_definitions[0], "P")]
    association_definitions = [rbm.AssociationDefinition("AssociationDomain1")]
    association_specifications = [rbm.AssociationSpecification(association_definitions[0], True)]
    localization_definitions = [rbm.LocalizationDefinition("Cytosole")]
    localization_specifications = [rbm.LocalizationSpecification(localization_definitions[0], True)]

    molecule_definition = rbm.MoleculeDefinition("A", modification_definitions, association_definitions,
                                                 localization_definitions)

    molecule_specification = rbm.MoleculeSpecification(molecule_definition, modification_specifications,
                                                       association_specifications, localization_specifications)

    molecule_reactant = rbm.MoleculeReactant(molecule_specification)

    assert molecule_reactant.molecule_specification == molecule_specification

    with pytest.raises(AssertionError):
        mod_def = rbm.ModificationDefinition("ModDomain1", ["U", "ATP"])
        modification_specifications = [rbm.ModificationSpecification(mod_def, "ATP")]
        molecule_specification = rbm.MoleculeSpecification(molecule_definition, modification_specifications,
                                                       association_specifications, localization_specifications)
        molecule_reactant = rbm.MoleculeReactant(molecule_specification)
        molecule_reactant.validate()

    with pytest.raises(AssertionError):
        association_def = rbm.AssociationDefinition("AssociationDomain3")
        association_specifications = [rbm.AssociationSpecification(association_def, True)]
        molecule_specification = rbm.MoleculeSpecification(molecule_definition, modification_specifications,
                                                       association_specifications, localization_specifications)
        molecule_reactant = rbm.MoleculeReactant(molecule_specification)
        molecule_reactant.validate()

    with pytest.raises(AssertionError):
        localization_def = rbm.LocalizationDefinition("Cell")
        localization_specifications = [rbm.LocalizationSpecification(localization_def, True)]
        molecule_specification = rbm.MoleculeSpecification(molecule_definition, modification_specifications,
                                                       association_specifications, localization_specifications)

        molecule_reactant = rbm.MoleculeReactant(molecule_specification)
        molecule_reactant.validate()


def test_binding():
    association_definitions_left = rbm.AssociationDefinition("AssociationDomain1")
    association_specifications_left = rbm.AssociationSpecification(association_definitions_left, True)

    association_definitions_right = rbm.AssociationDefinition("AssociationDomain2")
    association_specifications_right = rbm.AssociationSpecification(association_definitions_right, True)

    binding = rbm.Binding(left_partner=(0, association_specifications_left), right_partner=(1, association_specifications_right))

    assert binding.left_partner[0] == 0
    assert binding.left_partner[1].is_occupied == True

    assert binding.right_partner[0] == 1
    assert binding.right_partner[1].is_occupied == True

def test_complex_reactant():
    modification_definitions = [rbm.ModificationDefinition("ModDomain1", ["U", "P"])]
    modification_specifications = [rbm.ModificationSpecification(modification_definitions[0], "P")]
    association_definitions = [rbm.AssociationDefinition("AssociationDomain1"),
                               rbm.AssociationDefinition("AssociationDomain2")]
    association_specifications1 = [rbm.AssociationSpecification(association_definitions[0], True)]
    association_specifications2 = [rbm.AssociationSpecification(association_definitions[1], True)]

    localization_definitions = [rbm.LocalizationDefinition("Cytosole")]
    localization_specifications = [rbm.LocalizationSpecification(localization_definitions[0], True)]

    molecule_definition = rbm.MoleculeDefinition("A", modification_definitions, association_definitions,
                                                 localization_definitions)

    molecule_specification1 = rbm.MoleculeSpecification(molecule_definition, modification_specifications,
                                                       association_specifications1, localization_specifications)

    molecule_specification2 = rbm.MoleculeSpecification(molecule_definition, [], association_specifications2,
                                                        localization_specifications)

    mol_list = [molecule_specification1, molecule_specification2]
    binding_list = [rbm.Binding(left_partner=(0, association_specifications1[0]), right_partner=(1, association_specifications2[0]))]
    complex_reactant = rbm.ComplexReactant(mol_list, binding_list)

    assert len(complex_reactant.complex_parts) == 2

    assert complex_reactant.complex_parts[0].modification_specifications == modification_specifications
    assert complex_reactant.complex_parts[1].modification_specifications == []

    assert complex_reactant.complex_parts[0].association_specifications == association_specifications1
    assert complex_reactant.complex_parts[1].association_specifications == association_specifications2

    assert complex_reactant.complex_parts[0].localization_specifications == localization_specifications
    assert complex_reactant.complex_parts[1].localization_specifications == localization_specifications


def test_rule_ppi():
    # A(modDom~P,AssocB) + B(AssocA) -> A(modDom~P,AssocB!1).B(AssocA!1)

    assoc_def = [rbm.AssociationDefinition("AssocDom1"),
                   rbm.AssociationDefinition("AssocDom2")]

    localization_definitions = [rbm.LocalizationDefinition("Cytosole")]
    localization_specifications = [rbm.LocalizationSpecification(localization_definitions[0], True)]

    A_mod_def = [rbm.ModificationDefinition("ModDomain1", ["U", "P"])]
    A_mod_spec = [rbm.ModificationSpecification(A_mod_def[0], "P")]
    A_assoc_spec = [rbm.AssociationSpecification(assoc_def[0], True)]

    A_mol_def = rbm.MoleculeDefinition("A", A_mod_def, assoc_def, localization_definitions)
    A_mol_spec = rbm.MoleculeSpecification(A_mol_def, A_mod_spec, A_assoc_spec, localization_specifications)


    B_assoc_spec = [rbm.AssociationSpecification(assoc_def[1], True)]

    B_mol_def = rbm.MoleculeDefinition("B", [], assoc_def, localization_definitions)
    B_mol_spec = rbm.MoleculeSpecification(B_mol_def, [], B_assoc_spec, localization_specifications)
    A_reactant = rbm.MoleculeReactant(A_mol_spec)
    B_reactant = rbm.MoleculeReactant(B_mol_spec)
    left_hand_side = [A_reactant, B_reactant]

    binding = rbm.Binding((0, A_mol_spec.association_specifications[0]),
                          (1, B_mol_spec.association_specifications[0]))

    right_hand_side = [rbm.ComplexReactant([A_mol_spec, B_mol_spec],
                                           [binding])]

    rule = rbm.Rule(left_hand_side, right_hand_side, rbm.Arrow.reversible)
    assert len(rule.left_hand_side) == 2
    assert rule.left_hand_side[0] == A_reactant
    assert rule.left_hand_side[1] == B_reactant

    assert len(rule.right_hand_side) == 1
    assert len(rule.right_hand_side[0].complex_parts) == 2
    assert rule.right_hand_side[0].complex_parts == [A_mol_spec, B_mol_spec]
    assert len(rule.right_hand_side[0].complex_bindings) == 1
    assert rule.right_hand_side[0].complex_bindings == [binding]

    assert rule.arrow_type.name == "reversible"
    assert rule.arrow_type.value == "<->"
