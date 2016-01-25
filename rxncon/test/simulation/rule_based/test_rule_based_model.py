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


def test_molecule_specification_raise_double_domain():
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

    with pytest.raises(AssertionError):
        modification_specifications_err = [rbm.ModificationSpecification(modification_definitions[0], "P"),
                                           rbm.ModificationSpecification(modification_definitions[0], "P")]
        rbm.MoleculeSpecification(molecule_definition, modification_specifications_err,
                                                           association_specifications, localization_specifications)

    with pytest.raises(AssertionError):
        association_specifications_err = [rbm.AssociationSpecification(association_definitions[0], True),
                                          rbm.AssociationSpecification(association_definitions[0], True)]
        rbm.MoleculeSpecification(molecule_definition, modification_specifications,
                                  association_specifications_err, localization_specifications)

    with pytest.raises(AssertionError):
        # loc dom twice
        localization_specifications_err_dom = [rbm.LocalizationSpecification(localization_definitions[0], False),
                                               rbm.LocalizationSpecification(localization_definitions[0], False)]
        rbm.MoleculeSpecification(molecule_definition, modification_specifications,
                                  association_specifications, localization_specifications_err_dom)

    with pytest.raises(AssertionError):
        # in different compartments at the same time
        localization_specifications_err_loc = [rbm.LocalizationSpecification(localization_definitions[0], True),
                                               rbm.LocalizationSpecification(localization_definitions[1], True)]
        rbm.MoleculeSpecification(molecule_definition, modification_specifications,
                                  association_specifications, localization_specifications_err_loc)


def test_molecule_specification_raise_undefined_domain():
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

    #with pytest.raises(AssertionError):
    # modification_definition_raise = rbm.ModificationDefinition("ModDomain1", ["U", "ATP"])
    # modification_specifications_raise = [rbm.ModificationSpecification(modification_definition_raise, "ATP")]
    # rbm.MoleculeSpecification(molecule_definition, modification_specifications_raise,
    #                                                    association_specifications, localization_specifications)
    #
    #
    # with pytest.raises(AssertionError):
    #     association_def = rbm.AssociationDefinition("AssociationDomain3")
    #     association_specifications = [rbm.AssociationSpecification(association_def, True)]
    #     rbm.MoleculeSpecification(molecule_definition, modification_specifications,
    #                                                        association_specifications, localization_specifications)
    #
    # with pytest.raises(AssertionError):
    #     localization_def = rbm.LocalizationDefinition("Cell")
    #     localization_specifications = [rbm.LocalizationSpecification(localization_def, True)]
    #     rbm.MoleculeSpecification(molecule_definition, modification_specifications,
    #                                                        association_specifications, localization_specifications)

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


def test_binding():
    association_definitions_left = rbm.AssociationDefinition("AssociationDomain1")
    association_specifications_left = rbm.AssociationSpecification(association_definitions_left, True)

    association_definitions_right = rbm.AssociationDefinition("AssociationDomain2")
    association_specifications_right = rbm.AssociationSpecification(association_definitions_right, True)

    binding = rbm.Binding(left_partner=(0, association_specifications_left),
                          right_partner=(1, association_specifications_right))

    assert binding.left_partner[0] == 0
    assert binding.left_partner[1].is_occupied is True

    assert binding.right_partner[0] == 1
    assert binding.right_partner[1].is_occupied is True


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
    binding_list = [rbm.Binding(left_partner=(0, association_specifications1[0]),
                                right_partner=(1, association_specifications2[0]))]
    complex_reactant = rbm.ComplexReactant(mol_list, binding_list)

    assert len(complex_reactant.complex_parts) == 2

    assert complex_reactant.complex_parts[0].modification_specifications == modification_specifications
    assert complex_reactant.complex_parts[1].modification_specifications == []

    assert complex_reactant.complex_parts[0].association_specifications == association_specifications1
    assert complex_reactant.complex_parts[1].association_specifications == association_specifications2

    assert complex_reactant.complex_parts[0].localization_specifications == localization_specifications
    assert complex_reactant.complex_parts[1].localization_specifications == localization_specifications


def test_complex_reactant_parts_different_location():
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
    molecule_specification2 = rbm.MoleculeSpecification(molecule_definition, [], association_specifications2,
                                                        localization_specifications)

    binding_list = [rbm.Binding(left_partner=(0, association_specifications1[0]),
                                right_partner=(1, association_specifications2[0]))]
    with pytest.raises(AssertionError):
        localization_definitions_raise = [rbm.LocalizationDefinition("Nucleus")]
        localization_specifications_raise = [rbm.LocalizationSpecification(localization_definitions_raise[0], True)]

        molecule_definition_raise = rbm.MoleculeDefinition("A", modification_definitions, association_definitions,
                                                     localization_definitions_raise)

        molecule_specification1 = rbm.MoleculeSpecification(molecule_definition_raise, modification_specifications,
                                                            association_specifications1, localization_specifications_raise)

        mol_list = [molecule_specification1, molecule_specification2]
        rbm.ComplexReactant(mol_list, binding_list)


def test_rule_ppi():
    # A(modDom~P,AssocB) + B(AssocA) -> A(modDom~P,AssocB!1).B(AssocA!1)

    assoc_def = [rbm.AssociationDefinition("AssocDom1"), rbm.AssociationDefinition("AssocDom2")]

    localization_definitions = [rbm.LocalizationDefinition("Cytosole"), rbm.LocalizationDefinition("Nucleus")]
    localization_specifications = [rbm.LocalizationSpecification(localization_definitions[0], True)]

    A_modification_definition = [rbm.ModificationDefinition("ModDomain1", ["U", "P"])]
    A_modification_specification = [rbm.ModificationSpecification(A_modification_definition[0], "P")]
    A_association_specification = [rbm.AssociationSpecification(assoc_def[0], True)]

    A_molecule_definition = rbm.MoleculeDefinition("A", A_modification_definition, assoc_def, localization_definitions)
    A_molecule_specification = rbm.MoleculeSpecification(A_molecule_definition, A_modification_specification, A_association_specification, localization_specifications)

    B_assocociation_specification = [rbm.AssociationSpecification(assoc_def[1], True)]

    B_molecule_definition = rbm.MoleculeDefinition("B", [], assoc_def, localization_definitions)
    B_molecule_specification = rbm.MoleculeSpecification(B_molecule_definition, [], B_assocociation_specification, localization_specifications)
    A_reactant = rbm.MoleculeReactant(A_molecule_specification)
    B_reactant = rbm.MoleculeReactant(B_molecule_specification)

    left_hand_side = [A_reactant, B_reactant]

    binding = rbm.Binding((0, A_molecule_specification.association_specifications[0]),
                          (1, B_molecule_specification.association_specifications[0]))

    right_hand_side = [rbm.ComplexReactant([A_molecule_specification, B_molecule_specification],
                                           [binding])]

    rule = rbm.Rule(left_hand_side, right_hand_side, rbm.Arrow.reversible)
    assert len(rule.left_hand_side) == 2
    assert rule.left_hand_side[0] == A_reactant
    assert rule.left_hand_side[1] == B_reactant

    assert len(rule.right_hand_side) == 1
    assert len(rule.right_hand_side[0].complex_parts) == 2
    assert rule.right_hand_side[0].complex_parts == [A_molecule_specification, B_molecule_specification]
    assert len(rule.right_hand_side[0].complex_bindings) == 1
    assert rule.right_hand_side[0].complex_bindings == [binding]

    assert rule.arrow_type.name == "reversible"
    assert rule.arrow_type.value == "<->"


def test_test_rule_ppi_raise_reactant_different_location():
    assoc_def = [rbm.AssociationDefinition("AssocDom1"), rbm.AssociationDefinition("AssocDom2")]

    localization_definitions = [rbm.LocalizationDefinition("Cytosole"), rbm.LocalizationDefinition("Nucleus")]
    localization_specifications = [rbm.LocalizationSpecification(localization_definitions[0], True)]

    A_modification_definition = [rbm.ModificationDefinition("ModDomain1", ["U", "P"])]
    A_modification_specification = [rbm.ModificationSpecification(A_modification_definition[0], "P")]
    A_association_specification = [rbm.AssociationSpecification(assoc_def[0], True)]

    A_molecule_definition = rbm.MoleculeDefinition("A", A_modification_definition, assoc_def, localization_definitions)
    A_molecule_specification = rbm.MoleculeSpecification(A_molecule_definition, A_modification_specification, A_association_specification, localization_specifications)

    B_assocociation_specification = [rbm.AssociationSpecification(assoc_def[1], True)]

    A_reactant = rbm.MoleculeReactant(A_molecule_specification)

    with pytest.raises(AssertionError):
        # molecule should be at the same location to react with each other
        B_molecule_definition = rbm.MoleculeDefinition("B", [], assoc_def, localization_definitions)
        localization_specifications_raise = [rbm.LocalizationSpecification(localization_definitions[1], True)]
        B_molecule_specification = rbm.MoleculeSpecification(B_molecule_definition, [], B_assocociation_specification, localization_specifications_raise)

        B_reactant = rbm.MoleculeReactant(B_molecule_specification)

        binding = rbm.Binding((0, A_molecule_specification.association_specifications[0]),
                          (1, B_molecule_specification.association_specifications[0]))
        left_hand_side = [A_reactant, B_reactant]
        right_hand_side = [rbm.ComplexReactant([A_molecule_specification, B_molecule_specification],
                                           [binding])]
        rbm.Rule(left_hand_side, right_hand_side, rbm.Arrow.reversible)