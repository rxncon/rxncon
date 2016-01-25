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
    localization_definition = rbm.LocalizationDefinition(["Cell", "Cytoplasm", "Nucleus"])
    assert localization_definition.compartments == ["Cell", "Cytoplasm", "Nucleus"]


def test_localization_specification():
    localization_definition = rbm.LocalizationDefinition(["Cell", "Cytoplasm", "Nucleus"])
    localization_specification = rbm.LocalizationSpecification(localization_definition, "Cell")
    assert localization_specification.current_compartment == "Cell"

    #with pytest.raises(AssertionError):
    #    rbm.LocalizationSpecification(localization_definition, "Cell")


def test_molecule_definition():
    modification_definitions = [rbm.ModificationDefinition("ModDomain1", ["U", "P"]),
                                rbm.ModificationDefinition("ModDomain2", ["U", "GTP"])]
    association_definitions = [rbm.AssociationDefinition("AssociationDomain1"),
                               rbm.AssociationDefinition("AssociationDomain2")]
    localization_definitions = [rbm.LocalizationDefinition(["Cell", "Cytoplasm", "Nucleus"])]

    molecule_definition = rbm.MoleculeDefinition("A", modification_definitions, association_definitions,
                                                 localization_definitions)

    assert molecule_definition.name == "A"
    assert molecule_definition.modification_definitions == modification_definitions
    assert molecule_definition.association_definitions == association_definitions
    assert molecule_definition.localization_definitions == localization_definitions



def test_molecule_specification(modification_definitions, association_definitions, localization_definitions):
    modification_definitions = modification_definitions
    modification_specifications = [rbm.ModificationSpecification(modification_definitions[0], "P")]

    association_definitions = association_definitions
    association_specifications = [rbm.AssociationSpecification(association_definitions[0], True)]

    localization_definitions = localization_definitions
    localization_specifications = [rbm.LocalizationSpecification(localization_definitions[0], "Cytoplasm")]

    molecule_definition = rbm.MoleculeDefinition("A", modification_definitions, association_definitions,
                                                 localization_definitions)

    molecule_specification = rbm.MoleculeSpecification(molecule_definition, modification_specifications,
                                                       association_specifications, localization_specifications)

    assert molecule_specification.molecule_definition == molecule_definition



def test_molecule_reactant(modification_definitions, association_definitions, localization_definitions):
    modification_definitions = modification_definitions
    modification_specifications = [rbm.ModificationSpecification(modification_definitions[0], "P")]

    association_definitions = association_definitions
    association_specifications = [rbm.AssociationSpecification(association_definitions[0], True)]

    localization_definitions = localization_definitions
    localization_specifications = [rbm.LocalizationSpecification(localization_definitions[0], "Cytoplasm")]

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


def test_complex_reactant(molecules_bound):

    binding_list = [rbm.Binding(left_partner=(0, molecules_bound[0].association_specifications[0]),
                                right_partner=(1, molecules_bound[1].association_specifications[0]))]
    complex_reactant = rbm.ComplexReactant(molecules_bound, binding_list)

    assert len(complex_reactant.complex_parts) == 2

    assert complex_reactant.complex_parts[0].modification_specifications == molecules_bound[0].modification_specifications
    assert complex_reactant.complex_parts[1].modification_specifications == []

    assert complex_reactant.complex_parts[0].association_specifications == molecules_bound[0].association_specifications
    assert complex_reactant.complex_parts[1].association_specifications == molecules_bound[1].association_specifications

    assert complex_reactant.complex_parts[0].localization_specifications == molecules_bound[0].localization_specifications
    assert complex_reactant.complex_parts[1].localization_specifications == molecules_bound[1].localization_specifications


def test_complex_reactant_parts_different_location(molecules_bound):

    binding_list = [rbm.Binding(left_partner=(0, molecules_bound[0].association_specifications[0]),
                                right_partner=(1, molecules_bound[1].association_specifications[0]))]

    localization_specifications_raise = [rbm.LocalizationSpecification(localization_definitions, "Nucleus")]

    molecules_bound[1] = rbm.MoleculeSpecification(molecules_bound[1].molecule_definition, molecules_bound[1].modification_specifications,
                                                            molecules_bound[1].association_specifications, localization_specifications_raise)
    molecules_bound = [molecules_bound[0], molecules_bound[1]]

    with pytest.raises(AssertionError):
        rbm.ComplexReactant(molecules_bound, binding_list)


def test_rule_ppi(molecules_bound, modification_definitions, association_definitions, localization_definitions):
    # A(modDom~P,AssocB) + B(AssocA) -> A(modDom~P,AssocB!1).B(AssocA!1)

    molecule_specification_A_bound = molecules_bound[0]
    molecule_specification_B_bound = molecules_bound[1]


    association_specification_A_unbound = [rbm.AssociationSpecification(molecule_specification_A_bound.association_specifications[0], False)]
    assocociation_specification_B_unbound = [rbm.AssociationSpecification(molecule_specification_B_bound.association_specifications[0], False)]

    molecule_definition_A = rbm.MoleculeDefinition("A", modification_definitions, association_definitions , localization_definitions)
    molecule_specification_A_unbound = rbm.MoleculeSpecification(molecule_definition_A, molecule_specification_A_bound.modification_specifications, association_specification_A_unbound, molecule_specification_A_bound.localization_specifications)

    molecule_definition_B = rbm.MoleculeDefinition("B", [], association_definitions, localization_definitions)
    molecule_specification_B_unbound = rbm.MoleculeSpecification(molecule_definition_B, [], assocociation_specification_B_unbound, molecule_specification_B_bound.localization_specifications)

    reactant_A_undbound = rbm.MoleculeReactant(molecule_specification_A_unbound)
    reactant_B_unbound = rbm.MoleculeReactant(molecule_specification_B_unbound)

    left_hand_side = [reactant_A_undbound, reactant_B_unbound]

    binding = rbm.Binding((0, molecule_specification_A_bound.association_specifications[0]),
                          (1, molecule_specification_B_bound.association_specifications[0]))

    right_hand_side = [rbm.ComplexReactant([molecule_specification_A_bound, molecule_specification_B_bound],
                                           [binding])]

    rule = rbm.Rule(left_hand_side, right_hand_side, rbm.Arrow.reversible)
    assert len(rule.left_hand_side) == 2
    assert rule.left_hand_side[0] == reactant_A_undbound
    assert rule.left_hand_side[1] == reactant_B_unbound

    assert len(rule.right_hand_side) == 1
    assert len(rule.right_hand_side[0].complex_parts) == 2
    assert rule.right_hand_side[0].complex_parts == [molecule_specification_A_bound, molecule_specification_B_bound]
    assert len(rule.right_hand_side[0].complex_bindings) == 1
    assert rule.right_hand_side[0].complex_bindings == [binding]

    assert rule.arrow_type.name == "reversible"
    assert rule.arrow_type.value == "<->"


@pytest.fixture
def modification_definitions():
    return [rbm.ModificationDefinition("ModDomain1", ["U", "P"]),
            rbm.ModificationDefinition("ModDomain2", ["U", "GTP"])]


@pytest.fixture
def association_definitions():
    return [rbm.AssociationDefinition("AssociationDomain1"),
            rbm.AssociationDefinition("AssociationDomain2")]


@pytest.fixture
def localization_definitions():
    return [rbm.LocalizationDefinition(["Cell", "Cytoplasm", "Nucleus"])]


@pytest.fixture
def molecules_bound(modification_definitions, association_definitions, localization_definitions):
    modification_definitions = modification_definitions
    modification_specifications_A = [rbm.ModificationSpecification(modification_definitions[0], "P")]

    association_definitions = association_definitions

    association_specifications_A = [rbm.AssociationSpecification(association_definitions[0], True)]
    association_specifications_B = [rbm.AssociationSpecification(association_definitions[1], True)]

    localization_definitions = localization_definitions
    localization_specification_A = [rbm.LocalizationSpecification(localization_definitions[0], "Cytoplasm")]
    localization_specification_B = [rbm.LocalizationSpecification(localization_definitions[0], "Cytoplasm")]

    molecule_definition = rbm.MoleculeDefinition("A", modification_definitions, association_definitions,
                                                 localization_definitions)

    molecule_specification_A = rbm.MoleculeSpecification(molecule_definition, modification_specifications_A,
                                                         association_specifications_A, localization_specification_A)

    molecule_specification_B = rbm.MoleculeSpecification(molecule_definition, [], association_specifications_B,
                                                         localization_specification_B)

    return [molecule_specification_A, molecule_specification_B]


@pytest.fixture
def molecules_unbound(molecules_bound):
    molecule_specification_A_bound = molecules_bound[0]
    molecule_specification_B_bound = molecules_bound[1]


    association_specification_A_unbound = [rbm.AssociationSpecification(molecule_specification_A_bound.association_specifications[0], False)]
    assocociation_specification_B_unbound = [rbm.AssociationSpecification(molecule_specification_B_bound.association_specifications[0], False)]

    molecule_definition_A = rbm.MoleculeDefinition("A", modification_definitions, association_definitions , localization_definitions)
    molecule_specification_A_unbound = rbm.MoleculeSpecification(molecule_definition_A, molecule_specification_A_bound.modification_specifications, association_specification_A_unbound, molecule_specification_A_bound.localization_specifications)

    molecule_definition_B = rbm.MoleculeDefinition("B", [], association_definitions, localization_definitions)
    molecule_specification_B_unbound = rbm.MoleculeSpecification(molecule_definition_B, [], assocociation_specification_B_unbound, molecule_specification_B_bound.localization_specifications)

    return [molecule_specification_A_unbound, molecule_specification_B_unbound]


def test_ppi_rule_raise_reactant_different_location(molecules_unbound, molecules_bound, localization_definitions):
    molecule_specification_A_bound = molecules_bound[0]
    molecule_specification_B_bound = molecules_bound[1]

    molecule_specification_A_unbound = molecules_unbound[0]
    molecule_specification_B_unbound = molecules_unbound[1]

    molecule_specification_B_unbound.localization_specifications = [rbm.LocalizationSpecification(localization_definitions[0], "Nucleus")]
    reactant_A_unbound = rbm.MoleculeReactant(molecule_specification_A_unbound)
    reactant_B_unbound = rbm.MoleculeReactant(molecule_specification_B_unbound)

    binding = rbm.Binding((0, molecule_specification_A_bound.association_specifications[0]),
                          (1, molecule_specification_B_bound.association_specifications[0]))
    left_hand_side = [reactant_A_unbound, reactant_B_unbound]

    right_hand_side = [rbm.ComplexReactant([molecule_specification_A_bound, molecule_specification_B_bound],
                                           [binding])]

    with pytest.raises(AssertionError):
        # molecule should be at the same location to react with each other
        rbm.Rule(left_hand_side, right_hand_side, rbm.Arrow.reversible)


# def test_molecule_specification_raise_double_domain():
#     modification_definitions = [rbm.ModificationDefinition("ModDomain1", ["U", "P"]),
#                                 rbm.ModificationDefinition("ModDomain2", ["U", "GTP"])]
#     modification_specifications = [rbm.ModificationSpecification(modification_definitions[0], "P")]
#     association_definitions = [rbm.AssociationDefinition("AssociationDomain1"),
#                                rbm.AssociationDefinition("AssociationDomain2")]
#     association_specifications = [rbm.AssociationSpecification(association_definitions[0], True)]
#     localization_definitions = [rbm.LocalizationDefinition("Cytosole"),
#                                 rbm.LocalizationDefinition("Nucleus")]
#     localization_specifications = [rbm.LocalizationSpecification(localization_definitions[0], True)]
#
#     molecule_definition = rbm.MoleculeDefinition("A", modification_definitions, association_definitions,
#                                                  localization_definitions)
#
#     with pytest.raises(AssertionError):
#         modification_specifications_err = [rbm.ModificationSpecification(modification_definitions[0], "P"),
#                                            rbm.ModificationSpecification(modification_definitions[0], "P")]
#         rbm.MoleculeSpecification(molecule_definition, modification_specifications_err,
#                                                            association_specifications, localization_specifications)
#
#     with pytest.raises(AssertionError):
#         association_specifications_err = [rbm.AssociationSpecification(association_definitions[0], True),
#                                           rbm.AssociationSpecification(association_definitions[0], True)]
#         rbm.MoleculeSpecification(molecule_definition, modification_specifications,
#                                   association_specifications_err, localization_specifications)
#
#     with pytest.raises(AssertionError):
#         # loc dom twice
#         localization_specifications_err_dom = [rbm.LocalizationSpecification(localization_definitions[0], False),
#                                                rbm.LocalizationSpecification(localization_definitions[0], False)]
#         rbm.MoleculeSpecification(molecule_definition, modification_specifications,
#                                   association_specifications, localization_specifications_err_dom)
#
#     with pytest.raises(AssertionError):
#         # in different compartments at the same time
#         localization_specifications_err_loc = [rbm.LocalizationSpecification(localization_definitions[0], True),
#                                                rbm.LocalizationSpecification(localization_definitions[1], True)]
#         rbm.MoleculeSpecification(molecule_definition, modification_specifications,
#                                   association_specifications, localization_specifications_err_loc)


def test_molecule_specification_raise_undefined_domain(modification_definitions, association_definitions, localization_definitions):
    modification_definitions = modification_definitions
    modification_specifications = [rbm.ModificationSpecification(modification_definitions[0], "P")]

    association_definitions = association_definitions
    association_specifications = [rbm.AssociationSpecification(association_definitions[0], True)]

    localization_definitions = localization_definitions
    localization_specifications = [rbm.LocalizationSpecification(localization_definitions[0], "Cytoplasm")]

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

