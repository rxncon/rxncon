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
    localization_definition = rbm.LocalizationDefinition(["Cell","Nucleus"])
    loc_str = be.string_from_localization_definition(localization_definition)

    assert loc_str == "loc~Cell~Nucleus"


def test_string_from_localization_specification():
    localization_definition = rbm.LocalizationDefinition(["Cell"])
    localization_specification = rbm.LocalizationSpecification(localization_definition, "Cell")
    loc_str = be.string_from_localization_specification(localization_specification)

    assert loc_str == "loc~Cell"


def test_string_from_molecule_definition(molecule_A):

    mol_def_str = be.string_from_molecule_definition(molecule_A.molecule_definition)

    assert mol_def_str == "A(loc~Cytosole~Nucleus,ModDomain1~U~P,ModDomain2~U~GTP,AssociationDomain1,AssociationDomain2)"


def test_string_from_molecule_specification(molecule_A):

    molecule_specification_str = be.string_from_molecule_specification(molecule_A)

    assert molecule_specification_str == "A(loc~Cytosole,ModDomain1~P,AssociationDomain1)"


def test_string_from_molecule_reactant(molecule_A):

    molecule_reactant = rbm.MoleculeReactant(molecule_A)
    molecule_reactant_str = be.string_from_molecule_reactant(molecule_reactant)

    assert molecule_reactant_str == "A(loc~Cytosole,ModDomain1~P,AssociationDomain1)"


def test_string_from_complex_reactant(molecule_A):
    molecule_A.association_specifications[0].is_occupied = True
    molecule_specification_A_bound = molecule_A

    association_definitions_B_bound = [rbm.AssociationDefinition("AssociationDomain")]
    association_specifications_B_bound = [rbm.AssociationSpecification(association_definitions_B_bound[0], True)]
    localization_definitions_B_bound = [rbm.LocalizationDefinition(["Cytosole", "Nucleus"])]
    localization_specifications_B_bound = [rbm.LocalizationSpecification(localization_definitions_B_bound[0], "Cytosole")]

    molecule_definition_B_bound = rbm.MoleculeDefinition("B", [], association_definitions_B_bound,
                                                 localization_definitions_B_bound)

    molecule_specification_B_bound = rbm.MoleculeSpecification(molecule_definition_B_bound, [],
                                                       association_specifications_B_bound, localization_specifications_B_bound)

    molecules_bound = [molecule_specification_A_bound, molecule_specification_B_bound]

    binding_list = [rbm.Binding(left_partner=(0, molecule_specification_A_bound.association_specifications[0]),
                                right_partner=(1, molecule_specification_B_bound.association_specifications[0]))]
    complex_reactant = rbm.ComplexReactant(molecules_bound, binding_list)

    complex_reactant_str = be.string_from_complex_reactant(complex_reactant)

    assert complex_reactant_str == 'A(loc~Cytosole,ModDomain1~P,AssociationDomain1!0).B(loc~Cytosole,AssociationDomain!0)'

def test_simple_rule_based_model(simple_rule_based_model):
    bngl = be.BNGLSystem(simple_rule_based_model)

    expected_string = """begin model
begin parameters
A_total 100
B_total 100
NA 6.02214e23
V 1e-12
kf 1e6/(NA*V)
kr 0.1
end parameters

begin molecule types
A(b)
B(a)
end molecule types

begin seed species
A(b) A_total
B(a) B_total
end seed species

begin reaction rules
A(b) + B(a) <-> A(b!0).B(a!0) kf,kr
end reaction rules

end model

generate_network(max_iter=>1, max_agg=>4)
simulate({method=>"ode",t_end=>10,n_steps=>100})"""

    bngl.to_string()
    assert bngl.to_string() == expected_string

@pytest.fixture
def molecule_A():
    modification_definitions = [rbm.ModificationDefinition("ModDomain1", ["U", "P"]),
                                rbm.ModificationDefinition("ModDomain2", ["U", "GTP"])]
    modification_specifications = [rbm.ModificationSpecification(modification_definitions[0], "P")]
    association_definitions = [rbm.AssociationDefinition("AssociationDomain1"),
                               rbm.AssociationDefinition("AssociationDomain2")]
    association_specifications = [rbm.AssociationSpecification(association_definitions[0], False)]
    localization_definitions = [rbm.LocalizationDefinition(["Cytosole", "Nucleus"])]
    localization_specifications = [rbm.LocalizationSpecification(localization_definitions[0], "Cytosole")]

    molecule_definition = rbm.MoleculeDefinition("A", modification_definitions, association_definitions,
                                                 localization_definitions)

    molecule_specification = rbm.MoleculeSpecification(molecule_definition, modification_specifications,
                                                       association_specifications, localization_specifications)

    return molecule_specification

@pytest.fixture
def simple_rule_based_model():
    assoc_def_A = rbm.AssociationDefinition('b')
    assoc_def_B = rbm.AssociationDefinition('a')

    mol_def_A = rbm.MoleculeDefinition('A', [], [assoc_def_A], [])
    mol_def_B = rbm.MoleculeDefinition('B', [], [assoc_def_B], [])

    mol_spec_A_unbound = rbm.MoleculeSpecification(mol_def_A, [], [rbm.AssociationSpecification(assoc_def_A, False)], [])
    mol_spec_B_unbound = rbm.MoleculeSpecification(mol_def_B, [], [rbm.AssociationSpecification(assoc_def_B, False)], [])

    mol_spec_A_bound = rbm.MoleculeSpecification(mol_def_A, [], [rbm.AssociationSpecification(assoc_def_A, True)], [])
    mol_spec_B_bound = rbm.MoleculeSpecification(mol_def_B, [], [rbm.AssociationSpecification(assoc_def_B, True)], [])

    reactant_A_unbound = rbm.MoleculeReactant(mol_spec_A_unbound)
    reactant_B_unbound = rbm.MoleculeReactant(mol_spec_B_unbound)

    binding = rbm.Binding((0, mol_spec_A_bound.association_specifications[0]),
                          (1, mol_spec_B_bound.association_specifications[0]))

    reactant_AB = rbm.ComplexReactant([mol_spec_A_bound, mol_spec_B_bound], [binding])

    rule_kinetic_paramerters = [rbm.Parameter('kf', '1e6/(NA*V)'),
                                rbm.Parameter('kr', '0.1')]
    rule = rbm.Rule([reactant_A_unbound, reactant_B_unbound], [reactant_AB], rbm.Arrow.reversible, rule_kinetic_paramerters)

    parameters = [
        rbm.Parameter('A_total', '100'),
        rbm.Parameter('B_total', '100'),
        rbm.Parameter('NA', '6.02214e23'),
        rbm.Parameter('V', '1e-12'),
    ]

    parameters.extend(rule_kinetic_paramerters)

    initial_conditions = [
        rbm.InitialCondition(mol_spec_A_unbound, 'A_total'),
        rbm.InitialCondition(mol_spec_B_unbound, 'B_total')
    ]

    return rbm.RuleBasedModel([mol_def_A, mol_def_B], [rule], parameters, initial_conditions)








