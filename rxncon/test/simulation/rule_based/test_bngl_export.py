import pytest

import rxncon.simulation.rule_based.rule_based_model as rbm
import rxncon.simulation.rule_based.bngl_export as bex


def test_string_from_modification_definition():
    modification_definition = rbm.ModificationDefinition('ModDomain', ['U', 'P'])
    mod_def_str = bex.string_from_modification_definition(modification_definition)
    assert mod_def_str == 'ModDomain~U~P'


def test_string_from_modification_specification():
    modification_definition = rbm.ModificationDefinition('ModDomain', ['U', 'P'])
    modification_specification = rbm.ModificationSpecification(modification_definition, 'P')

    mod_spec_str = bex.string_from_modification_specification(modification_specification)
    assert mod_spec_str == 'ModDomain~P'


def test_string_from_association_definition():
    association_definition = rbm.AssociationDefinition('AssociationDomain')
    assoc_def_str = bex.string_from_association_definition(association_definition)
    assert assoc_def_str == 'AssociationDomain'


def test_string_from_association_specification():
    ass_def = rbm.AssociationDefinition('AssociationDomain')

    assert bex.string_from_association_specification(rbm.AssociationSpecification(ass_def, rbm.OccupationStatus.occupied_unknown_partner)) == \
        'AssociationDomain!+'

    assert bex.string_from_association_specification(rbm.AssociationSpecification(ass_def, rbm.OccupationStatus.unoccupied)) == \
        'AssociationDomain'

    assert bex.string_from_association_specification(rbm.AssociationSpecification(ass_def, rbm.OccupationStatus.not_specified)) == \
        'AssociationDomain!?'

    # This should not appear outside of a ComplexReactant.
    with pytest.raises(NotImplementedError):
        bex.string_from_association_specification(rbm.AssociationSpecification(ass_def, rbm.OccupationStatus.occupied_known_partner))


def test_string_from_localization_definition():
    localization_definition = rbm.LocalizationDefinition(['Cell', 'Nucleus'])
    loc_str = bex.string_from_localization_definition(localization_definition)

    assert loc_str == 'loc~Cell~Nucleus'


def test_string_from_localization_specification():
    localization_definition = rbm.LocalizationDefinition(['Cell'])
    localization_specification = rbm.LocalizationSpecification(localization_definition, 'Cell')
    loc_str = bex.string_from_localization_specification(localization_specification)

    assert loc_str == 'loc~Cell'


def test_string_from_molecule_definition(simple_molecule_spec):
    assert bex.string_from_molecule_definition(simple_molecule_spec.molecule_def) == \
        'A(loc~Cytosole~Nucleus,ModDomain1~U~P,ModDomain2~U~GTP,AssociationDomain1,AssociationDomain2)'


def test_string_from_molecule_specification(simple_molecule_spec):
    assert bex.string_from_molecule_specification(simple_molecule_spec) == \
        'A(loc~Cytosole,ModDomain1~P,AssociationDomain1)'


def test_string_from_molecule_reactant(simple_molecule_spec):
    assert bex.string_from_molecule_reactant(rbm.MoleculeReactant(simple_molecule_spec)) == \
        'A(loc~Cytosole,ModDomain1~P,AssociationDomain1)'


def test_string_from_complex_reactant(simple_molecule_spec):
    simple_molecule_spec.association_specs[0].occupation_status = rbm.OccupationStatus.occupied_known_partner
    bound_molecule_a = simple_molecule_spec

    bound_molecule_b = rbm.MoleculeSpecification(
            rbm.MoleculeDefinition(
                'B', [], [rbm.AssociationDefinition('AssociationDomain')], rbm.LocalizationDefinition(['Cytosole', 'Nucleus'])
            ),
            [],
            [rbm.AssociationSpecification(rbm.AssociationDefinition('AssociationDomain'), rbm.OccupationStatus.occupied_known_partner)],
            rbm.LocalizationSpecification(rbm.LocalizationDefinition(['Cytosole', 'Nucleus']), 'Cytosole')
    )

    complex_reactant = rbm.ComplexReactant(
            [bound_molecule_a, bound_molecule_b],
            [rbm.Binding(left_partner=(0, bound_molecule_a.association_specs[0]), right_partner=(1, bound_molecule_b.association_specs[0]))]
    )

    assert bex.string_from_complex_reactant(complex_reactant) == \
        'A(loc~Cytosole,ModDomain1~P,AssociationDomain1!0).B(loc~Cytosole,AssociationDomain!0)'


def test_simple_rule_based_model(simple_rule_based_model):
    bngl = bex.BNGLSystem(simple_rule_based_model, bex.BNGLSettings())

    expected_string = '''begin model
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
simulate({method=>"ode",t_end=>10,n_steps=>100})'''

    bngl.to_string()
    assert bngl.to_string() == expected_string


@pytest.fixture
def simple_molecule_spec():
    mod_defs  = [rbm.ModificationDefinition('ModDomain1', ['U', 'P']),
                 rbm.ModificationDefinition('ModDomain2', ['U', 'GTP'])]
    mod_specs = [rbm.ModificationSpecification(mod_defs[0], 'P')]

    ass_defs  = [rbm.AssociationDefinition('AssociationDomain1'),
                 rbm.AssociationDefinition('AssociationDomain2')]
    ass_specs = [rbm.AssociationSpecification(ass_defs[0], rbm.OccupationStatus.unoccupied)]

    loc_def   = rbm.LocalizationDefinition(['Cytosole', 'Nucleus'])
    loc_spec  = rbm.LocalizationSpecification(loc_def, 'Cytosole')

    molecule_definition = rbm.MoleculeDefinition('A', mod_defs, ass_defs, loc_def)

    molecule_specification = rbm.MoleculeSpecification(molecule_definition, mod_specs, ass_specs, loc_spec)

    return molecule_specification

@pytest.fixture
def simple_rule_based_model():
    assoc_def_a = rbm.AssociationDefinition('b')
    assoc_def_b = rbm.AssociationDefinition('a')

    mol_def_a = rbm.MoleculeDefinition('A', [], [assoc_def_a], None)
    mol_def_b = rbm.MoleculeDefinition('B', [], [assoc_def_b], None)

    mol_spec_a_unbound = rbm.MoleculeSpecification(mol_def_a, [], [rbm.AssociationSpecification(assoc_def_a, rbm.OccupationStatus.unoccupied)],
                                                   None)
    mol_spec_b_unbound = rbm.MoleculeSpecification(mol_def_b, [], [rbm.AssociationSpecification(assoc_def_b, rbm.OccupationStatus.unoccupied)],
                                                   None)

    mol_spec_a_bound = rbm.MoleculeSpecification(mol_def_a, [], [rbm.AssociationSpecification(assoc_def_a, rbm.OccupationStatus.occupied_known_partner)],
                                                 None)
    mol_spec_b_bound = rbm.MoleculeSpecification(mol_def_b, [], [rbm.AssociationSpecification(assoc_def_b, rbm.OccupationStatus.occupied_known_partner)],
                                                 None)

    reactant_a_unbound = rbm.MoleculeReactant(mol_spec_a_unbound)
    reactant_b_unbound = rbm.MoleculeReactant(mol_spec_b_unbound)

    binding = rbm.Binding((0, mol_spec_a_bound.association_specs[0]),
                          (1, mol_spec_b_bound.association_specs[0]))

    reactant_ab = rbm.ComplexReactant([mol_spec_a_bound, mol_spec_b_bound], [binding])

    rates = [rbm.Parameter('kf', '1e6/(NA*V)'), rbm.Parameter('kr', '0.1')]
    rule = rbm.Rule([reactant_a_unbound, reactant_b_unbound], [reactant_ab], rbm.Arrow.reversible, rates)

    parameters = [
        rbm.Parameter('A_total', '100'),
        rbm.Parameter('B_total', '100'),
        rbm.Parameter('NA', '6.02214e23'),
        rbm.Parameter('V', '1e-12'),
    ]

    parameters.extend(rates)

    initial_conditions = [
        rbm.InitialCondition(mol_spec_a_unbound, 'A_total'),
        rbm.InitialCondition(mol_spec_b_unbound, 'B_total')
    ]

    return rbm.RuleBasedModel([mol_def_a, mol_def_b], [rule], parameters, initial_conditions)








