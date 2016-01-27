import pytest

import rxncon.simulation.rule_based.rule_based_model as rbm


### BASIC BUILDING BLOCKS: modification, association, localization ###
def test_modification_definition():
    modification_definition = rbm.ModificationDefinition('ModDomain', ['U', 'P'])
    assert modification_definition.domain == 'ModDomain'
    assert modification_definition.valid_modifiers == ['U', 'P']


def test_modification_specification():
    modification_definition = rbm.ModificationDefinition('ModDomain', ['U', 'P'])
    modification_specification = rbm.ModificationSpecification(modification_definition, 'P')
    assert modification_specification.modifier == 'P'

    with pytest.raises(ValueError):
        rbm.ModificationSpecification(modification_definition, 'G')


def test_association_definition():
    association_definition = rbm.AssociationDefinition('AssociationDomain')
    assert association_definition.domain == 'AssociationDomain'


def test_association_specification():
    association_definition = rbm.AssociationDefinition('AssociationDomain')
    association_specification_occupied = rbm.AssociationSpecification(association_definition, rbm.OccupationStatus.occupied_unknown_partner)
    assert association_specification_occupied

    # Occupation status was Boolean in previous version, this should raise an error.
    with pytest.raises(AssertionError):
        rbm.AssociationSpecification(association_definition, True)

def test_localization_definition():
    localization_definition = rbm.LocalizationDefinition(['Cell', 'Cytoplasm', 'Nucleus'])
    assert localization_definition.valid_compartments == ['Cell', 'Cytoplasm', 'Nucleus']


def test_localization_specification():
    localization_definition = rbm.LocalizationDefinition(['Cell', 'Cytoplasm', 'Nucleus'])
    localization_specification = rbm.LocalizationSpecification(localization_definition, 'Cell')
    assert localization_specification.compartment == 'Cell'

    with pytest.raises(ValueError):
        rbm.LocalizationSpecification(localization_definition, 'Burger King')


### MOLECULE DEFINITION / SPECIFICATION ###
def test_molecule_definition():
    mod_defs = [rbm.ModificationDefinition('ModDomain1', ['U', 'P']), rbm.ModificationDefinition('ModDomain2', ['U', 'GTP'])]
    ass_defs = [rbm.AssociationDefinition('AssociationDomain1'), rbm.AssociationDefinition('AssociationDomain2')]
    loc_def  = rbm.LocalizationDefinition(['Cell', 'Cytoplasm', 'Nucleus'])

    molecule_definition = rbm.MoleculeDefinition('A', mod_defs, ass_defs, loc_def)

    assert molecule_definition.name == 'A'
    assert molecule_definition.modification_defs == mod_defs
    assert molecule_definition.association_defs == ass_defs
    assert molecule_definition.localization_def == loc_def


def test_molecule_specification(modification_defs, association_defs, localization_def):
    modification_specs = [rbm.ModificationSpecification(modification_defs[0], 'P')]
    association_specs  = [rbm.AssociationSpecification(association_defs[0], rbm.OccupationStatus.occupied_unknown_partner)]
    localization_spec  = rbm.LocalizationSpecification(localization_def, 'Cytoplasm')

    molecule_def = rbm.MoleculeDefinition('A', modification_defs, association_defs, localization_def)

    molecule_spec = rbm.MoleculeSpecification(molecule_def, modification_specs, association_specs, localization_spec)

    assert molecule_spec.molecule_def == molecule_def


### REACTANTS: MOLECULEREACTANT, COMPLEXREACTANT ###
def test_molecule_reactant(modification_defs, association_defs, localization_def):
    modification_specs = [rbm.ModificationSpecification(modification_defs[0], 'P')]
    association_specs  = [rbm.AssociationSpecification(association_defs[0], rbm.OccupationStatus.not_unoccupied)]
    localization_spec  = rbm.LocalizationSpecification(localization_def, 'Cytoplasm')

    molecule_definition = rbm.MoleculeDefinition('A', modification_defs, association_defs, localization_def)
    molecule_specification = rbm.MoleculeSpecification(molecule_definition, modification_specs, association_specs, localization_spec)
    molecule_reactant = rbm.MoleculeReactant(molecule_specification)

    assert molecule_reactant.molecule_specification == molecule_specification


def test_binding_valid():
    association_definitions_left = rbm.AssociationDefinition('AssociationDomain1')
    association_specifications_left = rbm.AssociationSpecification(association_definitions_left,
                                                                   rbm.OccupationStatus.occupied_known_partner)

    association_definitions_right = rbm.AssociationDefinition('AssociationDomain2')
    association_specifications_right = rbm.AssociationSpecification(association_definitions_right,
                                                                    rbm.OccupationStatus.occupied_known_partner)

    binding = rbm.Binding(left_partner=(0, association_specifications_left),
                          right_partner=(1, association_specifications_right))

    assert binding.left_partner[0] == 0
    assert binding.left_partner[1].occupation_status == rbm.OccupationStatus.occupied_known_partner

    assert binding.right_partner[0] == 1
    assert binding.right_partner[1].occupation_status == rbm.OccupationStatus.occupied_known_partner


def test_binding_raises_if_index_not_unique():
    association_definitions_left = rbm.AssociationDefinition('AssociationDomain1')
    association_specifications_left = rbm.AssociationSpecification(association_definitions_left,
                                                                   rbm.OccupationStatus.occupied_known_partner)

    association_definitions_right = rbm.AssociationDefinition('AssociationDomain2')
    association_specifications_right = rbm.AssociationSpecification(association_definitions_right,
                                                                    rbm.OccupationStatus.occupied_known_partner)

    with pytest.raises(ValueError):
        binding = rbm.Binding(left_partner=(0, association_specifications_left),
                              right_partner=(0, association_specifications_right))


def test_complex_reactant(molecules_bound):
    binding_list = [rbm.Binding(left_partner=(0, molecules_bound[0].association_specs[0]),
                                right_partner=(1, molecules_bound[1].association_specs[0]))]
    complex_reactant = rbm.ComplexReactant(molecules_bound, binding_list)

    assert len(complex_reactant.molecules) == 2

    assert complex_reactant.molecules[0].modification_specs == molecules_bound[0].modification_specs
    assert complex_reactant.molecules[1].modification_specs == []

    assert complex_reactant.molecules[0].association_specs == molecules_bound[0].association_specs
    assert complex_reactant.molecules[1].association_specs == molecules_bound[1].association_specs

    assert complex_reactant.molecules[0].localization_spec == molecules_bound[0].localization_spec
    assert complex_reactant.molecules[1].localization_spec == molecules_bound[1].localization_spec


def test_complex_reactant_raises_if_parts_in_different_location(molecules_bound, localization_def):
    binding_list = [rbm.Binding(left_partner=(0, molecules_bound[0].association_specs[0]),
                                right_partner=(1, molecules_bound[1].association_specs[0]))]

    localization_spec_invalid = rbm.LocalizationSpecification(localization_def, 'Nucleus')

    molecules_bound[1] = rbm.MoleculeSpecification(molecules_bound[1].molecule_def, molecules_bound[1].modification_specs,
                                                   molecules_bound[1].association_specs, localization_spec_invalid)

    with pytest.raises(ValueError):
        rbm.ComplexReactant(molecules_bound, binding_list)


### RULES ###
@pytest.fixture
def association_defs():
    return [rbm.AssociationDefinition('AssociationDomain1'), rbm.AssociationDefinition('AssociationDomain2')]


@pytest.fixture
def localization_def():
    return rbm.LocalizationDefinition(['Cell', 'Cytoplasm', 'Nucleus'])


@pytest.fixture
def molecules_bound(modification_defs, association_defs, localization_def):
    modification_specifications_A = [rbm.ModificationSpecification(modification_defs[0], 'P')]

    association_specifications_A = [rbm.AssociationSpecification(association_defs[0], rbm.OccupationStatus.occupied_known_partner)]
    association_specifications_B = [rbm.AssociationSpecification(association_defs[1], rbm.OccupationStatus.occupied_known_partner)]

    localization_specification_A = rbm.LocalizationSpecification(localization_def, 'Cytoplasm')
    localization_specification_B = rbm.LocalizationSpecification(localization_def, 'Cytoplasm')

    molecule_definition = rbm.MoleculeDefinition('A', modification_defs, association_defs, localization_def)

    molecule_specification_A = rbm.MoleculeSpecification(molecule_definition, modification_specifications_A,
                                                         association_specifications_A, localization_specification_A)

    molecule_specification_B = rbm.MoleculeSpecification(molecule_definition, [], association_specifications_B,
                                                         localization_specification_B)

    return [molecule_specification_A, molecule_specification_B]


@pytest.fixture
def molecules_unbound(molecules_bound):
    molecule_specification_A_bound = molecules_bound[0]
    molecule_specification_B_bound = molecules_bound[1]

    association_specification_A_unbound = [rbm.AssociationSpecification(molecule_specification_A_bound.association_specifications[0],
                                                                        rbm.OccupationStatus.not_unoccupied)]
    assocociation_specification_B_unbound = [rbm.AssociationSpecification(molecule_specification_B_bound.association_specifications[0],
                                                                          rbm.OccupationStatus.not_unoccupied)]

    molecule_definition_A = rbm.MoleculeDefinition('A', modification_defs, association_defs, localization_def)
    molecule_specification_A_unbound = rbm.MoleculeSpecification(molecule_definition_A,
                                                                 molecule_specification_A_bound.modification_specs,
                                                                 association_specification_A_unbound,
                                                                 molecule_specification_A_bound.localization_spec)

    molecule_definition_B = rbm.MoleculeDefinition('B', [], association_defs, localization_def)
    molecule_specification_B_unbound = rbm.MoleculeSpecification(molecule_definition_B,
                                                                 [],
                                                                 assocociation_specification_B_unbound,
                                                                 molecule_specification_B_bound.localization_spec)

    return [molecule_specification_A_unbound, molecule_specification_B_unbound]


@pytest.fixture
def modification_defs():
    return [rbm.ModificationDefinition('ModDomain1', ['U', 'P']), rbm.ModificationDefinition('ModDomain2', ['U', 'GTP'])]
