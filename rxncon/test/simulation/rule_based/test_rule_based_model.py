import pytest
import typecheck as tc

import rxncon.semantics.molecule
import rxncon.simulation.rule_based.rule_based_model as rbm


### BASIC BUILDING BLOCKS: modification, association, localization ###
def test_modification_definition():
    modification_definition = rxncon.semantics.molecule.ModificationDefinition('ModDomain', ['U', 'P'])
    assert modification_definition.spec == 'ModDomain'
    assert modification_definition.valid_modifiers == ['U', 'P']


def test_modification_specification():
    modification_definition = rxncon.semantics.molecule.ModificationDefinition('ModDomain', ['U', 'P'])
    modification_specification = rxncon.semantics.molecule.ModificationInstance(modification_definition, 'P')
    assert modification_specification.modifier == 'P'

    with pytest.raises(ValueError):
        rxncon.semantics.molecule.ModificationInstance(modification_definition, 'G')


def test_association_definition():
    association_definition = rxncon.semantics.molecule.AssociationDefinition('AssociationDomain')
    assert association_definition.spec == 'AssociationDomain'


def test_association_specification():
    association_definition = rxncon.semantics.molecule.AssociationDefinition('AssociationDomain')
    association_specification_occupied = rxncon.semantics.molecule.AssociationInstance(association_definition, rxncon.semantics.molecule.OccupationStatus.occupied_unknown_partner)
    assert association_specification_occupied

    # Occupation status was Boolean in previous version, this should raise an error.
    with pytest.raises(tc.InputParameterError):
        rxncon.semantics.molecule.AssociationInstance(association_definition, True)


def test_localization_definition():
    localization_definition = rxncon.semantics.molecule.LocalizationDefinition(['Cell', 'Cytoplasm', 'Nucleus'])
    assert localization_definition.valid_compartments == ['Cell', 'Cytoplasm', 'Nucleus']


def test_localization_specification():
    localization_definition = rxncon.semantics.molecule.LocalizationDefinition(['Cell', 'Cytoplasm', 'Nucleus'])
    localization_specification = rxncon.semantics.molecule.LocalizationInstance(localization_definition, 'Cell')
    assert localization_specification.compartment == 'Cell'

    with pytest.raises(ValueError):
        rxncon.semantics.molecule.LocalizationInstance(localization_definition, 'Burger King')


### MOLECULE DEFINITION / SPECIFICATION ###
def test_molecule_definition():
    mod_defs = [rxncon.semantics.molecule.ModificationDefinition('ModDomain1', ['U', 'P']), rxncon.semantics.molecule.ModificationDefinition('ModDomain2', ['U', 'GTP'])]
    ass_defs = [rxncon.semantics.molecule.AssociationDefinition('AssociationDomain1'), rxncon.semantics.molecule.AssociationDefinition('AssociationDomain2')]
    loc_def  = rxncon.semantics.molecule.LocalizationDefinition(['Cell', 'Cytoplasm', 'Nucleus'])

    molecule_definition = rxncon.semantics.molecule.MoleculeDefinition('A', mod_defs, ass_defs, loc_def)

    assert molecule_definition.spec == 'A'
    assert molecule_definition.modification_defs == mod_defs
    assert molecule_definition.association_defs == ass_defs
    assert molecule_definition.localization_def == loc_def


def test_molecule_specification(modification_defs, association_defs, localization_def):
    modification_specs = [rxncon.semantics.molecule.ModificationInstance(modification_defs[0], 'P')]
    association_specs  = [rxncon.semantics.molecule.AssociationInstance(association_defs[0], rxncon.semantics.molecule.OccupationStatus.occupied_unknown_partner)]
    localization_spec  = rxncon.semantics.molecule.LocalizationInstance(localization_def, 'Cytoplasm')

    molecule_def = rxncon.semantics.molecule.MoleculeDefinition('A', modification_defs, association_defs, localization_def)

    molecule_spec = rxncon.semantics.molecule.MoleculeInstance(molecule_def, modification_specs, association_specs, localization_spec)

    assert molecule_spec.molecule_def == molecule_def


### REACTANTS: MOLECULEREACTANT, COMPLEXREACTANT ###
def test_molecule_reactant(modification_defs, association_defs, localization_def):
    modification_specs = [rxncon.semantics.molecule.ModificationInstance(modification_defs[0], 'P')]
    association_specs  = [rxncon.semantics.molecule.AssociationInstance(association_defs[0], rxncon.semantics.molecule.OccupationStatus.not_occupied)]
    localization_spec  = rxncon.semantics.molecule.LocalizationInstance(localization_def, 'Cytoplasm')

    molecule_definition = rxncon.semantics.molecule.MoleculeDefinition('A', modification_defs, association_defs, localization_def)
    molecule_specification = rxncon.semantics.molecule.MoleculeInstance(molecule_definition, modification_specs, association_specs, localization_spec)
    molecule_reactant = rbm.MoleculeReactant(molecule_specification)

    assert molecule_reactant.molecule_specification == molecule_specification


def test_binding_valid():
    association_definitions_left = rxncon.semantics.molecule.AssociationDefinition('AssociationDomain1')
    association_specifications_left = rxncon.semantics.molecule.AssociationInstance(association_definitions_left,
                                                                                    rxncon.semantics.molecule.OccupationStatus.occupied_known_partner)

    association_definitions_right = rxncon.semantics.molecule.AssociationDefinition('AssociationDomain2')
    association_specifications_right = rxncon.semantics.molecule.AssociationInstance(association_definitions_right,
                                                                                     rxncon.semantics.molecule.OccupationStatus.occupied_known_partner)

    binding = rxncon.semantics.molecule.Binding(left_partner=(0, association_specifications_left),
                                                right_partner=(1, association_specifications_right))

    assert binding.left_partner[0] == 0
    assert binding.left_partner[1].occupation_status == rxncon.semantics.molecule.OccupationStatus.occupied_known_partner

    assert binding.right_partner[0] == 1
    assert binding.right_partner[1].occupation_status == rxncon.semantics.molecule.OccupationStatus.occupied_known_partner


def test_binding_raises_if_index_not_unique():
    association_definitions_left = rxncon.semantics.molecule.AssociationDefinition('AssociationDomain1')
    association_specifications_left = rxncon.semantics.molecule.AssociationInstance(association_definitions_left,
                                                                                    rxncon.semantics.molecule.OccupationStatus.occupied_known_partner)

    association_definitions_right = rxncon.semantics.molecule.AssociationDefinition('AssociationDomain2')
    association_specifications_right = rxncon.semantics.molecule.AssociationInstance(association_definitions_right,
                                                                                     rxncon.semantics.molecule.OccupationStatus.occupied_known_partner)

    with pytest.raises(ValueError):
        binding = rxncon.semantics.molecule.Binding(left_partner=(0, association_specifications_left),
                                                    right_partner=(0, association_specifications_right))


def test_complex_reactant(molecules_bound):
    binding_list = [rxncon.semantics.molecule.Binding(left_partner=(0, molecules_bound[0].association_specs[0]),
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
    binding_list = [rxncon.semantics.molecule.Binding(left_partner=(0, molecules_bound[0].association_specs[0]),
                                                      right_partner=(1, molecules_bound[1].association_specs[0]))]

    localization_spec_invalid = rxncon.semantics.molecule.LocalizationInstance(localization_def, 'Nucleus')

    molecules_bound[1] = rxncon.semantics.molecule.MoleculeInstance(molecules_bound[1].molecule_def, molecules_bound[1].modification_specs,
                                                                    molecules_bound[1].association_specs, localization_spec_invalid)

    with pytest.raises(ValueError):
        rbm.ComplexReactant(molecules_bound, binding_list)


### RULES ###
@pytest.fixture
def association_defs():
    return [rxncon.semantics.molecule.AssociationDefinition('AssociationDomain1'), rxncon.semantics.molecule.AssociationDefinition('AssociationDomain2')]


@pytest.fixture
def localization_def():
    return rxncon.semantics.molecule.LocalizationDefinition(['Cell', 'Cytoplasm', 'Nucleus'])


@pytest.fixture
def molecules_bound(modification_defs, association_defs, localization_def):
    modification_specifications_A = [rxncon.semantics.molecule.ModificationInstance(modification_defs[0], 'P')]

    association_specifications_A = [rxncon.semantics.molecule.AssociationInstance(association_defs[0], rxncon.semantics.molecule.OccupationStatus.occupied_known_partner)]
    association_specifications_B = [rxncon.semantics.molecule.AssociationInstance(association_defs[1], rxncon.semantics.molecule.OccupationStatus.occupied_known_partner)]

    localization_specification_A = rxncon.semantics.molecule.LocalizationInstance(localization_def, 'Cytoplasm')
    localization_specification_B = rxncon.semantics.molecule.LocalizationInstance(localization_def, 'Cytoplasm')

    molecule_definition = rxncon.semantics.molecule.MoleculeDefinition('A', modification_defs, association_defs, localization_def)

    molecule_specification_A = rxncon.semantics.molecule.MoleculeInstance(molecule_definition, modification_specifications_A,
                                                                          association_specifications_A, localization_specification_A)

    molecule_specification_B = rxncon.semantics.molecule.MoleculeInstance(molecule_definition, [], association_specifications_B,
                                                                          localization_specification_B)

    return [molecule_specification_A, molecule_specification_B]


@pytest.fixture
def molecules_unbound(molecules_bound):
    molecule_specification_A_bound = molecules_bound[0]
    molecule_specification_B_bound = molecules_bound[1]

    association_specification_A_unbound = [
        rxncon.semantics.molecule.AssociationInstance(molecule_specification_A_bound.association_specifications[0],
                                                      rxncon.semantics.molecule.OccupationStatus.not_occupied)]
    assocociation_specification_B_unbound = [
        rxncon.semantics.molecule.AssociationInstance(molecule_specification_B_bound.association_specifications[0],
                                                      rxncon.semantics.molecule.OccupationStatus.not_occupied)]

    molecule_definition_A = rxncon.semantics.molecule.MoleculeDefinition('A', modification_defs, association_defs, localization_def)
    molecule_specification_A_unbound = rxncon.semantics.molecule.MoleculeInstance(molecule_definition_A,
                                                                                  molecule_specification_A_bound.modification_specs,
                                                                                  association_specification_A_unbound,
                                                                                  molecule_specification_A_bound.localization_spec)

    molecule_definition_B = rxncon.semantics.molecule.MoleculeDefinition('B', [], association_defs, localization_def)
    molecule_specification_B_unbound = rxncon.semantics.molecule.MoleculeInstance(molecule_definition_B,
                                                                                  [],
                                                                                  assocociation_specification_B_unbound,
                                                                                  molecule_specification_B_bound.localization_spec)

    return [molecule_specification_A_unbound, molecule_specification_B_unbound]


@pytest.fixture
def modification_defs():
    return [rxncon.semantics.molecule.ModificationDefinition('ModDomain1', ['U', 'P']), rxncon.semantics.molecule.ModificationDefinition('ModDomain2', ['U', 'GTP'])]
