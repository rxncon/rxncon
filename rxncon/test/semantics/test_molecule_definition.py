import pytest
from rxncon.simulation.rule_based.molecule_from_string import mol_def_from_string

@pytest.fixture
def molecule_definitions():
    return [
        [mol_def_from_string('C#'), mol_def_from_string('A#'), mol_def_from_string('B#')],
        [mol_def_from_string('C#ass/C_[a]:A[z]'), mol_def_from_string('A#ass/A_[z]:C[a]~B_[a]'),
         mol_def_from_string('B#ass/B_[a]:A[z]')]
    ]

@pytest.fixture
def expected_molecule_definition_order():
    return [
        [mol_def_from_string('A#'), mol_def_from_string('B#'), mol_def_from_string('C#')],
        [mol_def_from_string('A#ass/A_[z]:C[a]~B_[a]'), mol_def_from_string('B#ass/B_[a]:A[z]'),
         mol_def_from_string('C#ass/C_[a]:A[z]')]
    ]

def test_molecule_definition_ordering(molecule_definitions, expected_molecule_definition_order):
    for i, mol_def in enumerate(molecule_definitions):
        assert sorted(mol_def) == expected_molecule_definition_order[i]