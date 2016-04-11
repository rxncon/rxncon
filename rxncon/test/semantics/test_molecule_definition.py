import pytest
from rxncon.simulation.rule_based.molecule_from_string import mol_def_from_string

@pytest.fixture
def molecule_definitions():
    return [
        [mol_def_from_string('C#'), mol_def_from_string('A#'), mol_def_from_string('B#')],
        [mol_def_from_string('C#ass/C_[a]:A[z]'), mol_def_from_string('A#ass/A_[z]:C[a]~B_[a]'),

         mol_def_from_string('B#ass/B_[a]:A[z]')],
        [mol_def_from_string('CmRNA#'), mol_def_from_string('B#'), mol_def_from_string('A#')]

    ]

@pytest.fixture
def expected_molecule_definition_order():
    return [
        [mol_def_from_string('A#'), mol_def_from_string('B#'), mol_def_from_string('C#')],
        [mol_def_from_string('A#ass/A_[z]:C[a]~B_[a]'), mol_def_from_string('B#ass/B_[a]:A[z]'),

         mol_def_from_string('C#ass/C_[a]:A[z]')],
        [mol_def_from_string('CmRNA#'), mol_def_from_string('A#'), mol_def_from_string('B#')]

    ]

def test_molecule_definition_ordering(molecule_definitions, expected_molecule_definition_order):
    for i, mol_def in enumerate(molecule_definitions):
        assert sorted(mol_def) == expected_molecule_definition_order[i]


def test_molecule_definition_wrongly_defined():

    assert mol_def_from_string('A#ass/A_[d]:B_[Aassoc],mod/A_[d]:u~ub') != mol_def_from_string('A#')
    assert mol_def_from_string('A#ass/A_[d]:B_[Aassoc]~C_[Aassoc],mod/A_[d]:u~ub~p')

    with pytest.raises(AssertionError):
        mol_def_from_string('A#ass/A_[d]:B_[Aassoc],ass/A_[d]:C_[Aassoc]')
    with pytest.raises(AssertionError):
        mol_def_from_string('A#mod/A_[d]:u~p,mod/A_[d]:u~ub')

