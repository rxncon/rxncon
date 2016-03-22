import pytest

from rxncon.semantics.molecule_from_string import mol_def_from_string
from rxncon.semantics.molecule_definition import MoleculeDefinition, ModificationPropertyDefinition, \
    AssociationPropertyDefinition, LocalizationPropertyDefinition, Compartment, Modifier
from rxncon.syntax.rxncon_from_string import component_from_string


def test_mol_def_from_string(mol_defs):
    for mol_def_string, mol_def in mol_defs.items():
        assert mol_def_from_string(mol_def_string) == mol_def


@pytest.fixture
def mol_defs():
    return {
        'X#': MoleculeDefinition(component_from_string('X'), set(), set(), None),
        'A#ass/A_[x]:B_[y]~C_[z]': MoleculeDefinition(
            component_from_string('A'),
            set(),
            {AssociationPropertyDefinition(component_from_string('A_[x]'),
                                           {component_from_string('B_[y]'),
                                            component_from_string('C_[z]')})},
            None
        ),
        'A#loc/cell~nucleus': MoleculeDefinition(
            component_from_string('A'),
            set(),
            set(),
            LocalizationPropertyDefinition({Compartment('nucleus'), Compartment('cell')})
        ),
        'A#mod/A_[d(r)]:u~ub~p': MoleculeDefinition(
            component_from_string('A'),
            {ModificationPropertyDefinition(component_from_string('A_[d(r)]'),
                                            {Modifier('ub'), Modifier('p'), Modifier('u')})},
            set(),
            None
        ),
        'B#mod/B_[dd(r1)]:u~p,ass/B_[z]:A_[x]': MoleculeDefinition(
            component_from_string('B'),
            {ModificationPropertyDefinition(component_from_string('B_[dd(r1)]'),
                                            {Modifier('p'), Modifier('u')})},
            {AssociationPropertyDefinition(component_from_string('B_[z]'),
                                           {component_from_string('A_[x]')})},
            None
        )
    }

