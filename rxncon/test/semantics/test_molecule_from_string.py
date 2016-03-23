import pytest

from rxncon.semantics.molecule_from_string import mol_def_from_string, mol_instance_from_string
from rxncon.semantics.molecule_definition import MoleculeDefinition, ModificationPropertyDefinition, \
    AssociationPropertyDefinition, LocalizationPropertyDefinition, Compartment, Modifier, OccupationStatus
from rxncon.semantics.molecule_instance import MoleculeInstance, AssociationPropertyInstance
from rxncon.syntax.rxncon_from_string import component_from_string


def test_mol_def_from_string(mol_defs):
    for mol_def_string, mol_def in mol_defs.items():
        assert mol_def_from_string(mol_def_string) == mol_def


def test_mol_ins_from_string(mol_instances):
    for strings, mol_ins in mol_instances.items():
        def_string = strings[0]
        ins_string = strings[1]

        assert mol_instance_from_string(def_string, ins_string) == mol_ins


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


@pytest.fixture
def mol_instances():
    mol_instances = {}

    ass_def = list(mol_def_from_string('A#ass/A_[x]:B_[y]~C_[z]').association_defs)[0]

    mol_instances[('A#ass/A_[x]:B_[y]~C_[z]', 'A#ass/A_[x]:')] = MoleculeInstance(
        mol_def_from_string('A#ass/A_[x]:B_[y]~C_[z]'),
        set(),
        {AssociationPropertyInstance(ass_def, OccupationStatus.not_occupied, None)},
        None
    )

    ass_def = list(mol_def_from_string('A#ass/A_[x]:B_[y]~C_[z]').association_defs)[0]

    mol_instances[('A#ass/A_[x]:B_[y]~C_[z]', 'A#ass/A_[x]:C_[z]')] = MoleculeInstance(
        mol_def_from_string('A#ass/A_[x]:B_[y]~C_[z]'),
        set(),
        {AssociationPropertyInstance(ass_def, OccupationStatus.occupied_known_partner, component_from_string('C_[z]'))},
        None
    )

    return mol_instances
