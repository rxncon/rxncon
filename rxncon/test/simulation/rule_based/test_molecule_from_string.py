import pytest

from rxncon.simulation.rule_based.molecule_from_string import mol_def_from_string, mol_instance_from_string, \
    mol_instances_and_bindings_from_string, rule_from_string
from rxncon.semantics.molecule_definition import MoleculeDefinition, ModificationPropertyDefinition, \
    AssociationPropertyDefinition, LocalizationPropertyDefinition, Compartment, OccupationStatus
from rxncon.semantics.molecule import Modifier, OccupationStatus, Compartment, MoleculeDefinition, \
    ModificationPropertyDefinition, AssociationPropertyDefinition, LocalizationPropertyDefinition, Molecule, \
    ModificationProperty, AssociationProperty
from rxncon.semantics.molecule_instance import MoleculeInstance, AssociationPropertyInstance, ModificationPropertyInstance
from rxncon.syntax.rxncon_from_string import specification_from_string


def test_mol_def_from_string(mol_defs):
    for mol_def_string, mol_def in mol_defs.items():
        assert mol_def_from_string(mol_def_string) == mol_def


def test_mol_ins_from_string(mol_instances):
    for strings, mol_ins in mol_instances.items():
        def_string = strings[0]
        ins_string = strings[1]

        assert mol_instance_from_string(def_string, ins_string) == mol_ins


def test_mol_instances_and_binding():
    mol_defs = ['A#ass/A_[x]:B_[y]', 'B#ass/B_[y]:A_[x]']
    complex_string = 'A#ass/A_[x]:B_[y]~1.B#ass/B_[y]:A_[x]~1'

    instances, bindings = mol_instances_and_bindings_from_string(mol_defs, complex_string)

    assert mol_instance_from_string(mol_defs[0], 'A#ass/A_[x]:B_[y]') in instances
    assert mol_instance_from_string(mol_defs[1], 'B#ass/B_[y]:A_[x]') in instances

    assert len(bindings) == 1
    binding = bindings[0]
    assert binding.left_partner[0] == 0
    assert binding.right_partner[0] == 1
    assert binding.left_partner[1].association_def.spec == specification_from_string('A_[x]')
    assert binding.right_partner[1].association_def.spec == specification_from_string('B_[y]')


def test_rule():
    mol_defs = ['A#ass/A_[x]:B_[y]~C_[z],mod/A_[(r)]:u~p', 'B#ass/B_[y]:A_[x]']
    rule_string = 'A#ass/A_[x]:,mod/A_[(r)]:p + B#ass/B_[y]: <-> A#ass/A_[x]:B_[y]~0,mod/A_[(r)]:p.B#ass/B_[y]:A_[x]~0'

    rule = rule_from_string(mol_defs, rule_string)
    print(rule)


@pytest.fixture
def mol_defs():
    return {
        'X#': MoleculeDefinition(specification_from_string('X'), set(), set(), None),
        'A#ass/A_[x]:B_[y]~C_[z]': MoleculeDefinition(
            specification_from_string('A'),
            set(),
            {AssociationPropertyDefinition(specification_from_string('A_[x]'),
                                           {specification_from_string('B_[y]'),
                                            specification_from_string('C_[z]')})},
            None
        ),
        'A#loc/cell~nucleus': MoleculeDefinition(
            specification_from_string('A'),
            set(),
            set(),
            LocalizationPropertyDefinition({Compartment('nucleus'), Compartment('cell')})
        ),
        'A#mod/A_[d(r)]:u~ub~p': MoleculeDefinition(
            specification_from_string('A'),
            {ModificationPropertyDefinition(specification_from_string('A_[d(r)]'),
                                            {Modifier('ub'), Modifier('p'), Modifier('u')})},
            set(),
            None
        ),
        'B#mod/B_[dd(r1)]:u~p,ass/B_[z]:A_[x]': MoleculeDefinition(
            specification_from_string('B'),
            {ModificationPropertyDefinition(specification_from_string('B_[dd(r1)]'),
                                            {Modifier('p'), Modifier('u')})},
            {AssociationPropertyDefinition(specification_from_string('B_[z]'),
                                           {specification_from_string('A_[x]')})},
            None
        )
    }


@pytest.fixture
def mol_instances():
    mol_instances = {}

    ass_def = list(mol_def_from_string('A#ass/A_[x]:B_[y]~C_[z]').association_defs)[0]

    mol_instances[('A#ass/A_[x]:B_[y]~C_[z]', 'A#ass/A_[x]:')] = Molecule(
        mol_def_from_string('A#ass/A_[x]:B_[y]~C_[z]'),
        set(),
        {AssociationProperty(ass_def, OccupationStatus.not_occupied, None)},
        None
    )

    ass_def = list(mol_def_from_string('A#ass/A_[x]:B_[y]~C_[z]').association_defs)[0]

    mol_instances[('A#ass/A_[x]:B_[y]~C_[z]', 'A#ass/A_[x]:C_[z]')] = Molecule(
        mol_def_from_string('A#ass/A_[x]:B_[y]~C_[z]'),
        set(),
        {AssociationProperty(ass_def, OccupationStatus.occupied_known_partner, specification_from_string('C_[z]'))},
        None
    )

    mod_def = list(mol_def_from_string('A#mod/A_[(r)]:u~p~ub').modification_defs)[0]

    mol_instances[('A#mod/A_[(r)]:u~p~ub', 'A#mod/A_[(r)]:p')] = Molecule(
        mol_def_from_string('A#mod/A_[(r)]:u~p~ub'),
        {ModificationProperty(mod_def, Modifier.phosphorylated)},
        set(),
        None
    )

    return mol_instances


