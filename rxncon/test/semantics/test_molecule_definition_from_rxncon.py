import rxncon.core.contingency as con
import rxncon.core.effector as eff
import rxncon.core.rxncon_system as rxs
import rxncon.core.specification as spe
import rxncon.semantics.molecule_definition as mol
import rxncon.semantics.molecule_definition_from_rxncon as mdr
import rxncon.syntax.rxncon_from_string as rfs
from rxncon.simulation.rule_based.molecule_from_string import mol_def_from_string
from rxncon.input.quick.quick import Quick


def test_molecule_definitions_single_reaction():
    rxncon = rxs.RxnConSystem([rfs.reaction_from_string('A_ppi_B')], [])

    actual_mol_defs = mdr.mol_defs_from_rxncon_sys(rxncon)

    expected_mol_defs = {
        rfs.specification_from_string('A'): mol_def_from_string('A#ass/A_[Bassoc]:B_[Aassoc]'),
        rfs.specification_from_string('B'): mol_def_from_string('B#ass/B_[Aassoc]:A_[Bassoc]')
    }

    for spec, mol_def in expected_mol_defs.items():
        assert mol_def == actual_mol_defs[spec]


def test_molecule_definitions_no_contingencies():
    rxncon = rxs.RxnConSystem(
        [
            rfs.reaction_from_string('A_ppi_B'),
            rfs.reaction_from_string('A_ppi_C'),
            rfs.reaction_from_string('B_ppi_E'),
            rfs.reaction_from_string('B_p+_E')
        ], [])

    actual_mol_defs = mdr.mol_defs_from_rxncon_sys(rxncon)

    expected_mol_defs = {
        rfs.specification_from_string('A'): mol_def_from_string('A#ass/A_[Bassoc]:B_[Aassoc],ass/A_[Cassoc]:C_[Aassoc]'),
        rfs.specification_from_string('B'): mol_def_from_string('B#ass/B_[Eassoc]:E_[Bassoc],ass/B_[Aassoc]:A_[Bassoc]'),
        rfs.specification_from_string('C'): mol_def_from_string('C#ass/C_[Aassoc]:A_[Cassoc]'),
        rfs.specification_from_string('E'): mol_def_from_string('E#mod/E_[(Bsite)]:u~p,ass/E_[Bassoc]:B_[Eassoc]')
    }

    for spec, mol_def in expected_mol_defs.items():
        assert mol_def == actual_mol_defs[spec]


def test_molecule_definitions_with_contingencies():
    rxncon = Quick('''A_ppi_B; x A--C
    A_ppi_C
    B_ppi_E; ! A--B
    B_p+_E; ! B--E''').rxncon_system

    actual_mol_defs = mdr.mol_defs_from_rxncon_sys(rxncon)

    expected_mol_defs = {
        rfs.specification_from_string('A'): mol_def_from_string('A#ass/A_[Bassoc]:B_[Aassoc],ass/A_[Cassoc]:C_[Aassoc]'),
        rfs.specification_from_string('B'): mol_def_from_string('B#ass/B_[Eassoc]:E_[Bassoc],ass/B_[Aassoc]:A_[Bassoc]'),
        rfs.specification_from_string('C'): mol_def_from_string('C#ass/C_[Aassoc]:A_[Cassoc]'),
        rfs.specification_from_string('E'): mol_def_from_string('E#mod/E_[(Bsite)]:u~p,ass/E_[Bassoc]:B_[Eassoc]')
    }

    for spec, mol_def in expected_mol_defs.items():
        assert mol_def == actual_mol_defs[spec]


def test_molecule_definitions_multiple_kinases_same_modification_at_same_residue():
    rxncon = rxs.RxnConSystem([rfs.reaction_from_string('A_p+_B_[(x)]'), rfs.reaction_from_string('C_p+_B_[(x)]')], [])

    actual_mol_defs = mdr.mol_defs_from_rxncon_sys(rxncon)

    expected_mol_defs = {
        rfs.specification_from_string('A'): mol_def_from_string('A#'),
        rfs.specification_from_string('B'): mol_def_from_string('B#mod/B_[(x)]:u~p'),
        rfs.specification_from_string('C'): mol_def_from_string('C#')
    }

    for spec, mol_def in expected_mol_defs.items():
        assert mol_def == actual_mol_defs[spec]


def test_molecule_definitions_multiple_kinases_different_modifications_at_same_residue():
    rxncon = rxs.RxnConSystem([rfs.reaction_from_string('A_p+_B_[(x)]'), rfs.reaction_from_string('C_ub+_B_[(x)]')], [])

    actual_mol_defs = mdr.mol_defs_from_rxncon_sys(rxncon)

    expected_mol_defs = {
        rfs.specification_from_string('A'): mol_def_from_string('A#'),
        rfs.specification_from_string('B'): mol_def_from_string('B#mod/B_[(x)]:u~p~ub'),
        rfs.specification_from_string('C'): mol_def_from_string('C#')
    }

    for spec, mol_def in expected_mol_defs.items():
        assert mol_def == actual_mol_defs[spec]


def test_molecule_definitions_multiple_partners_binding_same_domain():
    rxncon = rxs.RxnConSystem([rfs.reaction_from_string('A_ppi_B_[x]'), rfs.reaction_from_string('C_ppi_B_[x]')], [])

    actual_mol_defs = mdr.mol_defs_from_rxncon_sys(rxncon)

    expected_mol_defs = {
        rfs.specification_from_string('A'): mol_def_from_string('A#ass/A_[Bassoc]:B_[x]'),
        rfs.specification_from_string('B'): mol_def_from_string('B#ass/B_[x]:A_[Bassoc]~C_[Bassoc]'),
        rfs.specification_from_string('C'): mol_def_from_string('C#ass/C_[Bassoc]:B_[x]')
    }

    for spec, mol_def in expected_mol_defs.items():
        assert mol_def == actual_mol_defs[spec]
