import pytest
from collections import namedtuple

import rxncon.semantics.molecule_definition_from_rxncon as mdr
import rxncon.syntax.rxncon_from_string as rfs
from rxncon.simulation.rule_based.molecule_from_string import mol_def_from_string
from rxncon.input.quick.quick import Quick

MoleculeDefinitionTestCase = namedtuple('MoleculeDefinitionTestCase', ['quick_string', 'expected_mol_def'])


def test_molecule_defintion(the_case_molecule_definition):
    for the_case in the_case_molecule_definition:
        is_molecule_definition_correct(the_case)


def is_molecule_definition_correct(the_case):
    quick_system = Quick(the_case.quick_string)
    actual_mol_defs = mdr.mol_defs_from_rxncon_sys(quick_system.rxncon_system)

    assert all([expected_mol_def == actual_mol_defs[expected_spec] for expected_spec, expected_mol_def in the_case.expected_mol_def.items()])
    assert all([actual_mol_def == the_case.expected_mol_def[actual_spec] for actual_spec, actual_mol_def in actual_mol_defs.items()])


@pytest.fixture
def the_case_molecule_definition():
    return [
        MoleculeDefinitionTestCase("""A_ppi_B""",
                                   {rfs.specification_from_string('A'): mol_def_from_string('A#ass/A_[Bassoc]:B_[Aassoc]'),
                                   rfs.specification_from_string('B'): mol_def_from_string('B#ass/B_[Aassoc]:A_[Bassoc]')}
                                   ),

        MoleculeDefinitionTestCase("""A_ppi_B
                                      A_ppi_C
                                      B_ppi_E
                                      B_p+_E""",
                                   {rfs.specification_from_string('A'): mol_def_from_string('A#ass/A_[Bassoc]:B_[Aassoc],ass/A_[Cassoc]:C_[Aassoc]'),
                                    rfs.specification_from_string('B'): mol_def_from_string('B#ass/B_[Eassoc]:E_[Bassoc],ass/B_[Aassoc]:A_[Bassoc]'),
                                    rfs.specification_from_string('C'): mol_def_from_string('C#ass/C_[Aassoc]:A_[Cassoc]'),
                                    rfs.specification_from_string('E'): mol_def_from_string('E#mod/E_[(Bsite)]:u~p,ass/E_[Bassoc]:B_[Eassoc]')}
                                   ),

        MoleculeDefinitionTestCase('''A_ppi_B; x A--C
                                      A_ppi_C
                                      B_ppi_E; ! A--B
                                      B_p+_E; ! B--E''',
                                   {rfs.specification_from_string('A'): mol_def_from_string('A#ass/A_[Bassoc]:B_[Aassoc],ass/A_[Cassoc]:C_[Aassoc]'),
                                    rfs.specification_from_string('B'): mol_def_from_string('B#ass/B_[Eassoc]:E_[Bassoc],ass/B_[Aassoc]:A_[Bassoc]'),
                                    rfs.specification_from_string('C'): mol_def_from_string('C#ass/C_[Aassoc]:A_[Cassoc]'),
                                    rfs.specification_from_string('E'): mol_def_from_string('E#mod/E_[(Bsite)]:u~p,ass/E_[Bassoc]:B_[Eassoc]')}
                                   ),

        MoleculeDefinitionTestCase("""A_p+_B_[(x)]
                                      C_p+_B_[(x)]""",
                                   {rfs.specification_from_string('A'): mol_def_from_string('A#'),
                                    rfs.specification_from_string('B'): mol_def_from_string('B#mod/B_[(x)]:u~p'),
                                    rfs.specification_from_string('C'): mol_def_from_string('C#')}
                                   ),

        MoleculeDefinitionTestCase("""A_p+_B_[(x)]
                                      C_ub+_B_[(x)]""",
                                   {rfs.specification_from_string('A'): mol_def_from_string('A#'),
                                    rfs.specification_from_string('B'): mol_def_from_string('B#mod/B_[(x)]:u~p~ub'),
                                    rfs.specification_from_string('C'): mol_def_from_string('C#')}
                                   ),

        MoleculeDefinitionTestCase("""A_ppi_B_[x]
                                      C_ppi_B_[x]""",
                                   {rfs.specification_from_string('A'): mol_def_from_string('A#ass/A_[Bassoc]:B_[x]'),
                                    rfs.specification_from_string('B'): mol_def_from_string('B#ass/B_[x]:A_[Bassoc]~C_[Bassoc]'),
                                    rfs.specification_from_string('C'): mol_def_from_string('C#ass/C_[Bassoc]:B_[x]')}
                                   )
    ]