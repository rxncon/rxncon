import pytest
from collections import namedtuple
from typing import List

from rxncon.venntastic.sets import PropertySet, Intersection, Union, Complement, nested_expression_from_list_and_binary_op
from rxncon.syntax.rxncon_from_string import state_from_string, specification_from_string, reaction_from_string
from rxncon.semantics.molecule_instance import MoleculeInstance
from rxncon.simulation.rule_based.molecule_from_string import mol_def_from_string, mol_instance_from_string
from rxncon.simulation.rule_based.molecule_from_rxncon import mol_instance_set_from_state_set, mol_instance_set_pair_from_reaction


MoleculeInstancesFromStateSetTestCase = namedtuple('MoleculeInstancesFromStateSetTestCase',
                                                   ['mol_defs', 'state_set', 'expected_mol_instances'])

MoleculeInstancesPairFromReactionTestCase = namedtuple('MoleculeInstancesPairFromReactionTestCase',
                                                       ['mol_defs', 'reaction', 'expected_mol_instances_pair'])


def test_molecule_instances_from_state_set(state_set_test_cases):
    for test_case in state_set_test_cases:
        assert is_state_set_test_case_correct(test_case)


def test_molecule_instance_pair_from_reaction(reaction_test_cases):
    for test_case in reaction_test_cases:
        assert is_reaction_test_case_correct(test_case)


@pytest.fixture
def molecule_definitions():
    return {
        specification_from_string('A'): mol_def_from_string('A#ass/A_[Bassoc]:B_[Aassoc],mod/A_[(r)]:u~p~ub'),
        specification_from_string('B'): mol_def_from_string('B#ass/B_[Aassoc]:A_[Bassoc]')
    }


@pytest.fixture
def state_set_test_cases(molecule_definitions):
    return [
        MoleculeInstancesFromStateSetTestCase(
            molecule_definitions,
            PropertySet(state_from_string('A--B')),
            [
                {'A': 'A#ass/A_[Bassoc]:B_[Aassoc]', 'B': 'B#ass/B_[Aassoc]:A_[Bassoc]'}
            ]
        ),
        MoleculeInstancesFromStateSetTestCase(
            molecule_definitions,
            Complement(PropertySet(state_from_string('A-{p}'))),
            [
                {'A': 'A#mod/A_[(r)]:u'},
                {'A': 'A#mod/A_[(r)]:ub'}
            ]
        )
    ]


@pytest.fixture
def reaction_test_cases(molecule_definitions):
    return [
        MoleculeInstancesPairFromReactionTestCase(
            molecule_definitions,
            reaction_from_string('A_ppi_B'),
            (
                {'A': 'A#ass/A_[Bassoc]:', 'B': 'B#ass/B_[Aassoc]:'},
                {'A': 'A#ass/A_[Bassoc]:B_[Aassoc]', 'B': 'B#ass/B_[Aassoc]:A_[Bassoc]'}
            )
        ),
        MoleculeInstancesPairFromReactionTestCase(
            molecule_definitions,
            reaction_from_string('B_p+_A_[(r)]'),
            (
                {'B': 'B#', 'A': 'A#mod/A_[(r)]:u'},
                {'B': 'B#', 'A': 'A#mod/A_[(r)]:p'}
            )
        )
    ]


def is_state_set_test_case_correct(test_case: MoleculeInstancesFromStateSetTestCase) -> bool:
    actual_wrapped_mol_instance_lists = mol_instance_set_from_state_set(test_case.mol_defs, test_case.state_set).to_nested_list_form()
    actual_mol_instance_lists = []

    for wrapped_list in actual_wrapped_mol_instance_lists:
        assert all(isinstance(x, PropertySet) for x in wrapped_list)
        actual_mol_instance_lists.append([x.value for x in wrapped_list])

    expected_mol_instance_lists = [[mol_instance_from_string(test_case.mol_defs[specification_from_string(k)], v)
                                    for k, v in expected_dict.items()] for expected_dict in test_case.expected_mol_instances]

    for actual_list in actual_mol_instance_lists:
        matches = 0
        for expected_list in expected_mol_instance_lists:
            if are_mol_instance_lists_equivalent(actual_list, expected_list):
                matches += 1

        if matches != 1:
            return False

    return True

def is_reaction_test_case_correct(test_case: MoleculeInstancesPairFromReactionTestCase) -> bool:
    actual_lhs_set, actual_rhs_set = mol_instance_set_pair_from_reaction(test_case.mol_defs, test_case.reaction)

    expected_lhs_set = nested_expression_from_list_and_binary_op(
        [PropertySet(mol_instance_from_string(test_case.mol_defs[specification_from_string(k)], v))
         for k, v in test_case.expected_mol_instances_pair[0].items()],
        Intersection
    )

    expected_rhs_set = nested_expression_from_list_and_binary_op(
        [PropertySet(mol_instance_from_string(test_case.mol_defs[specification_from_string(k)], v))
         for k, v in test_case.expected_mol_instances_pair[1].items()],
        Intersection
    )

    return actual_lhs_set.is_equivalent_to(expected_lhs_set) and actual_rhs_set.is_equivalent_to(expected_rhs_set)

def are_mol_instance_lists_equivalent(first_list: List[MoleculeInstance], second_list: List[MoleculeInstance]) -> bool:
    return all(x in second_list for x in first_list) and all(x in first_list for x in second_list) and \
        len(first_list) == len(second_list)


