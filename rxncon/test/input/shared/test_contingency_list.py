import pytest
from collections import namedtuple

import rxncon.core.contingency as con
import rxncon.core.effector as eff
import rxncon.core.state as sta
import rxncon.input.shared.contingency_list as cli
import rxncon.syntax.rxncon_from_string as fst


ContingencyListTestCase = namedtuple('ContingencyListTestCase', ['entry', 'expected_contingency_name',
                                                                 'expected_operator', 'expected_state', 'expected_agent',
                                                                 'expected_agent_string'])


def test_contingency_list_entry_boolean_subject_state_agent(the_case_contingency_list_entry_boolean_subject_state_agent):
    for the_case in the_case_contingency_list_entry_boolean_subject_state_agent:
        is_boolean_entry_correct(the_case)


def test_contingency_list_entry_reaction_subject_state_agent(the_case_contingency_list_entry_reaction_subject_state_agent):
    for the_case in the_case_contingency_list_entry_reaction_subject_state_agent:
        is_reaction_entry_correct(the_case)


def is_boolean_entry_correct(the_case):
    assert the_case.entry.is_boolean_entry
    assert the_case.entry.subject == the_case.expected_contingency_name
    assert the_case.entry.predicate == the_case.expected_operator
    assert isinstance(the_case.entry.agent, the_case.expected_state)
    assert the_case.entry.agent == the_case.expected_agent
    assert str(the_case.entry.agent) == the_case.expected_agent_string


def is_reaction_entry_correct(the_case):
    assert the_case.entry.is_reaction_entry
    assert the_case.entry.subject == the_case.expected_contingency_name
    assert the_case.entry.predicate == the_case.expected_operator
    assert isinstance(the_case.entry.agent, the_case.expected_state)
    assert the_case.entry.agent == the_case.expected_agent
    assert str(the_case.entry.agent) == the_case.expected_agent_string


@pytest.fixture
def the_case_contingency_list_entry_boolean_subject_state_agent():
    return [
        ContingencyListTestCase(cli.contingency_list_entry_from_subject_predicate_agent_strings('<Ste11^{M/5}>', 'AND', 'Ste5_[MEKK]--Ste11'),
                                cli.BooleanContingencyName('<Ste11^{M/5}>'),
                                cli.BooleanOperator.op_and,
                                sta.InteractionState,
                                fst.state_from_string('Ste5_[MEKK]--Ste11'),
                                'Ste5_[MEKK]--Ste11'),

        ContingencyListTestCase(cli.contingency_list_entry_from_subject_predicate_agent_strings('<Ste11^{M/5}>', 'NOT', 'Ste5_[MEKK]--Ste11'),
                                cli.BooleanContingencyName('<Ste11^{M/5}>'),
                                cli.BooleanOperator.op_not,
                                sta.InteractionState,
                                fst.state_from_string('Ste5_[MEKK]--Ste11'),
                                'Ste5_[MEKK]--Ste11'),

        ContingencyListTestCase(cli.contingency_list_entry_from_subject_predicate_agent_strings('<Cdc24^{M}>', 'OR', '<Cdc24^{M/4}>'),
                                cli.BooleanContingencyName('<Cdc24^{M}>'),
                                cli.BooleanOperator.op_or,
                                cli.BooleanContingencyName,
                                cli.BooleanContingencyName('<Cdc24^{M/4}>'),
                                '<Cdc24^{M/4}>'),
    ]


@pytest.fixture
def the_case_contingency_list_entry_reaction_subject_state_agent():
    return [
        ContingencyListTestCase(cli.contingency_list_entry_from_subject_predicate_agent_strings('A_ppi_B', '!', 'A-{P}'),
                                fst.reaction_from_string('A_ppi_B'),
                                con.ContingencyType.requirement,
                                sta.CovalentModificationState,
                                fst.state_from_string('A-{P}'),
                                'A-{p}'),

        ContingencyListTestCase(cli.contingency_list_entry_from_subject_predicate_agent_strings('A_ppi_B', 'x', 'A-{P}'),
                                fst.reaction_from_string('A_ppi_B'),
                                con.ContingencyType.inhibition,
                                sta.CovalentModificationState,
                                fst.state_from_string('A-{P}'),
                                'A-{p}'),

        ContingencyListTestCase(cli.contingency_list_entry_from_subject_predicate_agent_strings('A_ppi_B', 'k+', 'A-{P}'),
                                fst.reaction_from_string('A_ppi_B'),
                                con.ContingencyType.positive,
                                sta.CovalentModificationState,
                                fst.state_from_string('A-{P}'),
                                'A-{p}'),

        ContingencyListTestCase(cli.contingency_list_entry_from_subject_predicate_agent_strings('A_ppi_B', 'k-', 'A-{P}'),
                                fst.reaction_from_string('A_ppi_B'),
                                con.ContingencyType.negative,
                                sta.CovalentModificationState,
                                fst.state_from_string('A-{P}'),
                                'A-{p}'),

        ContingencyListTestCase(cli.contingency_list_entry_from_subject_predicate_agent_strings('A_ppi_B', '0', 'A-{P}'),
                                fst.reaction_from_string('A_ppi_B'),
                                con.ContingencyType.no_effect,
                                sta.CovalentModificationState,
                                fst.state_from_string('A-{P}'),
                                'A-{p}'),

        ContingencyListTestCase(cli.contingency_list_entry_from_subject_predicate_agent_strings('A_ppi_B', '?', 'A-{P}'),
                                fst.reaction_from_string('A_ppi_B'),
                                con.ContingencyType.unknown,
                                sta.CovalentModificationState,
                                fst.state_from_string('A-{P}'),
                                'A-{p}')
    ]


LookupTableTestCase = namedtuple('LookupTableTestCase', ['lookup_table', 'look_up', 'expected_effector', 'expected_lookup_table_len'])


def test_lookup_table_from_contingency_list_entries(the_case_lookup_table):
    for the_case in the_case_lookup_table:
        is_lookup_table_correct(the_case)


def is_lookup_table_correct(the_case):
    if the_case.look_up:
        assert the_case.lookup_table[the_case.look_up] == the_case.expected_effector
    else:
        assert the_case.lookup_table == the_case.expected_effector

    assert len(the_case.lookup_table) == the_case.expected_lookup_table_len


# {BooleanContingencyName -> Effector} lookup table from [ContingencyListEntry]
@pytest.fixture
def the_case_lookup_table():
    return [
        LookupTableTestCase(cli._create_boolean_contingency_lookup_table([]),
                            '',
                            {},
                            0),

        LookupTableTestCase(cli._create_boolean_contingency_lookup_table([cli.contingency_list_entry_from_subject_predicate_agent_strings('<X>', 'AND', 'A-{P}'),
                                                                          cli.contingency_list_entry_from_subject_predicate_agent_strings('<X>', 'AND', 'A--B')]),
                            cli.BooleanContingencyName('<X>'),
                            eff.AndEffector(eff.StateEffector(fst.state_from_string('A-{P}')),
                                            eff.StateEffector(fst.state_from_string('A--B'))),
                            1),

        LookupTableTestCase(cli._create_boolean_contingency_lookup_table([cli.contingency_list_entry_from_subject_predicate_agent_strings('<X>', 'AND', '<Y>'),
                                                                          cli.contingency_list_entry_from_subject_predicate_agent_strings('<X>', 'AND', '<Z>')]),
                            cli.BooleanContingencyName('<X>'),
                            eff.AndEffector(cli._BooleanContingencyEffector(cli.BooleanContingencyName('<Y>')),
                                            cli._BooleanContingencyEffector(cli.BooleanContingencyName('<Z>'))),
                            1)
    ]


# [Contingency] from [ContingencyListEntry]
ContingencyTestCase = namedtuple('ContingencyTestCase', ['contingency', 'expected_contingency'])


def test_contingencies_from_contingency_list_entries(the_case_contingencies_from_contingency_list_entries):
    for the_case in the_case_contingencies_from_contingency_list_entries:
        assert the_case.contingency == the_case.expected_contingency


@pytest.fixture
def the_case_contingencies_from_contingency_list_entries(boolean_contingencies):
    return [
        ContingencyTestCase(cli.contingencies_from_contingency_list_entries([cli.contingency_list_entry_from_subject_predicate_agent_strings('A_ppi_B', '!', 'A-{P}')]),
                            [con.Contingency(fst.reaction_from_string('A_ppi_B'),
                                             con.ContingencyType.requirement,
                                             eff.StateEffector(fst.state_from_string('A-{P}')))]
                            ),

        ContingencyTestCase(cli.contingencies_from_contingency_list_entries([cli.contingency_list_entry_from_subject_predicate_agent_strings('A_ppi_C', 'x', '<X>'),
                                                                             cli.contingency_list_entry_from_subject_predicate_agent_strings('<X>', 'AND', 'A-{P}'),
                                                                             cli.contingency_list_entry_from_subject_predicate_agent_strings('<X>', 'AND', 'A--B')]),
                            [con.Contingency(fst.reaction_from_string('A_ppi_C'),
                                             con.ContingencyType.inhibition,
                                             boolean_contingencies['expected_X'])]
                            ),

        ContingencyTestCase(cli.contingencies_from_contingency_list_entries([cli.contingency_list_entry_from_subject_predicate_agent_strings('A_ppi_C', 'x', '<Z>'),
                                                                             cli.contingency_list_entry_from_subject_predicate_agent_strings('<Z>', 'AND', 'A-{P}'),
                                                                             cli.contingency_list_entry_from_subject_predicate_agent_strings('<Z>', 'AND', '<Y>'),
                                                                             cli.contingency_list_entry_from_subject_predicate_agent_strings('<Y>', 'OR', 'A--B'),
                                                                             cli.contingency_list_entry_from_subject_predicate_agent_strings('<Y>', 'OR', 'A--D')]),
                            [con.Contingency(fst.reaction_from_string('A_ppi_C'),
                                             con.ContingencyType.inhibition,
                                             boolean_contingencies['expected_Z'])]
                            )

    ]


@pytest.fixture
def boolean_contingencies():
    bool_dict = {}
    bool_dict['expected_X'] = eff.AndEffector(eff.StateEffector(fst.state_from_string('A-{P}')),
                                                               eff.StateEffector(fst.state_from_string('A--B')))
    bool_dict['expected_X'].name = '<X>'

    bool_dict['expected_Y'] = eff.OrEffector(eff.StateEffector(fst.state_from_string('A--B')),
                                             eff.StateEffector(fst.state_from_string('A--D')))
    bool_dict['expected_Y'].name = '<Y>'

    bool_dict['expected_Z'] = eff.AndEffector(eff.StateEffector(fst.state_from_string('A-{P}')), bool_dict['expected_Y'])
    bool_dict['expected_Z'].name = '<Z>'

    return bool_dict
