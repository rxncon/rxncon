import rxncon.core.contingency as con
import rxncon.core.effector as eff
import rxncon.core.state as sta
import rxncon.input.shared.contingency_list as cli
import rxncon.input.shared.from_string as fst


# ContingencyListEntry from strings
def test_contingency_list_entry_boolean_subject_state_agent():
    entry = cli.contingency_list_entry_from_subject_predicate_agent_strings('<Ste11^{M/5}>', 'AND', 'Ste5_[MEKK]--Ste11')

    assert entry.is_boolean_entry
    assert entry.subject == cli.BooleanContingencyName('<Ste11^{M/5}>')
    assert entry.predicate == cli.BooleanOperator.op_and
    assert isinstance(entry.agent, sta.InteractionState)
    assert entry.agent.full_name == 'Ste5_[MEKK]--Ste11'


def test_contingency_list_entry_boolean_subject_boolean_agent():
    entry = cli.contingency_list_entry_from_subject_predicate_agent_strings('<Cdc24^{M}>', 'OR', '<Cdc24^{M/4}>')

    assert entry.is_boolean_entry
    assert entry.subject == cli.BooleanContingencyName('<Cdc24^{M}>')
    assert entry.predicate == cli.BooleanOperator.op_or
    assert entry.agent == cli.BooleanContingencyName('<Cdc24^{M/4}>')


def test_contingency_list_entry_reaction_subject_state_agent():
    entry = cli.contingency_list_entry_from_subject_predicate_agent_strings('A_ppi_B', '!', 'A-{P}')

    assert entry.is_reaction_entry
    assert entry.subject == fst.reaction_from_string('A_ppi_B')
    assert entry.predicate == con.ContingencyType.requirement
    assert entry.agent == fst.state_from_string('A-{P}')


# {BooleanContingencyName -> Effector} lookup table from [ContingencyListEntry]
def test_lookup_table_from_contingency_list_entries_empty():
    lookup_table = cli._create_boolean_contingency_lookup_table([])

    assert lookup_table == {}


def test_lookup_table_from_contingency_list_entries_boolean_subjects_state_agents():
    entries = [
        cli.contingency_list_entry_from_subject_predicate_agent_strings('<X>', 'AND', 'A-{P}'),
        cli.contingency_list_entry_from_subject_predicate_agent_strings('<X>', 'AND', 'A--B')
    ]
    lookup_table = cli._create_boolean_contingency_lookup_table(entries)

    expected_effector = eff.AndEffector(eff.StateEffector(fst.state_from_string('A-{P}')), eff.StateEffector(fst.state_from_string('A--B')))

    assert lookup_table[cli.BooleanContingencyName('<X>')] == expected_effector
    assert len(lookup_table) == 1


def test_lookup_table_from_contingency_list_entries_boolean_subjects_boolean_agents():
    entries = [
        cli.contingency_list_entry_from_subject_predicate_agent_strings('<X>', 'AND', '<Y>'),
        cli.contingency_list_entry_from_subject_predicate_agent_strings('<X>', 'AND', '<Z>')
    ]
    lookup_table = cli._create_boolean_contingency_lookup_table(entries)

    expected_effector = eff.AndEffector(cli._BooleanContingencyEffector(cli.BooleanContingencyName('<Y>')),
                                        cli._BooleanContingencyEffector(cli.BooleanContingencyName('<Z>')))

    assert lookup_table[cli.BooleanContingencyName('<X>')] == expected_effector
    assert len(lookup_table) == 1


# [Contingency] from [ContingencyListEntry]
def test_contingencies_from_contingency_list_entries_single():
    entry = cli.contingency_list_entry_from_subject_predicate_agent_strings('A_ppi_B', '!', 'A-{P}')
    contingencies = cli.contingencies_from_contingency_list_entries([entry])

    assert contingencies == [con.Contingency(fst.reaction_from_string('A_ppi_B'),
                                             con.ContingencyType.requirement,
                                             eff.StateEffector(fst.state_from_string('A-{P}')))]


def test_contingencies_from_contingency_list_entries_boolean_flat():
    entries = [
        cli.contingency_list_entry_from_subject_predicate_agent_strings('A_ppi_C', 'x', '<X>'),
        cli.contingency_list_entry_from_subject_predicate_agent_strings('<X>', 'AND', 'A-{P}'),
        cli.contingency_list_entry_from_subject_predicate_agent_strings('<X>', 'AND', 'A--B')
    ]
    contingencies = cli.contingencies_from_contingency_list_entries(entries)

    assert contingencies == [con.Contingency(fst.reaction_from_string('A_ppi_C'),
                                             con.ContingencyType.inhibition,
                                             eff.AndEffector(eff.StateEffector(fst.state_from_string('A-{P}')),
                                                             eff.StateEffector(fst.state_from_string('A--B'))))]


def test_contingencies_from_contingency_list_entries_boolean_nested():
    entries = [
        cli.contingency_list_entry_from_subject_predicate_agent_strings('A_ppi_C', 'x', '<X>'),
        cli.contingency_list_entry_from_subject_predicate_agent_strings('<X>', 'AND', 'A-{P}'),
        cli.contingency_list_entry_from_subject_predicate_agent_strings('<X>', 'AND', '<Y>'),
        cli.contingency_list_entry_from_subject_predicate_agent_strings('<Y>', 'OR', 'A--B'),
        cli.contingency_list_entry_from_subject_predicate_agent_strings('<Y>', 'OR', 'A--D')
    ]
    contingencies = cli.contingencies_from_contingency_list_entries(entries)

    assert contingencies == [con.Contingency(fst.reaction_from_string('A_ppi_C'),
                                             con.ContingencyType.inhibition,
                                             eff.AndEffector(eff.StateEffector(fst.state_from_string('A-{P}')),
                                                             eff.OrEffector(eff.StateEffector(fst.state_from_string('A--B')),
                                                                            eff.StateEffector(fst.state_from_string('A--D')))))]

