import pytest
from collections import namedtuple

import rxncon.input.quick.quick as qui

import rxncon.syntax.rxncon_from_string as rfs
import rxncon.core.contingency as con
import rxncon.core.effector as eff
import rxncon.core.reaction as rxn
import rxncon.core.state as sta
QuickTestCase = namedtuple('QuickTestCase', ["string", "reactions", "contingencies"])

expected_effector1 = eff.AndEffector(eff.AndEffector(eff.StateEffector(sta.state_from_string('A--C')),
                                                     eff.StateEffector(sta.state_from_string('A--D'))),
                                     eff.StateEffector(sta.state_from_string('B--C')))

expected_effector1.name = '<bool>'

expected_effector2 = eff.AndEffector(eff.StateEffector(sta.state_from_string('B--E')),
                                     eff.StateEffector(sta.state_from_string('B--F')))
expected_effector2.name = '<bool1>'

expected_effector3 = eff.OrEffector(eff.StateEffector(sta.state_from_string('B--D')), expected_effector2)
expected_effector3.name = '<bool>'


def test_quick_input(the_case_quick):
    for the_case in the_case_quick:
        is_quick_output_correct(the_case)


def is_quick_output_correct(the_case):
    quick_system = qui.Quick(the_case.string)

    assert all(reaction in the_case.reactions for reaction in quick_system.rxncon_system.reactions)
    assert all(reaction in quick_system.rxncon_system.reactions for reaction in the_case.reactions)
    assert all(contingency in the_case.contingencies for contingency in quick_system.rxncon_system.contingencies)
    assert all(contingency in quick_system.rxncon_system.contingencies for contingency in the_case.contingencies)


@pytest.fixture
def the_case_quick():
    return [
        QuickTestCase('A_ppi_B',
                      [rxn.reaction_from_string('A_ppi_B')],
                      []),

        QuickTestCase("""A_ppi_C
                      A_ppi_B; ! A--C""",
                      [rxn.reaction_from_string('A_ppi_B'), rxn.reaction_from_string('A_ppi_C')],
                      [con.Contingency(rxn.reaction_from_string('A_ppi_B'),con.ContingencyType.requirement,
                                       eff.StateEffector(sta.state_from_string('A--C')))]),

        QuickTestCase("""C_p+_A
                        C_ppi_A
                        A_ppi_B; ! A-{P}; ! A--C""",
                      [rxn.reaction_from_string('C_p+_A'), rxn.reaction_from_string('C_ppi_A'),
                       rxn.reaction_from_string('A_ppi_B')],
                      [con.Contingency(rxn.reaction_from_string('A_ppi_B'), con.ContingencyType.requirement,
                                       eff.StateEffector(sta.state_from_string('A--C'))),
                       con.Contingency(rxn.reaction_from_string('A_ppi_B'), con.ContingencyType.requirement,
                                       eff.StateEffector(sta.state_from_string('A-{p}')))]),

        QuickTestCase("""A_ppi_B; ! <bool>
                        <bool>; AND A--C
                        <bool>; AND A--D
                        <bool>; AND B--C""",
                      [rxn.reaction_from_string('A_ppi_B')],
                      [con.Contingency(rxn.reaction_from_string('A_ppi_B'), con.ContingencyType.requirement, expected_effector1)]),

        QuickTestCase("""A_ppi_B; ! A_[n]--[m]
                         B_p+_C; ! <bool>
                         <bool>; OR B--D
                         <bool>; OR <bool1>
                         <bool1>; AND B--E
                         <bool1>; AND B--F
                        """,
                      [rxn.reaction_from_string('A_ppi_B'), rxn.reaction_from_string('B_p+_C')],
                      [con.Contingency(rxn.reaction_from_string('A_ppi_B'), con.ContingencyType.requirement,
                                       eff.StateEffector(sta.state_from_string('A_[n]--[m]'))),
                       con.Contingency(rxn.reaction_from_string('B_p+_C'), con.ContingencyType.requirement,
                                       expected_effector3)]),

        QuickTestCase("""A_ppi_B; ! A-{P}; ! [Input]
                        C_p+_A
                        [Output]; ! A--B""",
                      [rxn.reaction_from_string('A_ppi_B'), rxn.reaction_from_string('C_p+_A')],
                      [con.Contingency(rxn.reaction_from_string('A_ppi_B'), con.ContingencyType.requirement,
                                       eff.StateEffector(sta.state_from_string('A-{p}'))),
                       con.Contingency(rxn.reaction_from_string('A_ppi_B'), con.ContingencyType.requirement,
                                       eff.StateEffector(sta.state_from_string('[Input]'))),
                       con.Contingency(rxn.reaction_from_string('[Output]'), con.ContingencyType.requirement,
                                       eff.StateEffector(sta.state_from_string('A--B')))])
    ]
