import rxncon.input.quick.quick as qui
import rxncon.syntax.rxncon_from_string as rfs
import rxncon.core.contingency as con
import rxncon.core.effector as eff


def test_quick_single_reaction():
    quick = qui.Quick("A_ppi_B")

    assert isinstance(quick, qui.Quick)

    assert len(quick.rxncon_system.reactions) == 1
    assert len(quick.rxncon_system.contingencies) == 0

    assert quick.rxncon_system.reactions[0] == rfs.reaction_from_string('A_ppi_B')


def test_quick_single_reaction_simple_contingency():
    quick = qui.Quick("""
                      A_ppi_C
                      A_ppi_B; ! A--C""")

    assert isinstance(quick, qui.Quick)

    assert len(quick.rxncon_system.reactions) == 2
    assert quick.rxncon_system.reactions[0] != quick.rxncon_system.reactions[1]
    assert rfs.reaction_from_string('A_ppi_B') in quick.rxncon_system.reactions
    assert rfs.reaction_from_string('A_ppi_C') in quick.rxncon_system.reactions

    assert len(quick.rxncon_system.contingencies) == 1
    assert quick.rxncon_system.contingencies[0] == con.Contingency(
            rfs.reaction_from_string('A_ppi_B'),
            con.ContingencyType.requirement,
            eff.StateEffector(rfs.state_from_string('A--C'))
    )


def test_quick_single_reaction_two_contingencies():
    quick = qui.Quick("""
                        C_p+_A
                        C_ppi_A
                        A_ppi_B; ! A-{P}; ! A--C""")

    assert isinstance(quick, qui.Quick)
    assert len(quick.rxncon_system.reactions) == 3
    assert rfs.reaction_from_string('C_p+_A') in quick.rxncon_system.reactions
    assert rfs.reaction_from_string('C_ppi_A') in quick.rxncon_system.reactions
    assert rfs.reaction_from_string('A_ppi_B') in quick.rxncon_system.reactions

    assert len(quick.rxncon_system.contingencies) == 2
    expected_contingencies = [con.Contingency(rfs.reaction_from_string('A_ppi_B'), con.ContingencyType.requirement,
                                              eff.StateEffector(rfs.state_from_string('A--C'))),
                              con.Contingency(rfs.reaction_from_string('A_ppi_B'), con.ContingencyType.requirement,
                                              eff.StateEffector(rfs.state_from_string('A-{p}')))
                              ]
    assert quick.rxncon_system.contingencies[0] != quick.rxncon_system.contingencies[1]
    assert quick.rxncon_system.contingencies[0] in expected_contingencies
    assert quick.rxncon_system.contingencies[1] in expected_contingencies


def test_quick_single_reaction_simple_boolean_contingency():
    quick = qui.Quick("""A_ppi_B; ! <bool>
                        <bool>; AND A--C
                        <bool>; AND A--D
                        <bool>; AND B--C""")

    assert isinstance(quick, qui.Quick)

    assert len(quick.rxncon_system.reactions) == 1
    assert len(quick.rxncon_system.contingencies) == 1

    expected_effector = eff.AndEffector(eff.AndEffector(eff.StateEffector(rfs.state_from_string('A--C')),
                                                        eff.StateEffector(rfs.state_from_string('A--D'))),
                                        eff.StateEffector(rfs.state_from_string('B--C')))

    expected_effector.name = '<bool>'

    assert quick.rxncon_system.contingencies[0] == con.Contingency(
        rfs.reaction_from_string('A_ppi_B'),
        con.ContingencyType.requirement,
        expected_effector
    )


def test_quick_multiple_reactions_contingencies():
    quick = qui.Quick("""A_ppi_B; ! A_[n]--[m]
                         B_p+_C; ! <bool>
                         <bool>; OR B--D
                         <bool>; OR <bool1>
                         <bool1>; AND B--E
                         <bool1>; AND B--F
                        """)

    assert isinstance(quick, qui.Quick)

    assert len(quick.rxncon_system.reactions) == 2
    assert len(quick.rxncon_system.contingencies) == 2

    bool1_effector = eff.AndEffector(eff.StateEffector(rfs.state_from_string('B--E')),
                                     eff.StateEffector(rfs.state_from_string('B--F')))
    bool1_effector.name = '<bool1>'

    bool_effector = eff.OrEffector(eff.StateEffector(rfs.state_from_string('B--D')), bool1_effector)
    bool_effector.name = '<bool>'

    assert con.Contingency(
        rfs.reaction_from_string('A_ppi_B'),
        con.ContingencyType.requirement,
        eff.StateEffector(rfs.state_from_string('A_[n]--[m]'))
    ) in quick.rxncon_system.contingencies

    assert con.Contingency(
        rfs.reaction_from_string('B_p+_C'),
        con.ContingencyType.requirement,
        bool_effector
    ) in quick.rxncon_system.contingencies


