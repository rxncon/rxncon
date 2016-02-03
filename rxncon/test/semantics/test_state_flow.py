import rxncon.syntax.rxncon_from_string as rfs
import rxncon.core.effector as eff
import rxncon.core.contingency as con
import rxncon.semantics.state_flow as flo
import rxncon.venntastic.sets as venn


### SET STUFF ###
def test_set_from_contingencies():
    reaction = rfs.reaction_from_string('A_ppi_B')
    state1 = rfs.state_from_string('A_[x]-{p}')
    state2 = rfs.state_from_string('B_[y]-{p}')
    state3 = rfs.state_from_string('A_[z]-{p}')
    state4 = rfs.state_from_string('B_[w]-{p}')

    boundstate = rfs.state_from_string('A--B')

    contingency = con.Contingency(reaction,
                                  con.ContingencyType.requirement,
                                  eff.AndEffector(eff.OrEffector(eff.StateEffector(state1),
                                                                 eff.StateEffector(state2)),
                                                  eff.OrEffector(eff.StateEffector(state3),
                                                                 eff.StateEffector(state4))))

    source_contingency = con.Contingency(reaction,
                                         con.ContingencyType.inhibition,
                                         eff.StateEffector(boundstate))

    state_flows = flo.boolean_state_flows([contingency], [source_contingency])

    expected_pairs = [
        (
            [[venn.PropertySet(state1), venn.PropertySet(state3), venn.Complement(venn.PropertySet(boundstate))]],
            [[venn.PropertySet(state1), venn.PropertySet(state3), venn.PropertySet(boundstate)]]
        ),
        (
            [[venn.PropertySet(state1), venn.PropertySet(state4), venn.Complement(venn.PropertySet(boundstate))]],
            [[venn.PropertySet(state1), venn.PropertySet(state4), venn.PropertySet(boundstate)]]
        ),
        (
            [[venn.PropertySet(state2), venn.PropertySet(state3), venn.Complement(venn.PropertySet(boundstate))]],
            [[venn.PropertySet(state2), venn.PropertySet(state3), venn.PropertySet(boundstate)]]
        ),
        (
            [[venn.PropertySet(state2), venn.PropertySet(state4), venn.Complement(venn.PropertySet(boundstate))]],
            [[venn.PropertySet(state2), venn.PropertySet(state4), venn.PropertySet(boundstate)]]
        )
    ]

    assert len(state_flows) == 4

    for flow in state_flows:
        if flow.source.to_nested_list_form() == expected_pairs[0][0]:
            assert flow.target.to_nested_list_form() == expected_pairs[0][1]

        elif flow.source.to_nested_list_form() == expected_pairs[1][0]:
            assert flow.target.to_nested_list_form() == expected_pairs[1][1]

        elif flow.source.to_nested_list_form() == expected_pairs[2][0]:
            assert flow.target.to_nested_list_form() == expected_pairs[2][1]

        elif flow.source.to_nested_list_form() == expected_pairs[3][0]:
            assert flow.target.to_nested_list_form() == expected_pairs[3][1]

        else:
            raise AssertionError



