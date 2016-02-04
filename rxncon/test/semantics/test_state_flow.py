import rxncon.syntax.rxncon_from_string as rfs
import rxncon.core.effector as eff
import rxncon.core.contingency as con
import rxncon.semantics.state_flow as flo
import rxncon.venntastic.sets as venn


def test_boolean_state_flows():
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

    state_flows = flo.boolean_state_flows(reaction, [contingency], [source_contingency])

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


def test_quantified_state_flows_single_quantitative_contingency():
    reaction = rfs.reaction_from_string('A_ppi_B')
    a_b = rfs.state_from_string('A--B')
    a_phos = rfs.state_from_string('A-{p}')
    b_phos = rfs.state_from_string('B-{p}')

    reqd_cont = con.Contingency(reaction, con.ContingencyType.requirement, eff.StateEffector(a_phos))
    src_cont  = con.Contingency(reaction, con.ContingencyType.requirement, eff.NotEffector(eff.StateEffector(a_b)))
    quan_cont = con.Contingency(reaction, con.ContingencyType.positive, eff.StateEffector(b_phos))

    bool_state_flows = flo.boolean_state_flows(reaction, [reqd_cont], [src_cont])
    assert len(bool_state_flows) == 1

    quant_state_flows = flo.quantified_state_flows(bool_state_flows[0], [quan_cont])

    expected_pairs = [
        (
            [[venn.PropertySet(a_phos), venn.Complement(venn.PropertySet(a_b)), venn.PropertySet(b_phos)]],
            [[venn.PropertySet(a_phos), venn.PropertySet(a_b), venn.PropertySet(b_phos)]],
        ),
        (
            [[venn.PropertySet(a_phos), venn.Complement(venn.PropertySet(a_b)), venn.Complement(venn.PropertySet(b_phos))]],
            [[venn.PropertySet(a_phos), venn.PropertySet(a_b), venn.Complement(venn.PropertySet(b_phos))]],
        )
    ]

    for flow in quant_state_flows:
        if flow.source.to_nested_list_form() == expected_pairs[0][0]:
            assert flow.target.to_nested_list_form() == expected_pairs[0][1]

        elif flow.source.to_nested_list_form() == expected_pairs[1][0]:
            assert flow.target.to_nested_list_form() == expected_pairs[1][1]

        else:
            raise AssertionError


def test_disjunctified_state_flows():
    reaction = rfs.reaction_from_string('A_ppi_B')
    a_phos = rfs.state_from_string('A-{p}')
    b_phos = rfs.state_from_string('B-{p}')
    a_b = rfs.state_from_string('A--B')

    contingency = con.Contingency(reaction,
                                  con.ContingencyType.requirement,
                                  eff.OrEffector(eff.StateEffector(a_phos),
                                                 eff.StateEffector(b_phos)))

    source_contingency = con.Contingency(reaction, con.ContingencyType.inhibition, eff.StateEffector(a_b))

    state_flows = flo.boolean_state_flows(reaction, [contingency], [source_contingency])
    disjunct_state_flows = flo.disjunctified_state_flows(state_flows)

    assert len(disjunct_state_flows) == 2

    assert disjunct_state_flows[0].source.is_equivalent_to(venn.Intersection(venn.PropertySet(a_phos),
                                                                             venn.Complement(venn.PropertySet(a_b))))
    assert disjunct_state_flows[0].target.is_equivalent_to(venn.Intersection(venn.PropertySet(a_phos),
                                                                             venn.PropertySet(a_b)))

    assert disjunct_state_flows[1].source.is_equivalent_to(venn.Intersection(venn.Intersection(venn.PropertySet(b_phos),
                                                                                               venn.Complement(venn.PropertySet(a_b))),
                                                                             venn.Complement(venn.PropertySet(a_phos))))
    assert disjunct_state_flows[1].target.is_equivalent_to(venn.Intersection(venn.Intersection(venn.PropertySet(b_phos),
                                                                                               venn.PropertySet(a_b)),
                                                                             venn.Complement(venn.PropertySet(a_phos))))
