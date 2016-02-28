import pytest

import rxncon.core.rxncon_system as rxs
import rxncon.simulation
import rxncon.syntax.rxncon_from_string as rfs
import rxncon.core.contingency as con
import rxncon.core.effector as eff
import rxncon.simulation.rule_based.rbm_from_rxncon as rfr
from rxncon.core import contingency as con, effector, effector, effector, effector, rxncon_system
from rxncon.input.quick import quick as qui
from rxncon.syntax import rxncon_from_string as rfs
from rxncon.venntastic import sets as venn


def test_single_rule():
    c_phos_a = rfs.reaction_from_string('C_p+_A_[x]')
    a_phos = rfs.state_from_string('A_[x]-{p}')
    a_phos_b = rfs.reaction_from_string('A_p+_B')

    cont = con.Contingency(a_phos_b, con.ContingencyType.requirement, eff.StateEffector(a_phos))

    rxncon = rxs.RxnConSystem([c_phos_a, a_phos_b], [cont])
    rule_supervisor = rfr.RuleBasedModelSupervisor(rxncon)

    for rule in rule_supervisor.rules:
        print()
        print(rule)



# def test_single_ppi_rule():
#     a_ppi_b = rfs.reaction_from_string('A_ppi_B')
#     a_ppi_c = rfs.reaction_from_string('A_ppi_C')
#     #c_ppi_f = rfs.reaction_from_string('C_ppi_F')
#     b_ppi_e = rfs.reaction_from_string('B_ppi_E')
#     a_dash_dash_c = rfs.state_from_string('A--C')
#     b_dash_dash_e = rfs.state_from_string('B--E')
#     #c_dash_dash_f = rfs.state_from_string('C--F')
#     # A--B
#     cont_A = con.Contingency(a_ppi_b, con.ContingencyType.requirement, eff.StateEffector(a_dash_dash_c))
#     cont_B = con.Contingency(a_ppi_b, con.ContingencyType.requirement, eff.StateEffector(b_dash_dash_e))
#     #cont_C = con.Contingency(a_ppi_c, con.ContingencyType.requirement, eff.StateEffector(c_dash_dash_f))
#     rxncon = rxs.RxnConSystem([a_ppi_b, a_ppi_c, b_ppi_e], [cont_A, cont_B])
#
#     rules = rfr.rules_from_rxncon(rxncon)

@pytest.fixture
def simple_system():
    phosphorylation_reaction_1 = rfs.reaction_from_string('A_p+_X')
    phosphorylation_reaction_2 = rfs.reaction_from_string('B_p+_X')
    phosphorylation_reaction_3 = rfs.reaction_from_string('C_p+_X')

    binding_reaction = rfs.reaction_from_string('X_ppi_Y')
    phosphorylated_state = rfs.state_from_string('X-{p}')

    binding_contingency = con.Contingency(binding_reaction,
                                          con.ContingencyType.requirement,
                                          eff.StateEffector(phosphorylated_state))

    reactions = [phosphorylation_reaction_1, phosphorylation_reaction_2, phosphorylation_reaction_3, binding_reaction]
    contingencies = [binding_contingency]

    return rxs.RxnConSystem(reactions, contingencies)


def test_set_of_states_from_single_state_effector():
    a_ppi_c = rfs.reaction_from_string("A_ppi_C")
    a_dash_d = rfs.state_from_string('A--D')

    cont = con.Contingency(a_ppi_c, con.ContingencyType.requirement, eff.StateEffector(a_dash_d))

    set_of_state_effector_a_dash_d = rxncon.simulation.rule_based.rbm_from_rxncon.state_set_from_effector(cont.effector)
    assert set_of_state_effector_a_dash_d.is_equivalent_to(venn.PropertySet(a_dash_d))


def test_set_of_states_from_nested_AND_effectors():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; AND A--C
                        <comp>; AND C--E
                        <comp>; AND B--F""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_f = rfs.state_from_string("B--F")

    rxncon = quick.rxncon_system

    set_of_state_AND_effector = rxncon.simulation.rule_based.rbm_from_rxncon.set_of_states_from_effector(rxncon.contingencies[0].effector)

    assert set_of_state_AND_effector.is_equivalent_to(venn.Intersection(venn.Intersection(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.PropertySet(expected_b_dash_f)))


def test_set_of_states_from_nested_OR_effectors():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; OR A--C
                        <comp>; OR C--E
                        <comp>; OR B--F""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_f = rfs.state_from_string("B--F")

    rxncon = quick.rxncon_system

    set_of_state_AND_effector = rxncon.simulation.rule_based.rbm_from_rxncon.set_of_states_from_effector(rxncon.contingencies[0].effector)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Union(venn.Union(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.PropertySet(expected_b_dash_f)))


def test_set_of_states_from_nested_AND_OR_effector():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; AND <c1>
                        <comp>; AND <c2>
                        <c1>; OR A--C
                        <c1>; OR C--E
                        <c2>; OR B--F
                        <c2>; OR B--D""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_f = rfs.state_from_string("B--F")
    expected_b_dash_d = rfs.state_from_string("B--D")

    rxncon = quick.rxncon_system

    set_of_state_AND_effector = rxncon.simulation.rule_based.rbm_from_rxncon.set_of_states_from_effector(rxncon.contingencies[0].effector)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Intersection(venn.Union(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.Union(venn.PropertySet(expected_b_dash_f),venn.PropertySet(expected_b_dash_d))))


def test_set_of_states_from_nested_OR_AND_effector():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; OR <c1>
                        <comp>; OR <c2>
                        <c1>; AND A--C
                        <c1>; AND C--E
                        <c2>; AND B--F
                        <c2>; AND B--D""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_f = rfs.state_from_string("B--F")
    expected_b_dash_d = rfs.state_from_string("B--D")


    rxncon = quick.rxncon_system

    set_of_state_AND_effector = rxncon.simulation.rule_based.rbm_from_rxncon.set_of_states_from_effector(rxncon.contingencies[0].effector)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Union(venn.Intersection(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.Intersection(venn.PropertySet(expected_b_dash_f),venn.PropertySet(expected_b_dash_d))))


def test_set_of_states_from_complex_effector():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; OR <c1>
                        <comp>; OR <c2>
                        <c1>; AND A--C
                        <c1>; AND C--E
                        <c2>; NOT <c3>
                        <c3>; AND B--F
                        <c3>; AND B--D""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_f = rfs.state_from_string("B--F")
    expected_b_dash_d = rfs.state_from_string("B--D")

    rxncon = quick.rxncon_system

    set_of_state_AND_effector = rxncon.simulation.rule_based.rbm_from_rxncon.set_of_states_from_effector(rxncon.contingencies[0].effector)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Union(venn.Intersection(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.Complement(venn.Intersection(venn.PropertySet(expected_b_dash_d),
                                                                                                          venn.PropertySet(expected_b_dash_f)))))


def test_set_of_states_from_strict_contingencies():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_dash_b = rfs.state_from_string('A--B')

    a_ppi_c = rfs.reaction_from_string('A_ppi_C')
    a_dash_c = rfs.state_from_string('A--C')

    b_ppi_e = rfs.reaction_from_string('B_ppi_E')

    b_pplus_e = rfs.reaction_from_string('B_p+_E')
    e_pplus = rfs.state_from_string("E-{P}")

    cont_b_dash_e = con.Contingency(b_ppi_e, con.ContingencyType.requirement, eff.StateEffector(a_dash_b))  # B_ppi_E; ! A--B
    cont_e_pplus = con.Contingency(b_ppi_e, con.ContingencyType.requirement, eff.StateEffector(e_pplus))  # B_ppi_E; ! E-{P}
    cont_a_ppi_b = con.Contingency(a_ppi_b, con.ContingencyType.inhibition, eff.StateEffector(a_dash_c))  # A_ppi_B; x A--C

    rxncon = rxs.RxnConSystem([b_ppi_e, a_ppi_b, a_ppi_c, b_pplus_e], [cont_e_pplus, cont_b_dash_e, cont_a_ppi_b])

    strict_contingencies_state_set_b_ppi_e = rxncon.simulation.rule_based.rbm_from_rxncon.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(b_ppi_e))
    strict_contingencies_state_set_a_ppi_b = rxncon.simulation.rule_based.rbm_from_rxncon.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(a_ppi_b))
    strict_contingencies_state_set_a_ppi_c = rxncon.simulation.rule_based.rbm_from_rxncon.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(a_ppi_c))
    strict_contingencies_state_set_b_pplus_e = rxncon.simulation.rule_based.rbm_from_rxncon.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(b_pplus_e))

    expected_b_ppi_e_strict_cont = venn.Intersection(venn.PropertySet(e_pplus), venn.PropertySet(a_dash_b))
    assert strict_contingencies_state_set_b_ppi_e.is_equivalent_to(expected_b_ppi_e_strict_cont)

    expected_a_ppi_b_strict_cont = venn.Complement(venn.PropertySet(a_dash_c))
    assert strict_contingencies_state_set_a_ppi_b.is_equivalent_to(expected_a_ppi_b_strict_cont)

    expected_a_ppi_c_strict_cont = venn.UniversalSet()
    assert strict_contingencies_state_set_a_ppi_c.is_equivalent_to(expected_a_ppi_c_strict_cont)

    expected_b_pplus_e_strict_cont = venn.UniversalSet()
    assert strict_contingencies_state_set_b_pplus_e.is_equivalent_to(expected_b_pplus_e_strict_cont)


def test_set_of_states_from_contingencies_FOR_quant():
    # todo
    pass


def test_set_of_states_from_strict_contingencies_AND():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; AND A--C
                        <comp>; AND C--E
                        <comp>; AND B--F""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_f = rfs.state_from_string("B--F")

    rxncon = quick.rxncon_system

    strict_cont_state_set = rxncon.simulation.rule_based.rbm_from_rxncon.set_of_states_from_contingencies(rxncon.strict_contingencies_for_reaction(rfs.reaction_from_string('A_ppi_B')))

    assert strict_cont_state_set.is_equivalent_to(venn.Intersection(venn.Intersection(venn.PropertySet(expected_a_dash_c), venn.PropertySet(expected_c_dash_e)),
                                                                    venn.PropertySet(expected_b_dash_f)))


def test_source_set_of_states_for_reaction():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_dash_b = rfs.state_from_string('A--B')
    b_pplus_e = rfs.reaction_from_string('B_p+_E')
    e_pplus = rfs.state_from_string('E-{P}')
    b_pminus_e = rfs.reaction_from_string('B_p-_E')
    b_pt_e = rfs.reaction_from_string('B_pt_E')
    b_pplus = rfs.state_from_string('B-{P}')

    set_a_ppi_b = rxncon.simulation.rule_based.rbm_from_rxncon.source_state_set_from_reaction(a_ppi_b)
    set_b_pplus_e = rxncon.simulation.rule_based.rbm_from_rxncon.source_state_set_from_reaction(b_pplus_e)
    set_b_pminus_e = rxncon.simulation.rule_based.rbm_from_rxncon.source_state_set_from_reaction(b_pminus_e)
    set_b_pt_e = rxncon.simulation.rule_based.rbm_from_rxncon.source_state_set_from_reaction(b_pt_e)

    # todo: B_pt_E are two reactions in one B_p+_E -> E_[Bside] and E_p-_B -> B_[Eside]
    # todo: B_[n]_apt_B_[m] auto phosphortransfer B is the same molecule B_[n]_p+_B_[m] -> B_[m] and B_[m]_p-_B_[n] -> B_B[n]
    # todo: B_apt_B auto phosphortransfer B is the same molecule B_p+_B -> B_[Bsite1] and B_p-_B -> B_B[Site2]

    assert set_a_ppi_b.is_equivalent_to(venn.Complement(venn.PropertySet(a_dash_b)))
    assert set_b_pplus_e.is_equivalent_to(venn.Complement(venn.PropertySet(e_pplus)))
    assert set_b_pminus_e.is_equivalent_to(venn.PropertySet(e_pplus))
    assert set_b_pt_e.is_equivalent_to(venn.Intersection(venn.Complement(venn.PropertySet(e_pplus)), venn.PropertySet(b_pplus)))