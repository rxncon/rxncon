import os
import time
import tempfile
import pytest

import rxncon.simulation.boolean.bipartite_boolean_model as bbm
import rxncon.simulation.boolean.bBM_boolnet_exporter as bbe
import rxncon.syntax.rxncon_from_string as rfs
import rxncon.venntastic.sets as venn
import rxncon.core.state as sta



def test_generate_name():
    a_pplus_b = bbm.Node(rfs.reaction_from_string("a_p+_b"))
    assert bbe.string_from_reaction(a_pplus_b.value) == "a_pplus_b"

    a_dash_dash_b = bbm.Node(rfs.state_from_string("A--B"))
    assert bbe.replace_invalid_chars(str(a_dash_dash_b.value)) == "A__B"

    b_intra = bbm.Node(rfs.state_from_string("b_[n]--[m]"))

    assert bbe.replace_invalid_chars(str(b_intra.value)) == "b_n__m"

    A_ppi_B = bbm.Node(rfs.reaction_from_string("A_[n]_ppi_B_[d/s(r)]"))

    assert bbe.string_from_reaction(A_ppi_B.value) == "A_n_ppi_B_dsr"

    A_p = bbm.Node(rfs.state_from_string("A_[d/s(r)]-{p}"))
    assert bbe.replace_invalid_chars(str(A_p)) == 'A_dsr_p'

    A_p = bbm.Node(rfs.state_from_string("A-{p}"))
    assert bbe.replace_invalid_chars(str(A_p)) == 'A_p'

    output = bbm.Node(rfs.reaction_from_string('[Output]'))
    assert bbe.replace_invalid_chars(str(output)) == "_Output_"



def test_boolnet_string(rule_A__B, rule_A_ppi_B, rule_A_p, rule_C_pplus_A, initialConditions_boolnet_string):
    bbm_system = bbm.BipartiteBooleanModel([rule_A__B, rule_A_ppi_B, rule_A_p, rule_C_pplus_A],
                                           initialConditions_boolnet_string)

    bbe_system = bbe.BoolNetSystem(bbm_system)
    bbe_str = bbe_system.to_string()

    expected_str = """target, factors
A, A
B, B
C, C
A__B, A_ppi_B
A_ppi_B, ((A & B) & A_p)
A_p, C_pplus_A
C_pplus_A, (C & A)"""

    assert bbe_str == expected_str

def test_boolnet_input_output_string(rule_A__B, rule_A_ppi_B, rule_A_p, rule_C_pplus_A, rule_output, initial_conditions_boolnet_input_output_string):
    bbm_system = bbm.BipartiteBooleanModel([rule_A__B, rule_A_ppi_B, rule_A_p, rule_C_pplus_A, rule_output],
                                           initial_conditions_boolnet_input_output_string)

    bbe_system = bbe.BoolNetSystem(bbm_system)
    bbe_str = bbe_system.to_string()

    expected_str = """target, factors
A, A
B, B
C, C
_Input_, _Input_
A__B, A_ppi_B
A_ppi_B, ((A & B) & A_p)
A_p, C_pplus_A
C_pplus_A, (C & A)
_Output_, (_Output_ | A__B)"""

    assert bbe_str == expected_str

def test_boolnet_string_with_complement(rule_A__B, rule_A_ppi_B, rule_A_p_deg, rule_C_pplus_A, rule_D_pminus_A, rule_E_pminus_A, initialConditions):
    bbm_system = bbm.BipartiteBooleanModel([rule_A__B, rule_A_ppi_B, rule_A_p_deg, rule_C_pplus_A, rule_D_pminus_A, rule_E_pminus_A],
                                           initialConditions)

    bbe_system = bbe.BoolNetSystem(bbm_system)
    bbe_str = bbe_system.to_string()
    expected_str = """target, factors
A, A
B, B
C, C
D, D
E, E
A__B, A_ppi_B
A_ppi_B, ((A & B) & A_p)
A_p, (C_pplus_A | (A_p & (! D_pminus_A & ! E_pminus_A)))
C_pplus_A, (C & A)
D_pminus_A, (D & A)
E_pminus_A, (E & A)"""

    assert bbe_str == expected_str


def test_to_file(rule_A__B, rule_A_ppi_B, rule_A_p_deg, rule_C_pplus_A, rule_D_pminus_A, rule_E_pminus_A, initialConditions):
    bbm_system = bbm.BipartiteBooleanModel([rule_A__B, rule_A_ppi_B, rule_A_p_deg, rule_C_pplus_A, rule_D_pminus_A, rule_E_pminus_A], initialConditions)

    bbe_system = bbe.BoolNetSystem(bbm_system)
    bbe_str = bbe_system.to_string()
    path = "{0}/test{1}.bool".format(tempfile.gettempdir(), time.time())
    bbe_system.to_file(path)
    assert os.path.exists(path)
    os.remove(path)

@pytest.fixture
def rule_output():
    value_ouput = venn.Union(venn.ValueSet(bbm.Node(rfs.reaction_from_string('[Output]'))),
                             venn.ValueSet(bbm.Node(rfs.state_from_string("A--B"))))
    return bbm.Rule(bbm.Node(rfs.reaction_from_string('[Output]')), bbm.Factor(value_ouput))

@pytest.fixture
def rule_A__B():

    value_A__B = venn.ValueSet(bbm.Node(rfs.reaction_from_string("A_ppi_B")))

    return bbm.Rule(bbm.Node(rfs.state_from_string("A--B")), bbm.Factor(value_A__B))


@pytest.fixture
def rule_A_ppi_B():
    value_A_ppi_B = venn.Intersection(venn.Intersection(venn.ValueSet(bbm.Node(sta.ComponentState(rfs.reaction_from_string("A_ppi_B").components_lhs[0].to_component_spec()))),
                                                        venn.ValueSet(bbm.Node(sta.ComponentState(rfs.reaction_from_string("A_ppi_B").components_lhs[1].to_component_spec())))),
                                      venn.ValueSet(bbm.Node(rfs.state_from_string("A-{P}"))))
    return bbm.Rule(bbm.Node(rfs.reaction_from_string("A_ppi_B")), bbm.Factor(value_A_ppi_B))


@pytest.fixture
def rule_A_p():
    value_A_p = venn.ValueSet(bbm.Node(rfs.reaction_from_string("C_p+_A")))


    return bbm.Rule(bbm.Node(rfs.state_from_string("A-{P}")), bbm.Factor(value_A_p))


@pytest.fixture
def rule_A_p_deg():
    value_A_p_deg = venn.Union(venn.ValueSet(bbm.Node(rfs.reaction_from_string("C_p+_A"))),
                               venn.Intersection(venn.ValueSet(bbm.Node(rfs.state_from_string("A-{P}"))),
                                                 venn.Intersection(venn.Complement(venn.ValueSet(bbm.Node(rfs.reaction_from_string("D_p-_A")))),
                                                                   venn.Complement(venn.ValueSet(bbm.Node(rfs.reaction_from_string("E_p-_A"))))
                                                                   )
                                                 )
                               )


    return bbm.Rule(bbm.Node(rfs.state_from_string("A-{P}")), bbm.Factor(value_A_p_deg))


@pytest.fixture
def rule_C_pplus_A():
    value_C_pplus_A = venn.Intersection(venn.ValueSet(bbm.Node(sta.ComponentState(rfs.reaction_from_string("C_p+_A").components_lhs[0].to_component_spec()))),
                                        venn.ValueSet(bbm.Node(sta.ComponentState(rfs.reaction_from_string("C_p+_A").components_lhs[1].to_component_spec()))))
    return bbm.Rule(bbm.Node(rfs.reaction_from_string("C_p+_A")), bbm.Factor(value_C_pplus_A))


@pytest.fixture
def rule_D_pminus_A():
    value_D_pminus_A = venn.Intersection(venn.ValueSet(bbm.Node(sta.ComponentState(rfs.reaction_from_string("D_p-_A").components_lhs[0].to_component_spec()))),
                                         venn.ValueSet(bbm.Node(sta.ComponentState(rfs.reaction_from_string("D_p-_A").components_lhs[1].to_component_spec()))))
    return bbm.Rule(bbm.Node(rfs.reaction_from_string("D_p-_A")), bbm.Factor(value_D_pminus_A))


@pytest.fixture
def rule_E_pminus_A():
    value_E_pminus_A = venn.Intersection(venn.ValueSet(bbm.Node(sta.ComponentState(rfs.reaction_from_string("E_p-_A").components_lhs[0].to_component_spec()))),
                                         venn.ValueSet(bbm.Node(sta.ComponentState(rfs.reaction_from_string("E_p-_A").components_lhs[1].to_component_spec()))))

    return bbm.Rule(bbm.Node(rfs.reaction_from_string("E_p-_A")), bbm.Factor(value_E_pminus_A))

@pytest.fixture
def initial_conditions_boolnet_input_output_string():
    return [bbm.InitCondition(bbm.Node(sta.ComponentState(rfs.reaction_from_string("A_ppi_B").components_lhs[0].to_component_spec())), None),
            bbm.InitCondition(bbm.Node(sta.ComponentState(rfs.reaction_from_string("A_ppi_B").components_lhs[1].to_component_spec())), None),
            bbm.InitCondition(bbm.Node(sta.ComponentState(rfs.reaction_from_string("C_p+_A").components_lhs[0].to_component_spec())), None),
            bbm.InitCondition(bbm.Node(rfs.state_from_string("[Input]")),None)
            ]

@pytest.fixture
def initialConditions_boolnet_string():
    return [bbm.InitCondition(bbm.Node(sta.ComponentState(rfs.reaction_from_string("A_ppi_B").components_lhs[0].to_component_spec())), None),
            bbm.InitCondition(bbm.Node(sta.ComponentState(rfs.reaction_from_string("A_ppi_B").components_lhs[1].to_component_spec())), None),
            bbm.InitCondition(bbm.Node(sta.ComponentState(rfs.reaction_from_string("C_p+_A").components_lhs[0].to_component_spec())), None)
            ]


@pytest.fixture
def initialConditions():
    return [bbm.InitCondition(bbm.Node(sta.ComponentState(rfs.reaction_from_string("A_ppi_B").components_lhs[0].to_component_spec())), None),
                      bbm.InitCondition(bbm.Node(sta.ComponentState(rfs.reaction_from_string("A_ppi_B").components_lhs[1].to_component_spec())), None),
                      bbm.InitCondition(bbm.Node(sta.ComponentState(rfs.reaction_from_string("C_p+_A").components_lhs[0].to_component_spec())), None),
                      bbm.InitCondition(bbm.Node(sta.ComponentState(rfs.reaction_from_string("D_p-_A").components_lhs[0].to_component_spec())), None),
                      bbm.InitCondition(bbm.Node(sta.ComponentState(rfs.reaction_from_string("E_p-_A").components_lhs[0].to_component_spec())), None)
           ]


