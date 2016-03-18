import pytest

import rxncon.simulation.bBM.bBM_boolnet_exporter as bbe
import rxncon.simulation.bBM.bipartite_boolean_model as bbm

import rxncon.syntax.rxncon_from_string as rfs
import rxncon.venntastic.sets as venn


def test_generate_name():
    a_pplus_b = bbm.Node(rfs.reaction_from_string("a_p+_b"))
    assert bbe.string_from_reaction(a_pplus_b.value) == "a_pplus_b"

    a_dash_dash_b = bbm.Node(rfs.state_from_string("A--B"))
    assert bbe.string_from_inter_protein_interaction_state(a_dash_dash_b.value) == "A..B"

    b_intra = bbm.Node(rfs.state_from_string("b_[n]--[m]"))

    assert bbe.string_from_intra_protein_interaction_state(b_intra.value) == "b__n_.._m_"

    A_ppi_B = bbm.Node(rfs.reaction_from_string("A_[n]_ppi_B[d/s(r)]"))

    assert bbe.string_from_reaction(A_ppi_B.value) == "A__n__ppi_B_d_s_r__"


def test_test(rule_A__B, rule_A_ppi_B, rule_A_p, rule_C_pplus_A, initialConditions):
    bbm_system = bbm.Bipartite_Boolean_Model([rule_A__B, rule_A_ppi_B, rule_A_p, rule_C_pplus_A],
                                             initialConditions)

    bbe_system = bbe.BoolNet_System(bbm_system)
    bbe_str = bbe_system.to_string()

    expected_str = """target, factors
A, A
B, B
C, C
A..B, ((A_ppi_B & A.p) | (A..B & A.p))
A_ppi_B, ((A & B) & A.p)
A.p, (C_pplus_A | A.p)
C_pplus_A, (C & A)"""

    assert bbe_str == expected_str


@pytest.fixture
def rule_A__B():

    value_A__B = venn.Union(venn.Intersection(venn.PropertySet(bbm.Node(rfs.reaction_from_string("A_ppi_B"))),
                               venn.PropertySet(bbm.Node(rfs.state_from_string("A-{P}")))),
                       venn.Intersection(venn.PropertySet(bbm.Node(rfs.state_from_string("A--B"))),
                               venn.PropertySet(bbm.Node(rfs.state_from_string("A-{P}")))))

    return bbm.Rule(bbm.Node(rfs.state_from_string("A--B")), bbm.Factor(value_A__B))


@pytest.fixture
def rule_A_ppi_B():
    value_A_ppi_B = venn.Intersection(venn.Intersection(venn.PropertySet(bbm.Node(rfs.reaction_from_string("A_ppi_B").components[0])),
                                                            venn.PropertySet(bbm.Node(rfs.reaction_from_string("A_ppi_B").components[1]))),
                                          venn.PropertySet(bbm.Node(rfs.state_from_string("A-{P}"))))
    return bbm.Rule(bbm.Node(rfs.reaction_from_string("A_ppi_B")), bbm.Factor(value_A_ppi_B))


@pytest.fixture
def rule_A_p():
    value_A_p = venn.Union(venn.PropertySet(bbm.Node(rfs.reaction_from_string("C_p+_A"))),
                           venn.PropertySet(bbm.Node(rfs.state_from_string("A-{P}"))))

    return bbm.Rule(bbm.Node(rfs.state_from_string("A-{P}")), bbm.Factor(value_A_p))


@pytest.fixture
def rule_C_pplus_A():
    value_C_pplus_A = venn.Intersection(venn.PropertySet(bbm.Node(rfs.reaction_from_string("C_p+_A").components[0])),
                                        venn.PropertySet(bbm.Node(rfs.reaction_from_string("C_p+_A").components[1])))
    return bbm.Rule(bbm.Node(rfs.reaction_from_string("C_p+_A")), bbm.Factor(value_C_pplus_A))


@pytest.fixture
def initialConditions():
    return [bbm.InitConditions(bbm.Node(rfs.reaction_from_string("A_ppi_B").components[0]), None),
                      bbm.InitConditions(bbm.Node(rfs.reaction_from_string("A_ppi_B").components[1]), None),
                      bbm.InitConditions(bbm.Node(rfs.reaction_from_string("C_p+_A").components[0]), None)
           ]


