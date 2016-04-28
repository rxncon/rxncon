import pytest
import rxncon.simulation.bBM.bbm_from_rxncon as bfr
import rxncon.input.quick.quick as quick
import rxncon.simulation.bBM.bBM_boolnet_exporter as bbe


def test_rule():
    rxncon_sys = quick.Quick("""A_ppi_B; ! <comp>
                               <comp>; AND A--C; AND A--D; AND A-{p}
                               A_ppi_C
                               A_ppi_D
                               C_p+_A
                               D_p-_A
                               E_p-_A
                            """)
    bbm_sys = bfr.bipartite_boolean_model_from_rxncon(rxncon_sys.rxncon_system)
    bbe_system = bbe.BoolNetSystem(bbm_sys)
    expected_str = """target, factors
A, A
B, B
C, C
D, D
E, E
A_ppi_B, (((A__C & A__D) & A_.p.) & (A & B))
A__B, (A__B | A_ppi_B)
A_ppi_C, (A & C)
A__C, (A__C | A_ppi_C)
A_ppi_D, (A & D)
A__D, (A__D | A_ppi_D)
C_pplus_A, (C & A)
A_.p., (((A_.p. & ! D_pminus_A) & ! E_pminus_A) | ((C_pplus_A & ! D_pminus_A) & ! E_pminus_A))
D_pminus_A, (D & A)
E_pminus_A, (E & A)"""

    assert bbe_system.to_string() == expected_str


def test_convert_quantitative_contingencies_into_strict_contingencies():
    rxncon_sys = quick.Quick("""A_ppi_B; k- A-{P}; k+ A--C
                               A_ppi_C
                               D_p+_A""")
    bbm_sys = bfr.bipartite_boolean_model_from_rxncon(rxncon_sys.rxncon_system)
    bbe_system = bbe.BoolNetSystem(bbm_sys)
    expected_expression = """target, factors
A, A
B, B
C, C
D, D
A_ppi_B, ((A & B) & (A__C & ! A_.p.))
A__B, (A__B | A_ppi_B)
A_ppi_C, (A & C)
A__C, (A__C | A_ppi_C)
D_pplus_A, (D & A)
A_.p., (A_.p. | D_pplus_A)"""
    assert bbe_system.to_string() == expected_expression


def test_complex_multiple_boolean_expression():
    quick_sys = quick.Quick("""
    A_trsc_B; ! <comp1>; x <comp2>; x <comp3>; x <comp4>
    <comp1>; AND C--D; AND C--E
    <comp2>; AND G--F; AND <comp4>
    <comp3>;AND H--F; AND <comp4>
    <comp4>; AND F--E; AND <comp1>
    C_ppi_D
    C_ppi_E
    G_ppi_F
    H_ppi_F
    F_ppi_E
    """)

# Original by Basti:
#     expected_str = """target, factors
# A, A
# A_trsc_B, (((C__D & C__E) & ! F__E) & (A & B))
# B, (B | A_trsc_B)"""

    expected_str= """target, factors
A, A
C, C
D, D
E, E
G, G
F, F
H, H
A_trsc_B, (((C__D & C__E) & ! F__E) & (A & B))
B, (B | A_trsc_B)
C_ppi_D, (C & D)
C__D, (C__D | C_ppi_D)
C_ppi_E, (C & E)
C__E, (C__E | C_ppi_E)
G_ppi_F, (G & F)
G__F, (G__F | G_ppi_F)
H_ppi_F, (H & F)
H__F, (H__F | H_ppi_F)
F_ppi_E, (F & E)
F__E, (F__E | F_ppi_E)"""

    bbm_sys = bfr.bipartite_boolean_model_from_rxncon(quick_sys.rxncon_system)
    bbe_system = bbe.BoolNetSystem(bbm_sys)
    assert bbe_system.to_string() == expected_str


def test_input_output_system():
    rxncon_sys = quick.Quick("""A_ppi_B; ! A-{P}; ! [Input]
                               C_p+_A
                               [Output]; ! A--B
                            """)

    bbm_sys = bfr.bipartite_boolean_model_from_rxncon(rxncon_sys.rxncon_system)
    bbe_system = bbe.BoolNetSystem(bbm_sys)

    expected_str = """target, factors
A, A
B, B
C, C
.Input., .Input.
A__B, (A_ppi_B | A__B)
A_ppi_B, ((A & B) & A_.p.)
A_.p., (C_pplus_A | A_.p.)
C_pplus_A, (C & A)
.Output., (.Output. | A__B)"""

    assert bbe_system.to_string() == expected_str


def test_contradictory_expression():
    quick_sys = quick.Quick("""C_p+_A
                               A_ppi_B; x A-{p}; ! A-{P}""")
    with pytest.raises(AssertionError):
        bfr.bipartite_boolean_model_from_rxncon(quick_sys.rxncon_system)


def test_rule_validation():
    with pytest.raises(AssertionError):
        rxncon_sys = quick.Quick("""A_ppi_B; ! A-{P}; ! A--C
                                   D_p+_A""")
        bfr.bipartite_boolean_model_from_rxncon(rxncon_sys.rxncon_system)
