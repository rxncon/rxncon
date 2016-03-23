import rxncon.simulation.bBM.bbm_from_rxncon as bfr
import rxncon.input.quick.quick as quick
import rxncon.simulation.bBM.bBM_boolnet_exporter as bbe


def test_rule():
    rxncon_sys= quick.Quick("""A_ppi_B; ! <comp>
                               <comp>; AND A--C
                               <comp>; AND A--D
                               <comp>; AND A-{p}
                               A_ppi_C
                               A_ppi_D
                               C_p+_A
                               D_p-_A
                               E_p-_A
                            """)
    bbm_sys = bfr.bipartite_boolean_model_from_rxncon(rxncon_sys.rxncon_system)
    bbe_system = bbe.BoolNet_System(bbm_sys)
    expected_str = """target, factors
A, A
B, B
C, C
D, D
E, E
A_ppi_B, (((A..C & A..D) & A._p_) & (A & B))
A..B, (A..B | A_ppi_B)
A_ppi_C, (A & C)
A..C, (A..C | A_ppi_C)
A_ppi_D, (A & D)
A..D, (A..D | A_ppi_D)
C_pplus_A, (C & A)
A._p_, ((A._p_ & (! D_pminus_A & ! E_pminus_A)) | (C_pplus_A & (! D_pminus_A & ! E_pminus_A)))
D_pminus_A, (D & A)
E_pminus_A, (E & A)"""
    assert bbe_system.to_string() == expected_str


def test_test():
    rxncon_sys= quick.Quick("""A_ppi_B; k- A-{P}; k+ A--C
                               A_ppi_C
                               D_p+_A""")
    bbm_sys = bfr.bipartite_boolean_model_from_rxncon(rxncon_sys.rxncon_system)
    bbe_system = bbe.BoolNet_System(bbm_sys)
    expected_expression = """target, factors
A, A
B, B
C, C
D, D
A_ppi_B, ((A & B) & (A..C & ! A._p_))
A..B, (A..B | A_ppi_B)
A_ppi_C, (A & C)
A..C, (A..C | A_ppi_C)
D_pplus_A, (D & A)
A._p_, (A._p_ | D_pplus_A)"""
    assert bbe_system.to_string() == expected_expression


