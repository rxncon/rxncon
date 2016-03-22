import rxncon.simulation.bBM.bbm_from_rxncon as bfr
import rxncon.input.quick.quick as quick
import rxncon.venntastic.sets as venn
import rxncon.simulation.bBM.bBM_boolnet_exporter as bbe


def test_rule():
    rxncon_sys= quick.Quick("""A_ppi_B; ! A-{P}
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
A_ppi_B, (A._p_ & (A & B))
A..B, (A..B | A_ppi_B)
C_pplus_A, (C & A)
A._p_, ((A._p_ | C_pplus_A) & ! (D_pminus_A | E_pminus_A))
D_pminus_A, (D & A)
E_pminus_A, (E & A)"""
    assert bbe_system.to_string() == expected_str
