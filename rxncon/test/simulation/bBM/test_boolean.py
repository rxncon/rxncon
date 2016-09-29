import pytest
import rxncon.input.quick.quick as quick
import rxncon.simulation.bBM.boolean_model as bm


def test_simple_system():
    rxncon_sys = quick.Quick("""A_[b]_ppi_B_[a]; ! A_[(r)]-{p}
                                C_p+_A_[(r)]
                                D_p-_A_[(r)]""")
    boolmodel = bm.boolean_model_from_rxncon(rxncon_sys.rxncon_system)

