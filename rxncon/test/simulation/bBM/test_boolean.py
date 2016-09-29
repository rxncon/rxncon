import pytest
import rxncon.input.quick.quick as quick
import rxncon.simulation.bBM.boolean_model as bm


def test_simple_system():
    rxncon_sys = quick.Quick("""A_[b]_ppi+_B_[a]; ! A_[(r)]-{p}
                                A_[b]_ppi-_B_[a]
                                C_p+_A_[(r)]
                                D_p-_A_[(r)]""")
    boolmodel = bm.boolean_model_from_rxncon(rxncon_sys.rxncon_system)
    for rule in boolmodel.update_rules:
        print(rule)
    boolmodel

