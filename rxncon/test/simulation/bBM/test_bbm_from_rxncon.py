import rxncon.simulation.bBM.bbm_from_rxncon as bfr
import rxncon.input.quick.quick as quick
import rxncon.venntastic.sets as venn


def test_rule():
    rxncon_sys= quick.Quick("A_ppi_B")
    bbm_sys = bfr.rules_from_rxncon(rxncon_sys.rxncon_system)

    pass


