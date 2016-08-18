from rxncon.input.quick.quick import *
from rxncon.core.state import state_from_string

def test_simple_quick():
    rxncon_sys = Quick('A_p+_B_[(r)]').rxncon_system
    assert rxncon_sys.product_states == [state_from_string('B_[(r)]-{p}')]


def test_with_cont():
    rxncon_sys = Quick('''A_p+_B_[(r)]; ! A_[(x)]-{p}
                       C_p+_A_[(x)]''').rxncon_system
    assert rxncon_sys.product_states == [state_from_string('A_[(x)]-{p}'), state_from_string('B_[(r)]-{p}')]


def test_ppi():
    rxncon_sys = Quick('''A_[x]_ppi_B_[y]''').rxncon_system

    print(rxncon_sys)
