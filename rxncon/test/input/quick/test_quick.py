from rxncon.input.quick.quick import *

def test_simple_quick():
    rxncon_sys = Quick('A_p+_B_[(r)]').rxncon_system
    print(rxncon_sys.reactions)


def test_with_cont():
    rxncon_sys = Quick('''A_p+_B_[(r)]; ! A_[(x)]-{p}
                       C_p+_A_[(x)]''').rxncon_system
    print(rxncon_sys.reactions)
