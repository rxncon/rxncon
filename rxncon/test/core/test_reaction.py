from rxncon.core.reaction import *


def test_simple():
    rxn = reaction_from_string(REACTION_DEFINITIONS, 'A_[x]_ppi_B_[y]')
    print(rxn)
    print(rxn.reactants_pre)
    print(rxn.reactants_post)


