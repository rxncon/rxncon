from rxncon.core.reaction import *


def test_simple():
    rxn = reaction_from_string(REACTION_DEFINITIONS, 'A_p+_B_[(y)]')
    print(rxn)
    print(rxn.reactants_pre)
    print(rxn.reactants_post)


