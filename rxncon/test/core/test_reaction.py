from rxncon.core.reaction import *


def test_simple():
    rxn = reaction_from_string(REACTION_DEFINITIONS, 'A_trsc_B')
    print(rxn)
    print(rxn.reactants_pre)
    print(rxn.reactants_post)


