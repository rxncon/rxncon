from rxncon.core.reaction import *

def test_simple_reaction():
    rxn = reaction_from_string(REACTION_DEFINITIONS, 'A_p+_B_[(r)]')
    print(rxn.reactants_post)

    assert rxn.reactants_post == [Reactant(spec_from_string('A'), None),
                                  Reactant(spec_from_string('B'), [state_from_string(STATE_DEFS, 'B_[(r)]-{p}')])]
