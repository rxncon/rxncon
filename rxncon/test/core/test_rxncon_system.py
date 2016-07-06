import pytest

import rxncon.syntax.rxncon_from_string as rfs
import rxncon.core.reaction as rxn
import rxncon.core.contingency as con
import rxncon.core.effector as eff
import rxncon.core.rxncon_system as rxs


def test_simple_rxncon_system(simple_system):
    pass



@pytest.fixture
def simple_system():
    phosphorylation_reaction = rxn.reaction_from_string('A_p+_B_[(r)]')
    binding_reaction = rxn.reaction_from_string('B_[c]_ppi_C_[b]')

    phosphorylated_state = rxn.state_from_string('B_[(r)]-{p}')

    binding_contingency = con.Contingency(binding_reaction,
                                          con.ContingencyType.requirement,
                                          eff.StateEffector(phosphorylated_state))

    reactions = [phosphorylation_reaction, binding_reaction]
    contingencies = [binding_contingency]

    return rxs.RxnConSystem(reactions, contingencies)
