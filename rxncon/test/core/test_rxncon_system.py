import pytest

import rxncon.syntax.rxncon_from_string as rfs
import rxncon.core.contingency as con
import rxncon.core.effector as eff
import rxncon.core.rxncon_system as rxs


def test_simple_rxncon_system(simple_system):
    pass



@pytest.fixture
def simple_system():
    phosphorylation_reaction = rfs.reaction_from_string('A_p+_B')
    binding_reaction = rfs.reaction_from_string('B_ppi_C')

    phosphorylated_state = rfs.state_from_string('B-{p}')

    binding_contingency = con.Contingency(binding_reaction,
                                          con.ContingencyType.requirement,
                                          eff.StateEffector(phosphorylated_state))

    reactions = [phosphorylation_reaction, binding_reaction]
    contingencies = [binding_contingency]

    return rxs.RxnConSystem(reactions, contingencies)
