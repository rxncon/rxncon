import pytest

import rxncon.core.rxncon_system as rxs
import rxncon.syntax.rxncon_from_string as rfs
import rxncon.core.contingency as con
import rxncon.core.effector as eff

import rxncon.simulation.rule_based.rbm_from_rxncon as rfr
import rxncon.simulation.rule_based.rule_based_model as rbm
import rxncon.venntastic.sets as venn


def test_simple_rxncon_system(simple_system):
    mol_defs = rfr.molecule_defs_from_rxncon(simple_system)

    mol_def_X = mol_defs['X']

    lhs_set = venn.Intersection(venn.PropertySet(rfs.state_from_string('X-{p}')),
                                venn.Complement(venn.PropertySet(rfs.state_from_string('X--Y'))))

    matches = mol_def_X.match_with_state_set(lhs_set)

    for match in matches:
        print()
        print(match)


@pytest.fixture
def simple_system():
    phosphorylation_reaction_1 = rfs.reaction_from_string('A_p+_X')
    phosphorylation_reaction_2 = rfs.reaction_from_string('B_p+_X')
    phosphorylation_reaction_3 = rfs.reaction_from_string('C_p+_X')

    binding_reaction = rfs.reaction_from_string('X_ppi_Y')
    phosphorylated_state = rfs.state_from_string('X-{p}')

    binding_contingency = con.Contingency(binding_reaction,
                                          con.ContingencyType.requirement,
                                          eff.StateEffector(phosphorylated_state))

    reactions = [phosphorylation_reaction_1, phosphorylation_reaction_2, phosphorylation_reaction_3, binding_reaction]
    contingencies = [binding_contingency]

    return rxs.RxnConSystem(reactions, contingencies)
