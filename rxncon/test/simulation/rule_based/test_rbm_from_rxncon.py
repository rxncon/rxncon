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

    spec_set = mol_def_X.specification_set_from_state_set(lhs_set)

    spec_sets_overlapping = spec_set.to_union_list_form()
    spec_sets_disjunct = venn.gram_schmidt_disjunctify(spec_sets_overlapping)

    for ss in spec_sets_disjunct:
        print(mol_def_X.specification_from_specification_set(ss))


def test_single_rule():
    c_phos_a = rfs.reaction_from_string('C_p+_A_[x]')
    a_phos = rfs.state_from_string('A_[x]-{p}')
    a_phos_b = rfs.reaction_from_string('A_p+_B')

    cont = con.Contingency(a_phos_b, con.ContingencyType.requirement, eff.StateEffector(a_phos))

    rxncon = rxs.RxnConSystem([c_phos_a, a_phos_b], [cont])

    # print(rfr.molecule_defs_from_rxncon(rxncon)['B'].modification_defs[0].matching_state)

    rules = rfr.rules_from_rxncon(rxncon)

    for rule in rules:
        print()
        print(rule)



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
