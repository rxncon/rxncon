from rxncon.venntastic.sets import PropertySet, Intersection, Union, Complement

from rxncon.syntax.rxncon_from_string import state_from_string, specification_from_string, reaction_from_string
from rxncon.simulation.rule_based.molecule_from_string import mol_def_from_string, mol_instance_from_string
from rxncon.simulation.rule_based.rbm_from_rxncon import mol_instance_set_from_state_set, reacting_mol_instance_lhs_rhs_sets_from_reaction


def test():
    mol_defs = {
        specification_from_string('A'): mol_def_from_string('A#ass/A_[Bassoc]:B_[Aassoc],mod/A_[(r)]:u~p~ub'),
        specification_from_string('B'): mol_def_from_string('B#ass/B_[Aassoc]:A_[Bassoc]')
    }
    x = Complement(PropertySet(state_from_string('A--B')))
    y = Complement(PropertySet(state_from_string('A-{p}')))
    z = Intersection(x, y)

    print()

    for x in reacting_mol_instance_lhs_rhs_sets_from_reaction(mol_defs, reaction_from_string('A_ppi_B')):
        print()
        print(x)

    # lhss = mol_instance_set_from_state_set(mol_defs, z).to_nested_list_form()
    #
    # for lhs in lhss:
    #     print()
    #     print(lhs)
