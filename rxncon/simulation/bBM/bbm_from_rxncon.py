import typing as tg
import rxncon.core.rxncon_system as rxs
import rxncon.core.reaction as rxn
import rxncon.simulation.bBM.bipartite_boolean_model as bbm
import rxncon.venntastic.sets as venn

from rxncon.semantics.molecule_from_rxncon import set_of_states_from_contingencies


def rules_from_rxncon(rxconsys: rxs.RxnConSystem):

    rules = []
    for reaction in rxconsys.reactions:
        rules += rule_for_reaction_from_rxnconsys_and_reaction(rxconsys, reaction)

        rules += rule_for_state_from_rxnconsys_and_reaction(rxconsys, reaction)
    return rules


def rule_for_reaction_from_rxnconsys_and_reaction(rxnconsys: rxs.RxnConSystem, reaction: rxn.Reaction) -> bbm.Rule:

    strict_contingency_state_set = set_of_states_from_contingencies(rxnconsys.strict_contingencies_for_reaction(reaction))


    return bbm.Rule(bbm.Node(reaction), bbm.Factor(strict_contingency_state_set))



def vennset_to_bbm_factor_vennset(vennset: venn.Set):
    # want to rewrite venn.Set into bbm.Factor like venn.PropertySet(A--B) -> venn.PropertySet(bbm.Node(A--B))
    if isinstance(vennset, venn.PropertySet):
        return venn.PropertySet(bbm.Node(vennset.value))
    elif isinstance(vennset, venn.Intersection):
        return venn.Intersection(vennset_to_bbm_factor_vennset(vennset.left_expr), vennset_to_bbm_factor_vennset(vennset.right_expr))
    elif isinstance(vennset, venn.Union):
        return venn.Union(vennset_to_bbm_factor_vennset(vennset.left_expr), vennset_to_bbm_factor_vennset(vennset.right_expr))
    else:
        raise NotImplementedError


def reaction_dependencies(rxnconsys: rxs.RxnConSystem, reaction: rxn.Reaction):
    strict_contingency_state_set = set_of_states_from_contingencies(rxnconsys.strict_contingencies_for_reaction(reaction))
    reaction_dependency = venn.Intersection(venn.PropertySet(bbm.Node(reaction)), vennset_to_bbm_factor_vennset(strict_contingency_state_set))

    return reaction_dependency

def rule_for_state_from_rxnconsys_and_reaction(rxnconsys: rxs.RxnConSystem, reaction: rxn.Reaction) -> bbm.Rule:

    rules = []
    boolean_list = []
    reaction_states = [state for state in [reaction.source, reaction.product] if state is not None]
    for state in reaction_states:
        boolean_list.append(state)
        for rxn in rxnconsys.reactions:
            if state in [rxn.source, rxn.source]:
                boolean_list.append(reaction_dependencies(rxnconsys, rxn))

        rules.append(bbm.Rule(bbm.Node(state), bbm.Factor(venn.nested_expression_from_list_and_binary_op(boolean_list, venn.Union))))


