import typing as tg
import rxncon.core.rxncon_system as rxs
import rxncon.core.reaction as rxn
import rxncon.simulation.bBM.bipartite_boolean_model as bbm
import rxncon.venntastic.sets as venn

from rxncon.semantics.molecule_from_rxncon import set_of_states_from_contingencies


def rules_from_rxncon(rxconsys: rxs.RxnConSystem):

    rules = []
    for reaction in rxconsys.reactions:
        rules.append(rule_for_reaction_from_rxnconsys_and_reaction(rxconsys, reaction))

        rules.append(rule_for_state_from_rxnconsys_and_reaction(rxconsys, reaction, rules))
    return rules


def rule_for_reaction_from_rxnconsys_and_reaction(rxnconsys: rxs.RxnConSystem, reaction: rxn.Reaction) -> bbm.Rule:

    strict_contingency_state_set = set_of_states_from_contingencies(rxnconsys.strict_contingencies_for_reaction(reaction))
    vennset= venn.Intersection(strict_contingency_state_set,
                                                    venn.Intersection(venn.PropertySet(reaction.subject),
                                                                      venn.PropertySet(reaction.object)))
    quantitative_contingency_state_set= set_of_states_from_contingencies(rxnconsys.quantitative_contingencies_for_reaction(reaction)) # todo: ersetze k+ durch ! und k- durch x. in venntastic: ! ^= PropertySet(), x ^= Complement(PropertySet())
    return bbm.Rule(bbm.Node(reaction), bbm.Factor(vennset_to_bbm_factor_vennset(vennset.simplified_form())))



def vennset_to_bbm_factor_vennset(vennset: venn.Set):
    # creates new vennset with states contained by Node objects, for compareability
    # want to rewrite venn.Set into bbm.Factor like venn.PropertySet(A--B) -> venn.PropertySet(bbm.Node(A--B))
    if isinstance(vennset, venn.PropertySet):
        return venn.PropertySet(bbm.Node(vennset.value))
    #elif isinstance(vennset, venn.UniversalSet()):
     #   return venn.PropertySet(bbm.Node(vennset.value))
    elif isinstance(vennset, venn.Intersection):
        return venn.Intersection(vennset_to_bbm_factor_vennset(vennset.left_expr), vennset_to_bbm_factor_vennset(vennset.right_expr))
    elif isinstance(vennset, venn.Union):
        return venn.Union(vennset_to_bbm_factor_vennset(vennset.left_expr), vennset_to_bbm_factor_vennset(vennset.right_expr))
    else:
        raise NotImplementedError


# def reaction_dependencies_for_states(rxnconsys: rxs.RxnConSystem, reaction: rxn.Reaction):
#     strict_contingency_state_set = set_of_states_from_contingencies(rxnconsys.strict_contingencies_for_reaction(reaction))
#     reaction_dependency = venn.Intersection(venn.PropertySet(bbm.Node(reaction)), vennset_to_bbm_factor_vennset(strict_contingency_state_set))
#
#     return reaction_dependency
def get_rule_targets(rules: tg.List[bbm.Rule]):

    all_visited_states = [rule.target.value for rule in rules]
    return all_visited_states

def rule_for_state_from_rxnconsys_and_reaction(rxnconsys: rxs.RxnConSystem, reaction: rxn.Reaction, system_rules: tg.List[bbm.Rule]) -> bbm.Rule:
    all_visited_states = get_rule_targets(system_rules)
    rules = []
    if reaction.product is None and reaction.product not in all_visited_states:
        return rules

    pos_bool_def=[venn.PropertySet(bbm.Node(reaction.product)), venn.PropertySet(bbm.Node(reaction))]
    neg_bool_def=[]

    for rxn in rxnconsys.reactions:
        if rxn != reaction and reaction.product in [rxn.product]:
            pos_bool_def.append(venn.PropertySet(bbm.Node(rxn)))
        if rxn != reaction and reaction.product in [rxn.source]:
            neg_bool_def.append(venn.PropertySet(bbm.Node(rxn)))


    pos_rules= venn.nested_expression_from_list_and_binary_op(pos_bool_def, venn.Union)
    neg_rules = venn.nested_expression_from_list_and_binary_op(neg_bool_def, venn.Union)

    bool_rules= venn.Intersection(pos_rules, venn.Complement(neg_rules))
    # bool_rules.simplified_form() InterProteinInteraction has no _complement_expanded
    rules.append(bbm.Rule(bbm.Node(reaction.product), bbm.Factor(bool_rules)))

    return rules


