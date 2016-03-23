import typing as tg
import rxncon.core.rxncon_system as rxs
import rxncon.core.reaction as rxn
import rxncon.core.contingency as con
import rxncon.simulation.bBM.bipartite_boolean_model as bbm
import rxncon.venntastic.sets as venn

from rxncon.simulation.rule_based.rbm_from_rxncon import state_set_from_contingencies


def bipartite_boolean_model_from_rxncon(rxconsys: rxs.RxnConSystem):
    return bbm.Bipartite_Boolean_Model(rules_from_rxncon(rxconsys), initial_states_from_rxncon(rxconsys))


def rules_from_rxncon(rxconsys: rxs.RxnConSystem):

    rules = []
    for reaction in rxconsys.reactions:
        rules.append(rule_for_reaction_from_rxnconsys_and_reaction(rxconsys, reaction, rules))

        rules.append(rule_for_state_from_rxnconsys_and_reaction(rxconsys, reaction, rules))
    rules = [rule for rule in rules if rule]
    return rules


def initial_states_from_rxncon(rxconsys: rxs.RxnConSystem):
    initial_states = []
    for reaction in rxconsys.reactions:
        if bbm.InitConditions(bbm.Node(reaction.subject), None) not in initial_states:
            initial_states.append(bbm.InitConditions(bbm.Node(reaction.subject), None))
        if bbm.InitConditions(bbm.Node(reaction.object), None) not in initial_states:
            initial_states.append(bbm.InitConditions(bbm.Node(reaction.object), None))
    return initial_states


def rule_for_reaction_from_rxnconsys_and_reaction(rxnconsys: rxs.RxnConSystem, reaction: rxn.Reaction,
                                                  system_rules: tg.List[bbm.Rule]) -> bbm.Rule:
    all_visited_nodes = get_rule_targets(system_rules)
    if bbm.Node(reaction) in all_visited_nodes:
        return None
    strict_contingency_state_set = state_set_from_contingencies(rxnconsys.strict_contingencies_for_reaction(reaction))

    vennset = venn.Intersection(strict_contingency_state_set,
                                   venn.Intersection(venn.PropertySet(reaction.subject),
                                                     venn.PropertySet(reaction.object)))
    additional_strict_cont = convert_quantitative_contingencies_into_strict_contingencies(rxnconsys.quantitative_contingencies_for_reaction(reaction))
    additional_contingency_state_set= state_set_from_contingencies(additional_strict_cont) # todo: ersetze k+ durch ! und k- durch x. in venntastic: ! ^= PropertySet(), x ^= Complement(PropertySet())
    vennset = venn.Intersection(vennset, additional_contingency_state_set)
    return bbm.Rule(bbm.Node(reaction), bbm.Factor(vennset_to_bbm_factor_vennset(vennset.simplified_form())))


def convert_quantitative_contingencies_into_strict_contingencies(contingencies: tg.List[con.Contingency]):
    converted_contingencies = []
    for contingency in contingencies:
        if contingency.type == con.ContingencyType.positive:
            converted_contingencies.append(con.Contingency(contingency.target, con.ContingencyType.requirement, contingency.effector))
        elif contingency.type == con.ContingencyType.negative:
            converted_contingencies.append(con.Contingency(contingency.target, con.ContingencyType.inhibition, contingency.effector))
        else:
            raise NotImplementedError

    return converted_contingencies


def vennset_to_bbm_factor_vennset(vennset: venn.Set):
    # creates new vennset with states contained by Node objects, for compareability
    # want to rewrite venn.Set into bbm.Factor like venn.PropertySet(A--B) -> venn.PropertySet(bbm.Node(A--B))
    if isinstance(vennset, venn.PropertySet):
        return venn.PropertySet(bbm.Node(vennset.value))
    if isinstance(vennset, venn.Complement):
        return venn.Complement(vennset_to_bbm_factor_vennset(vennset.expr))
    elif isinstance(vennset, venn.Intersection):
        return venn.Intersection(vennset_to_bbm_factor_vennset(vennset.left_expr), vennset_to_bbm_factor_vennset(vennset.right_expr))
    elif isinstance(vennset, venn.Union):
        return venn.Union(vennset_to_bbm_factor_vennset(vennset.left_expr), vennset_to_bbm_factor_vennset(vennset.right_expr))
    else:
        raise NotImplementedError


def get_rule_targets(rules: tg.List[bbm.Rule]):

    all_visited_states = [bbm.Node(rule.target.value) for rule in rules if rule]
    return all_visited_states


def rule_for_state_from_rxnconsys_and_reaction(rxnconsys: rxs.RxnConSystem, reaction: rxn.Reaction, system_rules: tg.List[bbm.Rule]) -> bbm.Rule:
    all_visited_nodes = get_rule_targets(system_rules)

    if reaction.product is None or bbm.Node(reaction.product) in all_visited_nodes:
        return None

    pos_bool_def=[venn.PropertySet(bbm.Node(reaction.product)), venn.PropertySet(bbm.Node(reaction))]
    neg_bool_def=[]

    for rxn in rxnconsys.reactions:
        if rxn.product is not None and rxn != reaction and reaction.product in [rxn.product]:
            pos_bool_def.append(venn.PropertySet(bbm.Node(rxn)))
        if rxn.source is not None and rxn != reaction and reaction.product in [rxn.source]:
            neg_bool_def.append(venn.PropertySet(bbm.Node(rxn)))

    pos_rules= venn.nested_expression_from_list_and_binary_op(pos_bool_def, venn.Union)
    neg_rules = venn.nested_expression_from_list_and_binary_op(neg_bool_def, venn.Union)
    vennset = venn.Intersection(pos_rules, venn.Complement(neg_rules))
    return bbm.Rule(bbm.Node(reaction.product),
                    bbm.Factor(vennset.simplified_form()))
