import typing as tg
import rxncon.core.rxncon_system as rxs
import rxncon.core.reaction as rxn
import rxncon.core.contingency as con
import rxncon.simulation.bBM.bipartite_boolean_model as bbm
import rxncon.venntastic.sets as venn
import rxncon.core.effector as eff
import rxncon.core.state as sta
import rxncon.core.specification as spec


def bipartite_boolean_model_from_rxncon(rxconsys: rxs.RxnConSystem):
    return bbm.BipartiteBooleanModel(rules_from_rxncon(rxconsys), initial_states_from_rxncon(rxconsys))


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
        # todo: change this later to a specific state
        if bbm.InitCondition(bbm.Node(sta.ComponentState(spec.Specification(reaction.subject.name,
                                                                            None, None, None))), None) not in initial_states:
            initial_states.append(bbm.InitCondition(bbm.Node(sta.ComponentState(spec.Specification(reaction.subject.name,
                                                                                                   None, None, None))), None))
        if bbm.InitCondition(bbm.Node(sta.ComponentState(spec.Specification(reaction.object.name,
                                                                            None, None, None))), None) not in initial_states:
            initial_states.append(bbm.InitCondition(bbm.Node(sta.ComponentState(spec.Specification(reaction.object.name,
                                                                                                   None, None, None))), None))
    return initial_states


def rule_for_reaction_from_rxnconsys_and_reaction(rxnconsys: rxs.RxnConSystem, reaction: rxn.Reaction,
                                                  system_rules: tg.List[bbm.Rule]) -> bbm.Rule:
    all_visited_nodes = get_rule_targets(system_rules)
    if bbm.Node(reaction) in all_visited_nodes:
        return None
    strict_contingency_state_set = _state_set_from_contingencies(rxnconsys.strict_contingencies_for_reaction(reaction))
    if isinstance(strict_contingency_state_set.to_full_simplified_form(), venn.EmptySet):
        raise AssertionError("There is no way to fulfill the contingencies: {}".format(strict_contingency_state_set))
    vennset = venn.Intersection(strict_contingency_state_set.to_full_simplified_form(),
                                # todo: change this later to a specific state
                                   venn.Intersection(venn.PropertySet(sta.ComponentState(spec.Specification(reaction.subject.name,
                                                                                                            None, None, None))),
                                                     venn.PropertySet(sta.ComponentState(spec.Specification(reaction.object.name,
                                                                                                            None, None, None)))))
    additional_strict_cont = convert_quantitative_contingencies_into_strict_contingencies(rxnconsys.quantitative_contingencies_for_reaction(reaction))
    additional_contingency_state_set = _state_set_from_contingencies(additional_strict_cont)

    if isinstance(additional_contingency_state_set.to_full_simplified_form(), venn.EmptySet):
        raise AssertionError("There is no way to fulfill the contingencies: {}".format(additional_contingency_state_set))

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
    if isinstance(vennset, venn.EmptySet):
        raise AssertionError

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
        if rxn.product is not None and rxn != reaction and reaction.product == rxn.product:
            pos_bool_def.append(venn.PropertySet(bbm.Node(rxn)))
        if rxn.source is not None and rxn != reaction and reaction.product == rxn.source:
            neg_bool_def.append(venn.PropertySet(bbm.Node(rxn)))

    pos_rules= venn.nested_expression_from_list_and_binary_op(pos_bool_def, venn.Union)
    neg_rules = venn.nested_expression_from_list_and_binary_op(neg_bool_def, venn.Union)
    vennset = venn.Intersection(pos_rules, venn.Complement(neg_rules))

    if isinstance(vennset.to_full_simplified_form(), venn.EmptySet):
        raise AssertionError("There is no way to fulfill this rule: {}".format(vennset))

    return bbm.Rule(bbm.Node(reaction.product),
                    bbm.Factor(vennset.to_full_simplified_form()))


def _state_set_from_contingencies(contingencies: tg.List[con.Contingency]) -> venn.Set:
    if not contingencies:
        return venn.UniversalSet()

    for contingency in contingencies:
        assert contingency.target == contingencies[0].target
        assert contingency.type in [con.ContingencyType.inhibition, con.ContingencyType.requirement]

    requirements = []
    inhibitions = []

    for contingency in contingencies:
        if contingency.type == con.ContingencyType.requirement:
            requirements.append(_state_set_from_effector(contingency.effector))

        elif contingency.type == con.ContingencyType.inhibition:
            inhibitions.append(_state_set_from_effector(contingency.effector))

    required_set = venn.nested_expression_from_list_and_binary_op(requirements, venn.Intersection)
    inhibited_set = venn.Complement(venn.nested_expression_from_list_and_binary_op(inhibitions, venn.Union))

    if requirements and inhibitions:
        return venn.Intersection(required_set, inhibited_set)

    elif inhibitions:
        return inhibited_set

    elif requirements:
        return required_set


def _state_set_from_effector(effector: eff.Effector) -> venn.Set:
    if isinstance(effector, eff.StateEffector):
        return venn.PropertySet(effector.expr)

    elif isinstance(effector, eff.NotEffector):
        return venn.Complement(_state_set_from_effector(effector.expr))

    elif isinstance(effector, eff.AndEffector):
        return venn.Intersection(_state_set_from_effector(effector.left_expr), _state_set_from_effector(effector.right_expr))

    elif isinstance(effector, eff.OrEffector):
        return venn.Union(_state_set_from_effector(effector.left_expr), _state_set_from_effector(effector.right_expr))

    else:
        raise AssertionError