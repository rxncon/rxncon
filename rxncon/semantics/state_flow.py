from typing import List

import rxncon.core.contingency as con
import rxncon.venntastic.sets as venn
import rxncon.core.effector as eff
import rxncon.core.reaction as rxn


class StateFlow:
    def __init__(self, reaction: rxn.Reaction, source: venn.Set, target: venn.Set):
        self.source = source
        self.target = target

    def __str__(self):
        return 'Source:{0}, Target:{1}'.format(self.source, self.target)


def boolean_state_flows(reaction: rxn.Reaction,
                        strict_contingencies: List[con.Contingency],
                        source_contingencies: List[con.Contingency]) -> List[StateFlow]:

    or_terms = set_from_contingencies(strict_contingencies).to_union_list_form()

    return [StateFlow(reaction,
                      venn.Intersection(or_term, set_from_contingencies(source_contingencies)),
                      venn.Intersection(or_term, venn.Complement(set_from_contingencies(source_contingencies))))
            for or_term in or_terms]


def quantified_state_flows(transfer_map: StateFlow, quantitative_contingencies: List[con.Contingency]) -> List[StateFlow]:
    pass


def disjunctified_state_flows(rule_conditions: List[StateFlow]) -> List[StateFlow]:
    pass


def set_from_effector(effector: eff.Effector) -> venn.Set:
    if isinstance(effector, eff.StateEffector):
        return venn.PropertySet(effector.expr)

    elif isinstance(effector, eff.AndEffector):
        return venn.Intersection(set_from_effector(effector.left_expr), set_from_effector(effector.right_expr))

    elif isinstance(effector, eff.OrEffector):
        return venn.Union(set_from_effector(effector.left_expr), set_from_effector(effector.right_expr))

    elif isinstance(effector, eff.NotEffector):
        return venn.Complement(set_from_effector(effector.expr))

    else:
        raise NotImplementedError


def set_from_contingencies(contingencies: List[con.Contingency]) -> venn.Set:
    intersections = []
    for contingency in contingencies:
        if contingency.type == con.ContingencyType.requirement:
            intersections.append(set_from_effector(contingency.effector))

        elif contingency.type == con.ContingencyType.inhibition:
            intersections.append(venn.Complement(set_from_effector(contingency.effector)))

        else:
            raise NotImplementedError

    return venn.nested_expression_from_list_and_binary_op(intersections, venn.Intersection)