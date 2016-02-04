from typing import List
import typecheck as tc
import copy

import rxncon.core.contingency as con
import rxncon.venntastic.sets as venn
import rxncon.core.effector as eff
import rxncon.core.reaction as rxn


class StateFlow:
    def __init__(self, reaction: rxn.Reaction, source: venn.Set, target: venn.Set):
        if len(source.to_nested_list_form()) != 1:
            raise AssertionError('Number of OR-terms in StateFlow source needs to be 1, is {0}'
                                 .format(len(source.to_nested_list_form())))

        if len(target.to_nested_list_form()) != 1:
            raise AssertionError('Number of OR-terms in StateFlow target needs to be 1, is {0}'
                                 .format(len(source.to_nested_list_form())))

        self.reaction = reaction
        self.source = source
        self.target = target

    def __str__(self):
        return 'Source:{0}, Target:{1}'.format(self.source, self.target)


@tc.typecheck
def boolean_state_flows(reaction: rxn.Reaction,
                        strict_contingencies: List[con.Contingency],
                        source_contingencies: List[con.Contingency]) -> List[StateFlow]:

    for contingency in strict_contingencies + source_contingencies:
        assert contingency.type in [con.ContingencyType.requirement, con.ContingencyType.inhibition]

    or_terms = set_from_contingencies(strict_contingencies).to_union_list_form()

    return [StateFlow(reaction,
                      venn.Intersection(or_term, set_from_contingencies(source_contingencies)),
                      venn.Intersection(or_term, venn.Complement(set_from_contingencies(source_contingencies))))
            for or_term in or_terms]


@tc.typecheck
def quantified_state_flows(state_flow: StateFlow, quantitative_contingencies: List[con.Contingency]) -> List[StateFlow]:
    quantified_flows = [state_flow]

    for contingency in quantitative_contingencies:
        assert contingency.type in [con.ContingencyType.positive, con.ContingencyType.negative]

        new_flows = []
        for flow in quantified_flows:
            # Split the flow in (possibly more than one) WITH and (possibly more than one) WITHOUT the effector.
            # The reason this could be more than one is that the states making up the effector can be joined by Unions.
            for effector in set_from_effector(contingency.effector).to_union_list_form():
                new_flow = copy.deepcopy(flow)
                new_flow.source = venn.Intersection(new_flow.source, effector)
                new_flow.target = venn.Intersection(new_flow.target, effector)

                new_flows.append(new_flow)

            for effector in venn.Complement(set_from_effector(contingency.effector)).to_union_list_form():
                new_flow = copy.deepcopy(flow)
                new_flow.source = venn.Intersection(new_flow.source, effector)
                new_flow.target = venn.Intersection(new_flow.target, effector)

                new_flows.append(new_flow)

        quantified_flows = new_flows

    return quantified_flows


def disjunctified_state_flows(state_flow: List[StateFlow]) -> List[StateFlow]:
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