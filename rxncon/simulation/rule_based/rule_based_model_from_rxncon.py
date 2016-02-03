from enum import Enum
from typing import List, Tuple
from collections import defaultdict
import typecheck as tc

import rxncon.core.rxncon_system as rxs
import rxncon.simulation.rule_based.rule_based_model as rbm
import rxncon.core.state as sta
import rxncon.core.reaction as rxn
import rxncon.core.contingency as con
import rxncon.core.effector as eff
import rxncon.venntastic.sets as venn


class Default(Enum):
    unmodified = 'u'
    association_prefix = 'Assoc'
    modification_prefix = 'Mod'


class RuleConditions:
    def __init__(self, lhs_conditions: venn.Set, rhs_conditions: venn.Set):
        self.lhs_conditions = lhs_conditions
        self.rhs_conditions = rhs_conditions

    def __str__(self):
        return 'LHS:{0}, RHS:{1}'.format(self.lhs_conditions, self.rhs_conditions)


@tc.typecheck
def rule_based_model_from_rxncon(rxnconsys: rxs.RxnConSystem) -> rbm.RuleBasedModel:
    molecule_defs = molecule_defs_from_rxncon(rxnconsys)

    rule_conditions = []

    for reaction in rxnconsys.reactions:
        base_rule_conditions = determine_base_rule_conditions(rxnconsys.strict_contingencies_for_reaction(reaction),
                                                              rxnconsys.source_contingencies_for_reaction(reaction))

        derived_rule_conditions = []

        for condition in base_rule_conditions:
            derived_rule_conditions += \
                expand_base_rule_condition(condition,
                                           rxnconsys.quantitative_contingencies_for_reaction(reaction))

        rule_conditions += derived_rule_conditions

    rules = []

    for disjunct_rule_condition in disjunct_rule_conditions_from_rule_conditions(rule_conditions):
        rules.append(rule_from_rule_condition_and_molecule_definitions(disjunct_rule_condition, molecule_defs))

    return rules


def determine_base_rule_conditions(strict_contingencies: List[con.Contingency],
                                   source_contingencies: List[con.Contingency]) -> List[RuleConditions]:

    set_of_requirements = set_from_contingencies(strict_contingencies)
    set_of_source_cont = set_from_contingencies(source_contingencies)

    or_terms = set_of_requirements.to_union_list_form()

    return [RuleConditions(venn.Intersection(or_term, set_of_source_cont),
                           venn.Intersection(or_term, venn.Complement(set_of_source_cont))) for or_term in or_terms]


def expand_base_rule_condition(rule_condition: RuleConditions, quantitative_contingencies: List[con.Contingency]) -> List[RuleConditions]:
    pass


def disjunct_rule_conditions_from_rule_conditions(rule_conditions: List[RuleConditions]) -> List[RuleConditions]:
    pass


def rule_from_rule_condition_and_molecule_definitions(rule_condition: RuleConditions, mol_defs: List[rbm.MoleculeDefinition]) -> rbm.Rule:
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


@tc.typecheck
def molecule_defs_from_rxncon(rxnconsys: rxs.RxnConSystem) -> List[rbm.MoleculeDefinition]:
    names, name_to_mod_defs, name_to_assoc_defs, name_to_valid_compartments = name_to_defs_from_rxncon(rxnconsys)

    molecule_definitions = []

    for name in names:
        association_definitions = list(name_to_assoc_defs[name]) if name in name_to_assoc_defs else None
        modification_definitions = list(name_to_mod_defs[name]) if name in name_to_mod_defs else None
        localisation_definition = rbm.LocalizationDefinition(list(name_to_valid_compartments[name])) if name in name_to_valid_compartments else None
        rbm.MoleculeDefinition(name, modification_definitions,association_definitions, localisation_definition)
        molecule_definitions.append(rbm.MoleculeDefinition(name, modification_definitions, association_definitions, localisation_definition))

    return molecule_definitions


@tc.typecheck
def name_to_mod_defs_from_state_or_reaction(name_to_modification_definitions: dict, state: sta.CovalentModificationState, reaction: rxs.rxn.Reaction) -> dict:
    modifiers = [Default.unmodified.value, state.modifier.value]
    domain = state.substrate.domain if state.substrate.domain else '{0}{1}'.format(Default.modification_prefix.value, reaction.subject.name)

    name_to_modification_definitions[state.substrate.name].add(rbm.ModificationDefinition(domain, modifiers))

    return name_to_modification_definitions


@tc.typecheck
def name_to_assoc_defs_from_state(name_to_association_definitions: dict, state: tc.any(sta.InterProteinInteractionState, sta.IntraProteinInteractionState)):

    first_domain = state.first_component.domain if state.first_component.domain \
        else '{0}{1}'.format(Default.association_prefix.value, state.second_component.name)

    second_domain = state.second_component.domain if state.second_component.domain \
        else '{0}{1}'.format(Default.association_prefix.value, state.first_component.name)

    name_to_association_definitions[state.first_component.name].add(rbm.AssociationDefinition(first_domain))
    name_to_association_definitions[state.second_component.name].add(rbm.AssociationDefinition(second_domain))

    return name_to_association_definitions


@tc.typecheck
def name_to_names_from_reaction(names: set, first_name: str, second_name: str) -> set:

    names.add(first_name)
    names.add(second_name)
    return names


def name_to_defs_from_rxncon(rxnconsys: rxs.RxnConSystem):
    names = set()
    name_to_association_definitions = defaultdict(set)
    name_to_modification_definitions = defaultdict(set)
    name_to_valid_compartments = defaultdict(set)

    for reaction in rxnconsys.reactions:
        for state in [x for x in (reaction.source, reaction.product) if x]:
            if isinstance(state, sta.CovalentModificationState):
                names = name_to_names_from_reaction(names, reaction.subject.name, state.substrate.name)
                name_to_modification_definitions = name_to_mod_defs_from_state_or_reaction(name_to_modification_definitions,
                                                                                           state, reaction)

            elif isinstance(state, sta.InterProteinInteractionState) or isinstance(state, sta.IntraProteinInteractionState):
                names = name_to_names_from_reaction(names, state.first_component.name, state.second_component.name)
                name_to_association_definitions = name_to_assoc_defs_from_state(name_to_association_definitions,
                                                                                state)

            elif isinstance(state, sta.SynthesisDegradationState):
                names = name_to_names_from_reaction(names, reaction.subject.name, state.component.name)

            elif isinstance(state, sta.TranslocationState):
                names = name_to_names_from_reaction(names, reaction.subject.name, state.substrate.name)
                name_to_valid_compartments[state.substrate.name].add(state.compartment)
    return names, name_to_modification_definitions, name_to_association_definitions, name_to_valid_compartments















