from collections import defaultdict
from enum import Enum
from typing import List, Set

import typecheck as tc

import rxncon.core.rxncon_system as rxs
import rxncon.core.state as sta
import rxncon.simulation.rule_based.rule_based_model as rbm
import rxncon.semantics.state_flow as flo
import rxncon.venntastic.sets as venn


class Default(Enum):
    unmodified = 'u'
    association_prefix = 'Assoc'
    modification_prefix = 'Mod'


#@tc.typecheck
def rule_based_model_from_rxncon(rxnconsys: rxs.RxnConSystem) -> rbm.RuleBasedModel:
    molecule_defs = molecule_defs_from_rxncon(rxnconsys)

    flows = []

    for reaction in rxnconsys.reactions:
        boolean_flow = flo.boolean_state_flows(reaction, rxnconsys.strict_contingencies_for_reaction(reaction),
                                                       rxnconsys.source_contingencies_for_reaction(reaction))

        derived_rule_conditions = []

        # for condition in base_rule_conditions:
        #     derived_rule_conditions += \
        #         expand_base_rule_condition(condition,
        #                                    rxnconsys.quantitative_contingencies_for_reaction(reaction))
        #
        # rule_conditions += derived_rule_conditions
        flows += boolean_flow
    # rules = []
    #
    # for disjunct_rule_condition in disjunct_rule_conditions_from_rule_conditions(rule_conditions):
    #     rules.append(rule_from_rule_condition_and_molecule_definitions(disjunct_rule_condition, molecule_defs))
    #
    # return rules
    rules = []
    for flow in flows:
        rules.append(rule_from_flow_and_molecule_definitions(flow, molecule_defs))

    return rules

#@tc.typecheck
def rule_from_flow_and_molecule_definitions(flow: flo.StateFlow, molecule_defs: List[rbm.MoleculeDefinition]) -> rbm.Rule:
    def __state_type_helper(source):
        if isinstance(source, venn.Complement):
            return source.expr.value
        elif isinstance(source, venn.PropertySet):
            return source.value

    name_to_association_specifications = defaultdict(set)
    name_to_modification_specifications = defaultdict(set)
    name_to_localisation_specification = defaultdict(lambda: None)

    names = set()

    assert len(flow.source.to_nested_list_form()) == 1

    for source in flow.source.to_nested_list_form()[0]:
        # todo: source.expr.value only for complement otherwise source.value
        if isinstance(__state_type_helper(source), sta.InterProteinInteractionState):
            name_to_association_specifications.update({key: (name_to_association_specifications[key].union(value)
                                                       if key in name_to_association_specifications else values)
                                                       for key, value in association_specification_from_set(source).items()})
            names.update(names_from_set(source))

    source_mol_specs = source_molecule_specifications(names, name_to_association_specifications,
                                   name_to_modification_specifications,
                                   name_to_localisation_specification,
                                   molecule_defs)


    # for target in flow.target.to_nested_list_form()[0]:
    #     if isinstance(target.expr.value, sta.InterProteinInteractionState):
    #         name_to_association_specifications = association_specification_from_set(name_to_association_specifications, target)


def source_reactant(source_mol_specs: defaultdict) -> List[rbm.Reactant]:
    reactants = []
    #rbm.MoleculeReactant(mol_spec_a_unbound)

def source_molecule_specifications(names: Set[str],
                                   name_to_association_specifications: defaultdict,
                                   name_to_modification_specifications: defaultdict,
                                   name_to_localisation_specification: defaultdict,
                                   molecule_defs: List[rbm.MoleculeDefinition]) -> defaultdict:

    molecule_specifications = defaultdict()
    for name in names:
        mol_def = [mol_def for mol_def in molecule_defs if mol_def.name == name]
        molecule_specifications[name] = rbm.MoleculeSpecification(mol_def[0],
                                                                  list(name_to_modification_specifications[name]),
                                                                  list(name_to_association_specifications[name]),
                                                                  name_to_localisation_specification[name])
    return molecule_specifications



def names_from_set(set_state):
    if isinstance(set_state, venn.Complement):
        assert isinstance(set_state.expr, venn.PropertySet)
        assert isinstance(set_state.expr.value, sta.State)
        return {set_state.expr.value.first_component.name, set_state.expr.value.second_component.name}
    elif isinstance(set_state, venn.PropertySet):
        assert isinstance(set_state.value, sta.State)
        return {set_state.value.first_component, set_state.value.second_component.name}
    else:
        raise NotImplementedError

def association_specification_from_set(set_state: venn.Set):
    name_to_association_specifications = defaultdict(set)
    if isinstance(set_state, venn.Complement):
        assert isinstance(set_state.expr, venn.PropertySet)
        assert isinstance(set_state.expr.value, sta.State)
        name_to_association_definition = name_to_assoc_defs_from_state(set_state.expr.value)
        first_comp_assoc_spec = rbm.AssociationSpecification(list(name_to_association_definition[set_state.expr.value.first_component.name])[0],
                                                             rbm.OccupationStatus.not_occupied)
        second_comp_assoc_spec =rbm.AssociationSpecification(list(name_to_association_definition[set_state.expr.value.second_component.name])[0],
                                                             rbm.OccupationStatus.not_occupied)
        name_to_association_specifications[set_state.expr.value.first_component.name].add(first_comp_assoc_spec)
        name_to_association_specifications[set_state.expr.value.second_component.name].add(second_comp_assoc_spec)

    elif isinstance(set_state, venn.PropertySet):
        assert isinstance(set_state.value, sta.State)
        name_to_association_definition = name_to_assoc_defs_from_state(set_state.expr.value)
        first_comp_assoc_spec = rbm.AssociationSpecification(list(name_to_association_definition[set_state.expr.value.first_component.name])[0],
                                                             rbm.OccupationStatus.not_occupied)
        second_comp_assoc_spec =rbm.AssociationSpecification(list(name_to_association_definition[set_state.expr.value.second_component.name])[0],
                                                             rbm.OccupationStatus.not_occupied)

    else:
        raise NotImplementedError
    return name_to_association_specifications


def modification_specification_from_set():
    pass

def localisation_specification_from_set():
    pass

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
def name_to_mod_defs_from_state_or_reaction(state: sta.CovalentModificationState, reaction: rxs.rxn.Reaction) -> dict:
    name_to_modification_definitions = defaultdict(set)
    modifiers = [Default.unmodified.value, state.modifier.value]
    domain = state.substrate.domain if state.substrate.domain else '{0}{1}'.format(Default.modification_prefix.value, reaction.subject.name)

    name_to_modification_definitions[state.substrate.name].add(rbm.ModificationDefinition(domain, modifiers))

    return name_to_modification_definitions


@tc.typecheck
def name_to_assoc_defs_from_state(state: tc.any(sta.InterProteinInteractionState, sta.IntraProteinInteractionState)):
    name_to_association_definitions = defaultdict(set)
    first_domain = state.first_component.domain if state.first_component.domain \
        else '{0}{1}'.format(Default.association_prefix.value, state.second_component.name)

    second_domain = state.second_component.domain if state.second_component.domain \
        else '{0}{1}'.format(Default.association_prefix.value, state.first_component.name)

    name_to_association_definitions[state.first_component.name].add(rbm.AssociationDefinition(first_domain))
    name_to_association_definitions[state.second_component.name].add(rbm.AssociationDefinition(second_domain))

    return name_to_association_definitions


def name_to_defs_from_rxncon(rxnconsys: rxs.RxnConSystem):
    names = set()
    name_to_association_definitions = defaultdict(set)
    name_to_modification_definitions = defaultdict(set)
    name_to_valid_compartments = defaultdict(set)

    for reaction in rxnconsys.reactions:
        for state in [x for x in (reaction.source, reaction.product) if x]:
            if isinstance(state, sta.CovalentModificationState):
                names.update(reaction.subject.name, state.substrate.name)
                name_to_modification_definitions.update({key: (name_to_modification_definitions[key].union(value)
                                                       if key in name_to_modification_definitions else value)
                                                       for key, value in name_to_mod_defs_from_state_or_reaction(state, reaction).items()})

            elif isinstance(state, sta.InterProteinInteractionState) or isinstance(state, sta.IntraProteinInteractionState):
                names.update({state.first_component.name, state.second_component.name})

                name_to_association_definitions.update({key: (name_to_association_definitions[key].union(value)
                                                       if key in name_to_association_definitions else value)
                                                       for key, value in name_to_assoc_defs_from_state(state).items()})

            elif isinstance(state, sta.SynthesisDegradationState):
                names.update({reaction.subject.name, state.component.name})

            elif isinstance(state, sta.TranslocationState):
                names.update({reaction.subject.name, state.substrate.name})
                name_to_valid_compartments[state.substrate.name].add(state.compartment)
    return names, name_to_modification_definitions, name_to_association_definitions, name_to_valid_compartments
