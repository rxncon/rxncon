from collections import defaultdict
from enum import Enum
from typing import List, Set, Dict

import typecheck as tc

import rxncon.core.rxncon_system as rxs
import rxncon.core.reaction as rxn
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

def molecule_specification_dict_from_state_and_molecule_def(molecule_def: List[rbm.MoleculeDefinition], state: sta.State, reaction: rxs.rxn.Reaction,
                                     take_complement: bool) -> defaultdict:
    name_to_association_specifications = defaultdict(set)
    name_to_modification_specifications = defaultdict(set)
    name_to_localisation_specification = defaultdict(lambda: None)

    names = set()

    if isinstance(state, sta.InterProteinInteractionState) or isinstance(state, sta.IntraProteinInteractionState):
        names.update(state.first_component.name, state.second_component.name)
        if take_complement:
            name_to_association_definition = name_to_assoc_defs_from_state(state)
            name_to_association_specifications[state.first_component.name].add(rbm.AssociationSpecification(list(name_to_association_definition[state.first_component.name])[0],
                                                                               rbm.OccupationStatus.not_occupied))
            name_to_association_specifications[state.second_component.name].add(rbm.AssociationSpecification(list(name_to_association_definition[state.second_component.name])[0],
                                                                                rbm.OccupationStatus.not_occupied))

        elif not take_complement:
            name_to_association_definition = name_to_assoc_defs_from_state(state)
            # todo: unknown partner?
            name_to_association_specifications[state.first_component.name].add(rbm.AssociationSpecification(list(name_to_association_definition[state.first_component.name])[0],
                                                                               rbm.OccupationStatus.occupied_known_partner))
            name_to_association_specifications[state.second_component.name].add(rbm.AssociationSpecification(list(name_to_association_definition[state.first_component.name])[0],
                                                                                rbm.OccupationStatus.occupied_known_partner))
        else:
            raise NotImplementedError
    elif isinstance(state, sta.CovalentModificationState):
        names.update(reaction.subject.name, state.substrate.name)
        modification_def = name_to_mod_defs_from_state_or_reaction(state,reaction)
        if take_complement:
            name_to_modification_specifications[state.substrate.name].add(rbm.ModificationSpecification(modification_def, Default.unmodified.value))
        else:
            name_to_modification_specifications[state.substrate.name].add(rbm.ModificationSpecification(modification_def, state.modifier.value))

    elif isinstance(state, sta.TranslocationState):
        names.update({reaction.subject.name, state.substrate.name})
        name_to_localisation_specification[state.substrate.name].add(state.compartment)

    elif isinstance(state, sta.SynthesisDegradationState):
        names.update({reaction.subject.name, state.component.name})


    molecule_specifications = molecule_specification_from_name_to_defaultdict(names, name_to_association_specifications,
                                                    name_to_modification_specifications,
                                                    name_to_localisation_specification,
                                                    molecule_def)

    return molecule_specifications

def reactants_from_specs_pair(mol_defs: Dict[str, rbm.MoleculeDefinition], specs_pair):
    source_specs = specs_pair[0]
    target_specs = specs_pair[1]




def specs_pair_from_flow(mol_defs: Dict[str, rbm.MoleculeDefinition], flow: flo.StateFlow):
    source_specs = {}
    assert len(flow.source.to_nested_list_form()) == 1
    for setstate in flow.source.to_nested_list_form()[0]:
        if isinstance(setstate, venn.Complement):
            source_specs.update(specs_from_state(mol_defs, setstate.expr.value, True))
        elif isinstance(setstate, venn.PropertySet):
            source_specs.update(specs_from_state(mol_defs, setstate.value, False))
        else:
            raise AssertionError

    target_specs = {}
    assert len(flow.target.to_nested_list_form()) == 1
    for setstate in flow.target.to_nested_list_form()[0]:
        if isinstance(setstate, venn.Complement):
            target_specs.update(specs_from_state(mol_defs, setstate.expr.value, True))
        elif isinstance(setstate, venn.PropertySet):
            target_specs.update(specs_from_state(mol_defs, setstate.value, False))
        else:
            raise AssertionError

    return source_specs, target_specs


def specs_from_state(mol_defs: Dict[str, rbm.MoleculeDefinition], state: sta.State, take_complement: bool) -> Dict:
    if isinstance(state, sta.CovalentModificationState) and not take_complement:
        mod_def = mol_defs[state.substrate.name].modification_def_by_domain_name(state.substrate.domain)
        return {state.substrate.name: rbm.ModificationSpecification(mod_def, state.modifier.value)}

    elif isinstance(state, sta.CovalentModificationState) and take_complement:
        mod_def = mol_defs[state.substrate.name].modification_def_by_domain_name(state.substrate.domain)
        return {state.substrate.name: rbm.ModificationSpecification(mod_def, Default.unmodified.value)}

    elif isinstance(state, sta.InterProteinInteractionState) and not take_complement:
        return






def molecule_specification_from_name_to_defaultdict(names: Set[str],
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
    pass


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
