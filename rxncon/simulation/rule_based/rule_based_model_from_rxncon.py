from collections import defaultdict
from enum import Enum
from typing import List, Set, Dict

import typecheck as tc

import rxncon.core.rxncon_system as rxs
import rxncon.core.component as com
import rxncon.core.state as sta
import rxncon.simulation.rule_based.rule_based_model as rbm
import rxncon.semantics.state_flow as flo
import rxncon.venntastic.sets as venn


class Default(Enum):
    unmodified = 'u'
    association_prefix = 'Assoc'
    modification_prefix = 'Mod'

def rule_based_model_from_rxncon(rxnconsys: rxs.RxnConSystem) -> rbm.RuleBasedModel:
    molecule_defs = molecule_defs_from_rxncon(rxnconsys)

    flows = []

    for reaction in rxnconsys.reactions:
        boolean_flow = flo.boolean_state_flows(reaction, rxnconsys.strict_contingencies_for_reaction(reaction),
                                                       rxnconsys.source_contingencies_for_reaction(reaction))

        # derived_rule_conditions
    # disjunct rule conditions

        flows += boolean_flow

        # rules
    rules = []

    for rxn_number, flow in enumerate(flows):
        rules.append(rule_from_flow_and_molecule_definitions(molecule_defs, flow , flow.reaction, rxn_number))

    return rules


def rule_from_flow_and_molecule_definitions_and_reaction(mol_def: Dict[str, rbm.MoleculeDefinition], flow: flo.StateFlow, reaction: rxs.rxn.Reaction, rxn_number: int):
    specification_pairs = specs_pair_from_flow(mol_def, flow, reaction)
    reactant_pairs = reactants_from_specs_pair(mol_def, specification_pairs)

    if reaction.directionality.irreversible:
        rates = [rbm.Parameter("k"+str(rxn_number),"1")]
        rbm.Rule(reactant_pairs[0], reactant_pairs[1], rbm.Arrow.irreversible, rates)
    elif reaction.directionality.reversible:
        rates = [rbm.Parameter("kf"+str(rxn_number),"1"), rbm.Parameter("kr"+str(rxn_number),"1")]
        rbm.Rule(reactant_pairs[0], reactant_pairs[1], rbm.Arrow.irreversible, rates)


def molecule_specification_from_spec(mol_defs: Dict[str, rbm.MoleculeDefinition], spec):
    molecule_specifications = []
    for mol, specs in spec.items():
        modification_specifcations = []
        association_specifications = []
        localisation_specification = None
        for spec in specs:
            if isinstance(spec, rbm.ModificationSpecification):
                modification_specifcations.append(spec)
            elif isinstance(spec, rbm.AssociationSpecification):
                association_specifications.append(spec)
            elif isinstance(spec, rbm.LocalizationSpecification):
                localisation_specification = spec
        molecule_specifications.append(rbm.MoleculeSpecification(mol_defs[mol], modification_specifcations,
                                                          association_specifications, localisation_specification))
    return molecule_specifications


def reactants_from_specs_pair(mol_defs: Dict[str, rbm.MoleculeDefinition], specs_pair):
    source_specs = specs_pair[0]
    target_specs = specs_pair[1]

    source_molecule_specifications = molecule_specification_from_spec(mol_defs, source_specs)
    source_reactants = []
    for mol_specification in source_molecule_specifications:
        source_reactants.append(rbm.MoleculeReactant(mol_specification))

    target_molecule_specifications = molecule_specification_from_spec(mol_defs, target_specs)
    target_reactants = []
    for mol_specification in target_molecule_specifications:
        target_reactants.append(rbm.MoleculeReactant(mol_specification))

    return source_reactants, target_reactants


def specs_pair_from_flow(mol_defs: Dict[str, rbm.MoleculeDefinition], flow: flo.StateFlow, reaction: rxs.rxn.Reaction):
    source_specs = defaultdict(set)
    assert len(flow.source.to_nested_list_form()) == 1
    for setstate in flow.source.to_nested_list_form()[0]:
        if isinstance(setstate, venn.Complement):
            for mol, spec in specs_from_state(mol_defs, setstate.expr.value, reaction, True).items():
                if spec:
                    source_specs[mol].add(spec)
                elif mol not in source_specs:
                    source_specs[mol]
        elif isinstance(setstate, venn.PropertySet):
            for mol, spec in specs_from_state(mol_defs, setstate.value, reaction, False).items():
                if spec:
                    source_specs[mol].add(spec)
                elif mol not in source_specs:
                    source_specs[mol]
        else:
            raise AssertionError

    target_specs = defaultdict(set)
    assert len(flow.target.to_nested_list_form()) == 1
    for setstate in flow.target.to_nested_list_form()[0]:
        if isinstance(setstate, venn.Complement):
            for mol, spec in specs_from_state(mol_defs, setstate.expr.value, reaction, True).items():
                if spec:
                    target_specs[mol].add(spec)
                elif mol not in target_specs:
                    target_specs[mol]
        elif isinstance(setstate, venn.PropertySet):
            for mol, spec in specs_from_state(mol_defs, setstate.value, reaction, False).items():
                if spec:
                    target_specs[mol].add(spec)
                elif mol not in target_specs:
                    target_specs[mol]
        else:
            raise AssertionError

    return source_specs, target_specs

 # todo: change default
def specs_from_state(mol_defs: Dict[str, rbm.MoleculeDefinition], state: sta.State, reaction: rxs.rxn.Reaction, take_complement: bool) -> Dict:
    # MODIFICATION
    if isinstance(state, sta.CovalentModificationState):
        return modification_specification_from_state(mol_defs, state, reaction, take_complement)
    # INTERACTION
    elif isinstance(state, sta.InterProteinInteractionState) or isinstance(state, sta.IntraProteinInteractionState):
        return association_specification_from_state(mol_defs, state, take_complement)

    # TRANSLOCATION
    elif isinstance(state, sta.TranslocationState) and not take_complement:
        # todo: what is the not of a TranslocationState
        #return {state.substrate.name: rbm.LocalizationSpecification(mol_defs[state.substrate.name].localization_def, state.compartment)}
        assert NotImplementedError
    elif isinstance(state, sta.TranslocationState) and take_complement:
        assert NotImplementedError

    # SYNTHESISDEGRADATION
    elif isinstance(state, sta.SynthesisDegradationState) and not take_complement:
        assert NotImplementedError
    elif isinstance(state, sta.SynthesisDegradationState) and take_complement:
        assert NotImplementedError
    else:
        assert AssertionError

def association_domain_name(reference_component: com.Component, alternativ_component: com.Component):
    return reference_component.domain if reference_component.domain \
        else '{0}{1}'.format(Default.association_prefix.value, alternativ_component.name)

def modification_domain_name(substrate: com.Component, reaction: rxs.rxn.Reaction):
    return substrate.domain if substrate.domain else '{0}{1}'.format(Default.modification_prefix.value, reaction.subject.name)


def modification_specification_from_state(mol_defs: Dict[str, rbm.MoleculeDefinition] ,state: venn.Set, reaction: rxs.rxn.Reaction, take_complement: bool):
    if not take_complement:
        mod_def = mol_defs[state.substrate.name].modification_def_by_domain_name(modification_domain_name(state.substrate, reaction))
        return {state.substrate.name: rbm.ModificationSpecification(mod_def, state.modifier.value),
                reaction.subject.name: None}

    elif take_complement:
        mod_def = mol_defs[state.substrate.name].modification_def_by_domain_name(modification_domain_name(state.substrate, reaction))
        return {state.substrate.name: rbm.ModificationSpecification(mod_def, Default.unmodified.value),
                reaction.subject.name: None}


def association_specification_from_state(mol_defs: Dict[str, rbm.MoleculeDefinition] ,state: venn.Set, take_complement: bool):
    if not take_complement:
            # todo: unknown partner?
        return {state.first_component.name: rbm.AssociationSpecification(mol_defs[state.first_component.name].association_def_by_domain_name(association_domain_name(state.first_component, state.second_component)),
                                                                         rbm.OccupationStatus.occupied_known_partner),
                state.second_component.name: rbm.AssociationSpecification(mol_defs[state.second_component.name].association_def_by_domain_name(association_domain_name(state.second_component, state.first_component)),
                                                                          rbm.OccupationStatus.occupied_known_partner)}
    elif take_complement:
        return {state.first_component.name: rbm.AssociationSpecification(mol_defs[state.first_component.name].association_def_by_domain_name(association_domain_name(state.first_component, state.second_component)),
                                                                         rbm.OccupationStatus.not_occupied),
                state.second_component.name: rbm.AssociationSpecification(mol_defs[state.second_component.name].association_def_by_domain_name(association_domain_name(state.second_component, state.first_component)),
                                                                          rbm.OccupationStatus.not_occupied)}


def modification_specification_from_set():
    pass

def localisation_specification_from_set():
    pass

@tc.typecheck
def molecule_defs_from_rxncon(rxnconsys: rxs.RxnConSystem) -> Dict[str, rbm.MoleculeDefinition]:
    names, name_to_mod_defs, name_to_assoc_defs, name_to_valid_compartments = name_to_defs_from_rxncon(rxnconsys)

    molecule_definitions = dict()

    for name in names:
        association_definitions = list(name_to_assoc_defs[name]) if name in name_to_assoc_defs else None
        modification_definitions = list(name_to_mod_defs[name]) if name in name_to_mod_defs else None
        localisation_definition = rbm.LocalizationDefinition(list(name_to_valid_compartments[name])) \
            if name in name_to_valid_compartments else None
        rbm.MoleculeDefinition(name, modification_definitions,association_definitions, localisation_definition)
        molecule_definitions[name] = rbm.MoleculeDefinition(name, modification_definitions, association_definitions, localisation_definition)

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
