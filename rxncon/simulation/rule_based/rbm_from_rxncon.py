from collections import defaultdict
import itertools as itt
import functools as ft
import typing as tg

import rxncon.core.rxncon_system as rxs
import rxncon.simulation.rule_based.rule_based_model as rbm
import rxncon.core.state as sta
import rxncon.core.reaction as rxn
import rxncon.core.component as com
import rxncon.core.contingency as con
import rxncon.core.effector as eff
import rxncon.venntastic.sets as venn


def rules_from_rxncon(rxconsys: rxs.RxnConSystem):
    mol_defs = molecule_defs_from_rxncon(rxconsys)

    rules = []
    for reaction in rxconsys.reactions:
        rules += rules_from_reaction(rxconsys, mol_defs, reaction)

    return rules


def rules_from_reaction(rxnconsys: rxs.RxnConSystem, mol_defs: tg.Dict[str, rbm.MoleculeDefinition],
                        reaction: rxn.Reaction) -> tg.List[rbm.Rule]:
    relevant_molecule_names = set(component.name for component in rxnconsys.components_for_reaction(reaction))

    strict_contingency_state_set = state_set_from_contingencies(rxnconsys.strict_contingencies_for_reaction(reaction))
    source_state_set = source_state_set_from_reaction(reaction)

    quant_contingency_configurations = \
        contingency_configurations_from_quantitative_contingencies(rxnconsys.quantitative_contingencies_for_reaction(reaction))

    mol_specs_lhs_to_rhs = defaultdict(list)

    for molecule_name in relevant_molecule_names:
        mol_def = mol_defs[molecule_name]

        strict_spec_set = mol_def.specification_set_from_state_set(strict_contingency_state_set)
        source_spec_set = mol_def.specification_set_from_state_set(source_state_set)

        lhs_sets = []

        for strict_term in strict_spec_set.to_union_list_form():
            for contingency in quant_contingency_configurations:
                quant_contingency_state_set = state_set_from_contingencies([contingency])
                quant_term = mol_def.specification_set_from_state_set(quant_contingency_state_set)

                lhs_sets.append(venn.Intersection(strict_term, quant_term))

            else:
                lhs_sets.append(strict_term)

        lhs_sets_disjunct = venn.gram_schmidt_disjunctify(lhs_sets)

        if lhs_sets_disjunct == [venn.EmptySet()]:
            continue

        for lhs in lhs_sets_disjunct:
            if source_spec_set == venn.EmptySet():
                left_mol_spec = mol_def.specification_from_specification_set(lhs)
                right_mol_spec = left_mol_spec

            else:
                left_mol_spec = mol_def.specification_from_specification_set(venn.Intersection(lhs, source_spec_set))
                right_mol_spec = mol_def.specification_from_specification_set(venn.Intersection(lhs, venn.Complement(source_spec_set)))

            mol_specs_lhs_to_rhs[left_mol_spec.molecule_def.name].append((left_mol_spec, right_mol_spec))

    rules = []

    # @todo
    arrow = rbm.Arrow.reversible
    rates = [rbm.Parameter('k1', '1.0'), rbm.Parameter('k2', '1.0')]

    for lhs_to_rhs_per_molecule in itt.product(*mol_specs_lhs_to_rhs.values()):
        lhs_to_rhs = list(zip(*lhs_to_rhs_per_molecule))
        lhs = reactants_from_specs(lhs_to_rhs[0])
        rhs = reactants_from_specs(lhs_to_rhs[1])

        rules.append(rbm.Rule(lhs, rhs, arrow, rates))

    return rules


def state_set_from_contingencies(contingencies: tg.List[con.Contingency]) -> venn.Set:
    if not contingencies:
        return venn.UniversalSet()

    for contingency in contingencies:
        assert contingency.target == contingencies[0].target
        assert contingency.type in [con.ContingencyType.inhibition, con.ContingencyType.requirement]

    requirements = []
    inhibitions = []

    for contingency in contingencies:
        if contingency.type == con.ContingencyType.requirement:
            requirements.append(state_set_from_effector(contingency.effector))

        elif contingency.type == con.ContingencyType.inhibition:
            inhibitions.append(state_set_from_effector(contingency.effector))

    required = venn.nested_expression_from_list_and_binary_op(requirements, venn.Intersection)
    inhibited = venn.nested_expression_from_list_and_binary_op(inhibitions, venn.Union)

    if required != venn.EmptySet() and inhibited != venn.EmptySet():
        return venn.Intersection(required, venn.Complement(inhibited))

    elif required == venn.EmptySet():
        return venn.Complement(inhibited)

    elif inhibited == venn.EmptySet():
        return required


def contingency_configurations_from_quantitative_contingencies(contingencies: tg.List[con.Contingency]) -> tg.List[con.Contingency]:
    if not contingencies:
        return []

    reaction = contingencies[0].target

    for contingency in contingencies:
        assert contingency.target == reaction
        assert contingency.type in [con.ContingencyType.positive, con.ContingencyType.negative]

    config_pairs = [(effector, eff.NotEffector(effector)) for effector in [contingency.effector for contingency in contingencies]]

    configs = itt.product(*config_pairs)

    contingency_configs = []

    for config in configs:
        if len(config) == 1:
            total_effector = config[0]
        else:
            total_effector = ft.reduce(eff.AndEffector, config[1:], config[0])

        contingency_configs.append(con.Contingency(reaction, con.ContingencyType.requirement, total_effector))

    return contingency_configs


def source_state_set_from_reaction(reaction: rxn.Reaction) -> venn.Set:
    source_state = reaction.source
    product_state = reaction.product

    if not source_state and product_state:
        return venn.Complement(venn.PropertySet(product_state))

    elif source_state and not product_state:
        return venn.PropertySet(source_state)

    elif source_state and product_state:
        return venn.Intersection(venn.Complement(venn.PropertySet(product_state)), venn.PropertySet(source_state))

    else:
        raise AssertionError


def reactants_from_specs(mol_specs: tg.List[rbm.MoleculeSpecification]) -> tg.List[rbm.Reactant]:
    reactants = []
    binding_spec = []
    binding = []
    for mol_spec in mol_specs:
        for assoc_spec in mol_spec.association_specs:
            if assoc_spec.occupation_status in [rbm.OccupationStatus.occupied_known_partner, rbm.OccupationStatus.occupied_unknown_partner]:
                binding_state = assoc_spec.association_def.matching_state
                for possible_assoc_mol_spec in mol_specs:
                    if possible_assoc_mol_spec != mol_spec:
                        for possible_binding in possible_assoc_mol_spec.association_specs:
                            if possible_binding.association_def.matching_state == binding_state:
                                if mol_spec not in binding_spec:
                                    binding_spec.append(mol_spec)
                                if possible_assoc_mol_spec not in binding_spec:
                                    binding_spec.append(possible_assoc_mol_spec)

                                binding_tuple = rbm.Binding((binding_spec.index(mol_spec),assoc_spec),
                                                       (binding_spec.index(possible_assoc_mol_spec),possible_binding))
                                check_tuple = rbm.Binding((binding_spec.index(possible_assoc_mol_spec),possible_binding),
                                                          (binding_spec.index(mol_spec),assoc_spec))
                                if  check_tuple not in binding:
                                    binding.append(binding_tuple)
                                else:
                                    binding_spec = binding_spec[:-2]


    if binding_spec:
        reactants.append((rbm.ComplexReactant(binding_spec, binding)))

    for mol_spec in mol_specs:
        if mol_spec not in binding_spec:
            reactants.append(rbm.MoleculeReactant(mol_spec))

    return reactants


def state_set_from_effector(effector: eff.Effector) -> venn.Set:
    if isinstance(effector, eff.StateEffector):
        return venn.PropertySet(effector.expr)

    elif isinstance(effector, eff.NotEffector):
        return venn.Complement(state_set_from_effector(effector.expr))

    elif isinstance(effector, eff.AndEffector):
        return venn.Intersection(state_set_from_effector(effector.left_expr), state_set_from_effector(effector.right_expr))

    elif isinstance(effector, eff.OrEffector):
        return venn.Union(state_set_from_effector(effector.left_expr), state_set_from_effector(effector.right_expr))

    else:
        raise AssertionError


def molecule_defs_from_rxncon(rxnconsys: rxs.RxnConSystem) -> tg.Dict[str, rbm.MoleculeDefinition]:
    names = set()

    name_to_ass_defs = defaultdict(list)
    name_to_locs = defaultdict(list)
    name_to_mod_defs = defaultdict(list)

    for reaction in rxnconsys.reactions:
        for state in [x for x in [reaction.source, reaction.product] if x]:
            if isinstance(state, sta.CovalentModificationState):
                names.add(state.substrate.name)
                names.add(reaction.subject.name)
                name_to_mod_defs[state.substrate.name].append(_mod_def_from_state_and_reaction(state, reaction))

            elif isinstance(state, sta.InterProteinInteractionState) or isinstance(state, sta.IntraProteinInteractionState):
                names.add(state.first_component.name)
                names.add(state.second_component.name)
                ass_defs = _state_to_ass_defs(state)
                name_to_ass_defs[state.first_component.name].append(ass_defs[0])
                name_to_ass_defs[state.second_component.name].append(ass_defs[1])

            elif isinstance(state, sta.TranslocationState):
                names.add(state.substrate.name)
                names.add(reaction.subject.name)
                name_to_locs[state.substrate.name].append(state.compartment)

            else:
                raise NotImplementedError

    mol_defs = {}

    for name in names:
        mol_defs[name] = rbm.MoleculeDefinition(name,
                                                name_to_mod_defs[name],
                                                name_to_ass_defs[name],
                                                rbm.LocalizationDefinition(name_to_locs[name]))

    return mol_defs


def _mod_def_from_state_and_reaction(state: sta.CovalentModificationState, reaction: rxn.Reaction):
    if not state.substrate.domain and not state.substrate.subdomain and not state.substrate.residue:
        domain_name = 'Mod{}'.format(reaction.subject.name)

    else:
        domain_name = state.substrate.domain

    mod_def = rbm.ModificationDefinition(domain_name, [sta.StateModifier.unmodified.value, state.modifier.value])
    mod_def.matching_state = state

    return mod_def


def _state_to_ass_defs(state: tg.Union[sta.InterProteinInteractionState, sta.IntraProteinInteractionState]):
    first_dom = state.first_component.domain if state.first_component.domain else 'Ass{}'.format(state.second_component.name)
    second_dom = state.second_component.domain if state.second_component.domain else 'Ass{}'.format(state.first_component.name)

    first_def = rbm.AssociationDefinition(first_dom)
    first_def.matching_state = state

    second_def = rbm.AssociationDefinition(second_dom)
    second_def.matching_state = state

    return first_def, second_def




