import itertools as itt
import typing as tg
from collections import defaultdict

import rxncon.core.reaction as rxn
import rxncon.core.rxncon_system as rxs
import rxncon.core.state as sta
import rxncon.semantics.molecule
import rxncon.simulation.rule_based.rule_based_model as rbm
import rxncon.venntastic.sets as venn
from rxncon.semantics.molecule_from_rxncon import set_of_states_from_contingencies, \
    contingency_configurations_from_quantitative_contingencies, source_set_of_states_from_reaction, \
    molecule_defs_from_rxncon


def rules_from_rxncon(rxconsys: rxs.RxnConSystem):
    mol_defs = molecule_defs_from_rxncon(rxconsys)

    rules = []
    for reaction in rxconsys.reactions:
        rules += rules_from_reaction(rxconsys, mol_defs, reaction)

    return rules


def rules_from_reaction(rxnconsys: rxs.RxnConSystem, mol_defs: tg.Dict[str, rxncon.semantics.molecule.MoleculeDefinition],
                        reaction: rxn.Reaction) -> tg.List[rbm.Rule]:
    relevant_molecule_names = set(component.name for component in rxnconsys.components_for_reaction(reaction))

    strict_contingency_state_set = set_of_states_from_contingencies(rxnconsys.strict_contingencies_for_reaction(reaction))
    source_state_set = source_set_of_states_from_reaction(reaction)

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
                quant_contingency_state_set = set_of_states_from_contingencies([contingency])
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
        lhs = reactants_from_specs(list(lhs_to_rhs[0]))
        rhs = reactants_from_specs(list(lhs_to_rhs[1]))

        rules.append(rbm.Rule(lhs, rhs, arrow, rates))

    return rules


def reactants_from_specs(molecule_specifications: tg.List[rxncon.semantics.molecule.MoleculeInstance]) -> tg.List[rbm.Reactant]:

    reactants = []

    specs_in_complexes = defaultdict(list)

    while molecule_specifications:
        mol_spec = molecule_specifications.pop()

        if not mol_spec.occupied_association_specs():
            reactants.append(rbm.MoleculeReactant(mol_spec))
        else:
            for mol_assoc_spec in mol_spec.occupied_association_specs():
                specs_in_complexes[mol_assoc_spec.association_def.matching_state].append((mol_spec, mol_assoc_spec))
    states = [state for state in specs_in_complexes.keys()]

    states_in_complex = complexes_by_state_connection(states)


    for states in states_in_complex:
        molecule_specifications = []
        bindings = []
        for state in states:
            if specs_in_complexes[state][0][0] not in molecule_specifications:
                molecule_specifications.append(specs_in_complexes[state][0][0])
            if specs_in_complexes[state][1][0] not in molecule_specifications:
                molecule_specifications.append(specs_in_complexes[state][1][0])
            binding = rxncon.semantics.molecule.Binding((molecule_specifications.index(specs_in_complexes[state][0][0]), specs_in_complexes[state][0][1]),
                                                        (molecule_specifications.index(specs_in_complexes[state][1][0]), specs_in_complexes[state][1][1]))
            bindings.append(binding)
        reactants.append(rbm.ComplexReactant(molecule_specifications, bindings))
    return reactants


def complexes_by_state_connection(states: tg.List[sta.State]):
    states_in_complex = []
    while states:
        state = states.pop()
        states, connected_states = _find_connectivity_of_state(state, states, connected_states)
        states_in_complex.append(connected_states)


    return states_in_complex


def _find_connectivity_of_state(state: sta.State, states: tg.List[sta.State]):
    connected_states = [state]
    not_connected_states = []
    new_states = []
    for connected_state in states:
        if {connected_state.first_component.name, connected_state.second_component.name} & \
                {state.first_component.name, state.second_component.name}:
            connected_states.append(connected_state)
        else:
            not_connected_states.append(connected_state)

    #for state in not_connected_states:
    #    states, connected_states = complexes_by_state_connection(new_states)
    return new_states, connected_states




