from typing import Dict, List

from rxncon.core.contingency import ContingencyType
from rxncon.core.effector import Effector, AndEffector, OrEffector, StateEffector, NotEffector
from rxncon.core.reaction import Reaction
from rxncon.core.specification import Specification
from rxncon.core.state import State
from rxncon.core.rxncon_system import RxnConSystem
from rxncon.semantics.molecule import MoleculeDef, Molecule
from rxncon.semantics.molecule_from_rxncon import molecule_defs_from_rxncon
from rxncon.semantics.rule import Rule, Arrow
from rxncon.venntastic.sets import UniversalSet, Intersection, MultiUnion, MultiIntersection, PropertySet, Complement, Union, Set


def rules_from_reaction(rxnconsys: RxnConSystem, reaction: Reaction):
    state_set = UniversalSet()

    for contingency in rxnconsys.strict_contingencies_for_reaction(reaction):
        if contingency.type == ContingencyType.requirement:
            state_set = Intersection(state_set, state_set_from_effector(contingency.effector))
        elif contingency.type == ContingencyType.inhibition:
            state_set = Intersection(state_set, Complement(state_set_from_effector(contingency.effector)))
        else:
            raise NotImplementedError

    state_solutions = state_set.to_union_list_form()
    consumed_solutions = []
    disjunct_solutions = []

    connected_components = [reactant.component for reactant in reaction.reactants_pre]

    while state_solutions:
        solution = state_solutions.pop(0)
        if is_connected(solution, connected_components):
            consumed_solutions.append(solution)
            connected_components = update_connected_components(connected_components, solution)

            for other_solution in consumed_solutions:
                solution = Intersection(solution, Complement(other_solution))

            for term in solution.to_union_list_form():
                disjunct_solutions.append(term)
        else:
            state_solutions.append(solution)

    rules = []

    for solution in disjunct_solutions:
        rules += rules_from_solution(rxnconsys, reaction, solution)


    return rules


def rules_from_solution(rxnconsys: RxnConSystem, reaction: Reaction, solution: Set):
    molecule_defs = molecule_defs_from_rxncon(rxnconsys)

    assert len(solution.to_nested_list_form()) == 1
    state_sets = solution.to_nested_list_form()[0]

    molecule_set = UniversalSet()

    for state_set in state_sets:
        if isinstance(state_set, PropertySet):
            molecule_set = Intersection(molecule_set, molecule_set_from_state(molecule_defs, state_set.value))
        elif isinstance(state_set, Complement):
            molecule_set = Intersection(molecule_set, molecule_set_from_state_complement(molecule_defs, state_set.expr.value))

    rules = []

    for molecules in molecule_set.to_nested_list_form():
        molecules = merge_identical_molecules(molecules)
        complexes = complexes_from_molecules(molecules)
        lhs = combine_complexes_reactants(complexes, reaction.reactants_pre)
        rhs = combine_complexes_reactants(complexes, reaction.reactants_post)

        rules.append(Rule(lhs, rhs, Arrow.reversible, set()))

    return rules


def merge_identical_molecules(molecules: List[Molecule]) -> List[Molecule]:
    merged_molecules = []

    while molecules:
        molecule = molecules.pop()
        for other in molecules:
            if molecule.component_matches(other):
                molecule.merge_with(other)
                molecules.remove(other)

        merged_molecules.append(molecule)

    return merged_molecules


def state_set_from_effector(effector: Effector):
    if isinstance(effector, StateEffector):
        return PropertySet(effector.expr)
    elif isinstance(effector, AndEffector):
        return Intersection(state_set_from_effector(effector.left_expr), state_set_from_effector(effector.right_expr))
    elif isinstance(effector, OrEffector):
        return Union(state_set_from_effector(effector.left_expr), state_set_from_effector(effector.right_expr))
    elif isinstance(effector, NotEffector):
        return Complement(state_set_from_effector(effector.expr))
    else:
        raise NotImplementedError


def molecule_set_from_state(molecule_defs: Dict[Specification, MoleculeDef], state: State):
    molecule_set = UniversalSet()
    for spec in state.specs:
        molecule_def = molecule_defs[spec.to_component_specification()]

        matching_molecules = []

        for matching_state in [x for x in molecule_def.states if x.is_subset_of(state)]:
            molecule = Molecule(molecule_def, spec)
            molecule.set_state(matching_state)
            matching_molecules.append(molecule)

        molecule_set = Intersection(molecule_set, MultiUnion(*[PropertySet(molecule) for molecule in matching_molecules]))

    return molecule_set.simplified_form()


def molecule_set_from_state_complement(molecule_defs: Dict[Specification, MoleculeDef], state: State):
    molecule_set = UniversalSet()

    for state_spec in state.specs:
        molecule_def = molecule_defs[state_spec.to_component_specification()]

        for molecule_spec in [x for x in molecule_def.specs if state_spec.is_superset_of(x)
                              and any(state.is_superset_of(y) for y in molecule_def.states_for_spec(x))]:
            complementary_states = [x for x in molecule_def.states_for_spec(molecule_spec) if not state.is_superset_of(x)]

            complementary_molecules = []

            for complementary_state in complementary_states:
                molecule = Molecule(molecule_def, state_spec)
                molecule.set_state(complementary_state)
                complementary_molecules.append(molecule)

            molecule_set = Intersection(molecule_set,
                                        MultiUnion(*[PropertySet(molecule) for molecule in complementary_molecules]))

    return molecule_set.simplified_form()
