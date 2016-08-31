from typing import Dict, List

from rxncon.core.contingency import ContingencyType
from rxncon.core.effector import Effector, AndEffector, OrEffector, StateEffector, NotEffector
from rxncon.core.reaction import Reaction, Reactant
from rxncon.core.specification import Spec
from rxncon.core.state import State
from rxncon.core.rxncon_system import RxnConSystem
from rxncon.semantics.molecule import MoleculeDef, Molecule
from rxncon.semantics.molecule_from_rxncon import molecule_defs_from_rxncon
from rxncon.semantics.rule import Rule, Arrow, Complex
from rxncon.venntastic.sets import UniversalSet, Intersection, MultiUnion, MultiIntersection, ValueSet, Complement, Union, Set


def rules_from_reaction(rxnconsys: RxnConSystem, reaction: Reaction):
    state_set = UniversalSet()

    for contingency in rxnconsys.strict_contingencies(reaction):
        if contingency.type is ContingencyType.requirement:
            state_set = Intersection(state_set, state_set_from_effector(contingency.effector))
        elif contingency.type is ContingencyType.inhibition:
            state_set = Intersection(state_set, Complement(state_set_from_effector(contingency.effector)))
        else:
            raise NotImplementedError

    state_solutions = state_set.to_intersection_terms()
    consumed_solutions = []
    disjunct_solutions = []

    connected_components = [reactant.spec for reactant in reaction.reactants_lhs]

    while state_solutions:
        solution = state_solutions.pop(0)
        if is_connected(solution, connected_components):
            consumed_solutions.append(solution)
            connected_components = update_connected_components(connected_components, solution)

            for other_solution in consumed_solutions:
                solution = Intersection(solution, Complement(other_solution))

            for term in solution.to_intersection_terms():
                disjunct_solutions.append(term)
        else:
            state_solutions.append(solution)

    rules = []

    for solution in disjunct_solutions:
        rules += rules_from_solution(rxnconsys, reaction, solution)

    return rules


def rules_from_solution(rxnconsys: RxnConSystem, reaction: Reaction, solution: Set):
    def molecule_set_from_state_set(molecule_defs, state_set):
        assert state_set.is_single_intersection()
        state_sets = state_set.to_nested_lists()[0]

        molecule_set = UniversalSet()

        for state_set in state_sets:
            if isinstance(state_set, ValueSet):
                molecule_set = Intersection(molecule_set, molecule_set_from_state(molecule_defs, state_set.value))
            elif isinstance(state_set, Complement):
                molecule_set = Intersection(molecule_set,
                                            molecule_set_from_state_complement(molecule_defs, state_set.expr.value))

        return molecule_set

    molecule_defs = molecule_defs_from_rxncon(rxnconsys)
    molecule_set = molecule_set_from_state_set(molecule_defs, solution)

    rules = []

    for molecules in molecule_set.values:
        molecules = merge_identical_molecules(molecules)

        lhs = complexify_molecules_reactants(molecule_defs, molecules, reaction.reactants_lhs)
        rhs = complexify_molecules_reactants(molecule_defs, molecules, reaction.reactants_rhs)

        rules.append(Rule(lhs, rhs, Arrow.reversible, set()))

    return rules


def complexify_molecules_reactants(molecule_defs: Dict[Spec, MoleculeDef], molecules: List[Molecule],
                                   reactants: List[Reactant]):
    def reactant_to_molecule(reactant: Reactant, structure_index: int) -> Molecule:
        molecule = Molecule(molecule_defs[reactant.spec], structure_index)
        [molecule.set_state(state) for state in reactant.value]
        return molecule

    complexes = [Complex(set(reactant_to_molecule(reactant, index))) for index, reactant in enumerate(reactants)]

    while molecules:
        molecule = molecules.pop(0)
        for complex in complexes:
            if complex.can_bind_molecule(molecule):
                complex.add_molecule(molecule)
                break


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
        return ValueSet(effector.expr)
    elif isinstance(effector, AndEffector):
        return Intersection(state_set_from_effector(effector.left_expr), state_set_from_effector(effector.right_expr))
    elif isinstance(effector, OrEffector):
        return Union(state_set_from_effector(effector.left_expr), state_set_from_effector(effector.right_expr))
    elif isinstance(effector, NotEffector):
        return Complement(state_set_from_effector(effector.expr))
    else:
        raise NotImplementedError


def molecule_set_from_state(molecule_defs: Dict[Spec, MoleculeDef], state: State):
    molecule_set = UniversalSet()
    for spec in state.specs:
        molecule_def = molecule_defs[spec.to_component_spec()]

        matching_molecules = []

        for matching_state in [x for x in molecule_def.states if x.is_subset_of(state)]:
            molecule = Molecule(molecule_def, spec)
            molecule.set_state(matching_state)
            matching_molecules.append(molecule)

        molecule_set = Intersection(molecule_set, MultiUnion(*[ValueSet(molecule) for molecule in matching_molecules]))

    return molecule_set.simplified_form()


def molecule_set_from_state_complement(molecule_defs: Dict[Spec, MoleculeDef], state: State):
    molecule_set = UniversalSet()

    for state_spec in state.specs:
        molecule_def = molecule_defs[state_spec.to_component_spec()]

        for molecule_spec in [x for x in molecule_def.specs if state_spec.is_superset_of(x)
                              and any(state.is_superset_of(y) for y in molecule_def.states_for_spec(x))]:
            complementary_states = [x for x in molecule_def.states_for_spec(molecule_spec) if not state.is_superset_of(x)]

            complementary_molecules = []

            for complementary_state in complementary_states:
                molecule = Molecule(molecule_def, state_spec)
                molecule.set_state(complementary_state)
                complementary_molecules.append(molecule)

            molecule_set = Intersection(molecule_set,
                                        MultiUnion(*[ValueSet(molecule) for molecule in complementary_molecules]))

    return molecule_set.simplified_form()
