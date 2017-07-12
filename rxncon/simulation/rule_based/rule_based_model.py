from typing import Dict, List, Optional, Tuple, Iterable, Iterator  # pylint: disable=unused-import
from itertools import combinations, product, chain, permutations, count
from copy import copy, deepcopy
from collections import defaultdict, OrderedDict
from re import match
import logging

from rxncon.core.rxncon_system import RxnConSystem
from rxncon.core.reaction import Reaction, ReactionTerm, OutputReaction
from rxncon.core.state import State, StateModifier, ModificationState, InteractionState, SelfInteractionState, \
    GlobalState, \
    EmptyBindingState
from rxncon.core.spec import Spec
from rxncon.core.contingency import Contingency, ContingencyType
from rxncon.venntastic.sets import Set as VennSet, Intersection, Union, Complement, ValueSet, UniversalSet, \
    DisjunctiveUnion
from rxncon.util.utils import current_function_name

LOGGER = logging.getLogger(__name__)

NEUTRAL_MOD = '0'
INITIAL_MOLECULE_COUNT = 1000
SITE_NAME_REGEX = r'^[a-zA-Z0-9]+$'


STATE_TO_COMPLEX_BUILDER_FN = {
    ModificationState: [
        lambda state, builder: builder.add_mol(state.spec),
        lambda state, builder: builder.set_mod(state.spec, state.modifier)
    ],
    InteractionState: [
        lambda state, builder: builder.add_mol(state.first),
        lambda state, builder: builder.add_mol(state.second),
        lambda state, builder: builder.set_bond(state.first, state.second)
    ],
    SelfInteractionState: [
        lambda state, builder: builder.add_mol(state.first),
        lambda state, builder: builder.set_bond(state.first, state.second)
    ],
    EmptyBindingState: [
        lambda state, builder: builder.add_mol(state.spec),
        lambda state, builder: builder.set_half_bond(state.spec, None)
    ],
    GlobalState: [
        lambda state, builder: LOGGER.warning(
            '{} : IGNORING INPUT STATE {}'.format(current_function_name(), str(state)))
    ]
}  # type: ignore

STATE_TO_MOL_DEF_BUILDER_FN = {
    ModificationState: [
        lambda state, builder: builder.add_site(state.spec),
        lambda state, builder: builder.add_mod(state.spec, state.modifier)
    ],
    InteractionState: [
        lambda state, builder: builder.add_site(state.first) if builder.spec == state.first.to_component_spec()
        else builder.add_site(state.second),
    ],
    SelfInteractionState: [
        lambda state, builder: builder.add_site(state.first),
        lambda state, builder: builder.add_site(state.second)
    ],
    EmptyBindingState: [
        lambda state, builder: builder.add_site(state.spec)
    ]
}  # type: ignore


class MolDef:
    def __init__(self, name: str, site_defs: Dict[str, List[str]]) -> None:
        self.name, self.site_defs = name, site_defs

    def __str__(self) -> str:
        site_strs = []
        for name in sorted(self.site_defs.keys()):
            site_def = self.site_defs[name]
            if site_def:
                site_strs.append('{0}:{1}'.format(name, '~'.join(x for x in sorted(site_def))))
            else:
                site_strs.append(name)
        return '{0}({1})'.format(self.name, ','.join(site_strs))

    def __repr__(self) -> str:
        return 'MolDef<{}>'.format(str(self))

    def mods_for_site(self, site: str) -> List[str]:
        return self.site_defs[site]

    @property
    def sites(self) -> List[str]:
        return list(self.site_defs.keys())

    def create_neutral_complex(self) -> 'Complex':
        site_to_mod = {}  # type: Dict[str, str]
        site_to_bond = {}  # type: Dict[str, Optional[int]]
        for site, mods in self.site_defs.items():
            if mods:
                site_to_mod[site] = NEUTRAL_MOD
            else:
                site_to_bond[site] = None

        return Complex([Mol(self.name, site_to_mod, site_to_bond, False)])


class MolDefBuilder:
    def __init__(self, spec: Spec) -> None:
        self.spec = spec
        self.name = str(spec.to_non_struct_spec())  # type: str
        self.site_defs = {}  # type: Dict[str, List[str]]

    def build(self) -> MolDef:
        return MolDef(self.name, self.site_defs)

    def add_site(self, site: Spec) -> None:
        if site_name(site) not in self.site_defs:
            self.site_defs[site_name(site)] = []

    def add_mod(self, site: Spec, mod: StateModifier) -> None:
        self.site_defs[site_name(site)].append(str(mod.value))


class Mol:
    def __init__(self, name: str, site_to_mod: Dict[str, str], site_to_bond: Dict[str, Optional[int]],
                 is_reactant: bool) -> None:
        self.name = name
        self.site_to_mod = site_to_mod
        self.site_to_bond = site_to_bond
        self.is_reactant = is_reactant
        self._assert_valid_site_names()

    def __str__(self) -> str:
        mod_str = ','.join('{}{}'.format(site, '~' + mod)
                           for site, mod in sorted(self.site_to_mod.items()))
        bond_str = ','.join('{}{}'.format(site, '!' + str(bond) if bond is not None else '')
                            for site, bond in sorted(self.site_to_bond.items()))

        strs = []
        if mod_str:
            strs.append(mod_str)
        if bond_str:
            strs.append(bond_str)

        return '{}({})'.format(self.name, ','.join(strs))

    def __repr__(self) -> str:
        return 'Mol<{}>'.format(str(self))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Mol):
            return NotImplemented
        # This equality skips checking the is_reactant property, since it does not appear
        # in the BNGL representation.
        return self.name == other.name and self.site_to_mod == other.site_to_mod and self.site_to_bond == other.site_to_bond

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, Mol):
            return NotImplemented
        return str(self) < str(other)

    def clone(self) -> 'Mol':
        return deepcopy(self)

    def has_bond(self, bond: int) -> bool:
        return bond in self.bonds

    @property
    def sites(self) -> List[str]:
        return list(set(list(self.site_to_mod.keys()) + list(self.site_to_bond.keys())))

    @property
    def bonds(self) -> List[int]:
        return [x for x in self.site_to_bond.values() if x is not None]

    def sites_by_bond(self, bond: int) -> List[str]:
        sites = [s for s, b in self.site_to_bond.items() if b is not None and bond == b]
        return sites

    def with_relabeled_bonds(self, bond_to_bond: Dict[int, int]) -> 'Mol':
        new_mol = self.clone()
        for site, old_bond in new_mol.site_to_bond.items():
            if old_bond is not None:
                new_mol.site_to_bond[site] = bond_to_bond[old_bond]

        return new_mol

    def _assert_valid_site_names(self) -> None:
        for site in self.sites:
            assert match(SITE_NAME_REGEX, site), 'Invalid site name: {}'.format(site)


def site_name(spec: Spec) -> str:
    bad_chars = ['[', ']', '/', '(', ')', ':', '-']
    spec_str = (spec.locus.domain + 'D' if spec.locus.domain else '') + \
               (spec.locus.subdomain + 'S' if spec.locus.subdomain else '') + \
               (spec.locus.residue + 'R' if spec.locus.residue else '')

    for bad_char in bad_chars:
        spec_str = spec_str.replace(bad_char, '')

    return spec_str


class MolBuilder:
    def __init__(self, spec: Spec, is_reactant: bool = False) -> None:
        self.name = str(spec.to_non_struct_spec())
        self.site_to_mod = {}  # type: Dict[str, str]
        self.site_to_bond = {}  # type: Dict[str, Optional[int]]
        self.is_reactant = is_reactant

    def build(self) -> Mol:
        return Mol(self.name, self.site_to_mod, self.site_to_bond, self.is_reactant)

    def set_bond_index(self, spec: Spec, bond_index: Optional[int]) -> None:
        self.site_to_bond[site_name(spec)] = bond_index

    def set_mod(self, spec: Spec, mod: StateModifier) -> None:
        self.site_to_mod[site_name(spec)] = str(mod.value)


class Complex:
    def __init__(self, mols: List[Mol]) -> None:
        self.mols = sorted(mols)

        self._assert_bonds_connected()
        self._assert_mols_connected()

    def __str__(self) -> str:
        return '.'.join(str(mol) for mol in self.mols)

    def __repr__(self) -> str:
        return 'Complex<{}>'.format(str(self))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Complex):
            return NotImplemented
        return self.mols == other.mols

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, Complex):
            return NotImplemented
        return self.mols < other.mols

    @property
    def bonds(self) -> List[int]:
        return sorted(set(bond for mol in self.mols for bond in mol.bonds))

    def neighbors(self, mol: Mol) -> List[Mol]:
        neighbors = []
        for bond in mol.bonds:
            mols = self.mols_by_bond(bond)
            if len(mols) == 1:
                # Intra-particle bond
                continue
            assert len(mols) == 2
            if mols[0] == mols[1]:
                # Homodimer
                neighbors.append(mols[0])
                continue
            if mols[0] == mol:
                neighbors.append(mols[1])
            elif mols[1] == mol:
                neighbors.append(mols[0])
            else:
                raise AssertionError

        return neighbors

    def mols_by_bond(self, bond: int) -> List[Mol]:
        return [mol for mol in self.mols if mol.has_bond(bond)]

    def is_equivalent_to(self, other: 'Complex') -> bool:
        if self == other:
            return True
        if len(self.mols) != len(other.mols) or \
                        sorted(mol.name for mol in self.mols) != sorted(mol.name for mol in other.mols) or \
                        len(self.bonds) != len(other.bonds):
            return False

        my_bonds = self.bonds

        for new_bonds in permutations(other.bonds):
            bond_to_bond = {other_bond: my_bond for other_bond, my_bond in zip(new_bonds, my_bonds)}
            if self == other.with_relabeled_bonds(bond_to_bond):
                return True

        return False

    def with_relabeled_bonds(self, bond_to_bond: Dict[int, int]) -> 'Complex':
        return Complex([mol.with_relabeled_bonds(bond_to_bond) for mol in self.mols])

    @property
    def is_reactant(self) -> bool:
        return any(mol.is_reactant for mol in self.mols)

    def _assert_mols_connected(self) -> None:
        connected = []  # type: List[Mol]
        to_visit = [self.mols[0]]
        while to_visit:
            current = to_visit.pop()
            neighbors = self.neighbors(current)
            for neighbor in neighbors:
                if neighbor not in connected and neighbor not in to_visit:
                    to_visit.append(neighbor)

            connected.append(current)

        assert sorted(connected) == self.mols

    def _assert_bonds_connected(self) -> None:
        bond_to_site_count = defaultdict(int)  # type: Dict[int, int]

        for mol in self.mols:
            for bond in mol.bonds:
                bond_to_site_count[bond] += len(mol.sites_by_bond(bond))

        non_doubly_connected_bonds = [bond for bond, site_count in bond_to_site_count.items() if site_count % 2 != 0]
        if non_doubly_connected_bonds:
            raise AssertionError(
                'Non doubly connected bonds {} in complex: {}'.format(non_doubly_connected_bonds, str(self)))


class ComplexExprBuilder:
    def __init__(self, reaction: Optional[Reaction] = None) -> None:
        self._mol_builders = {}  # type: Dict[Spec, MolBuilder]
        self._current_bond = 0
        self._bonds = []  # type: List[Tuple[Spec, Spec]]
        self.reaction = reaction

    def build(self, only_reactants: bool = True) -> List[Complex]:
        complexes = []

        LOGGER.debug('{} : Building complex with molecules {}'.format(current_function_name(),
                                                                      ' & '.join(str(spec) for spec in
                                                                                 self._mol_builders.keys())))
        LOGGER.debug('{} : Grouped specs are {}'.format(current_function_name(),
                                                        ', '.join(str(x) for x in self._grouped_specs())))
        for group in self._grouped_specs():
            possible_complex = Complex([self._mol_builders[spec].build() for spec in group])

            if not only_reactants or (only_reactants and possible_complex.is_reactant):
                complexes.append(possible_complex)
                LOGGER.debug('{} : Adding complex {}'.format(current_function_name(), possible_complex))
            else:
                LOGGER.info('{} : DISCONNECTED CONTINGENCY Reaction {} / Not adding complex {}'
                            .format(current_function_name(), self.reaction, possible_complex))

        return complexes

    def add_mol(self, spec: Spec, is_reactant: bool = False) -> None:
        if spec.to_component_spec() not in self._mol_builders:
            self._mol_builders[spec.to_component_spec()] = MolBuilder(spec.to_component_spec(), is_reactant)

    def set_bond(self, first: Spec, second: Spec) -> None:
        if (first, second) not in self._bonds and (second, first) not in self._bonds:
            self._bonds.append((first, second))
            self._current_bond += 1
            self.set_half_bond(first, self._current_bond)
            self.set_half_bond(second, self._current_bond)

    def set_half_bond(self, spec: Spec, value: Optional[int]) -> None:
        self._mol_builders[spec.to_component_spec()].set_bond_index(spec, value)

    def set_mod(self, spec: Spec, mod: StateModifier) -> None:
        self._mol_builders[spec.to_component_spec()].set_mod(spec, mod)

    def _grouped_specs(self) -> List[List[Spec]]:
        grouped_specs = [[spec] for spec in self._mol_builders.keys()]
        for bond in self._bonds:
            for i, group in enumerate(grouped_specs):
                if bond[0].to_component_spec() in group:
                    grouped_specs.pop(i)
                    try:
                        other_group = next(
                            other_group for other_group in grouped_specs if bond[1].to_component_spec() in other_group)
                        other_group += group
                        break
                    except StopIteration:
                        grouped_specs.insert(i, group)
                        break

        return grouped_specs


class Parameter:
    def __init__(self, name: str, value: Optional[str]) -> None:
        self.name, self.value = name, value

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Parameter):
            return NotImplemented
        return self.name == other.name and self.value == other.value

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return str(self)


class InitialCondition:
    def __init__(self, complex: Complex, value: Parameter) -> None:
        self.complex, self.value = complex, value

    def __str__(self) -> str:
        return '{} {}'.format(str(self.complex), str(self.value))

    def __repr__(self) -> str:
        return 'InitialCondition<{}>'.format(str(self))

    def is_equivalent_to(self, other: 'InitialCondition') -> bool:
        return self.complex.is_equivalent_to(other.complex) and self.value.name == other.value.name


class Observable:
    def __init__(self, name: str, complex: Complex) -> None:
        self.name, self.complex = name, complex

    def __str__(self) -> str:
        return '{} {}'.format(self.name, str(self.complex))

    def __repr__(self) -> str:
        return 'Observable<{}>'.format(str(self))


class Rule:
    def __init__(self, lhs: List[Complex], rhs: List[Complex], rate: Parameter,
                 parent_reaction: Reaction = None) -> None:
        self.lhs, self.rhs, self.rate = sorted(lhs), sorted(rhs), rate
        self.parent_reaction = parent_reaction

    def __str__(self) -> str:
        return ' + '.join(str(x) for x in self.lhs) + ' -> ' + ' + '.join(str(x) for x in self.rhs) + \
               ' ' + str(self.rate) + ' ' + str(self.parent_reaction)

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Rule):
            return NotImplemented
        return self.lhs == other.lhs and self.rhs == other.lhs and \
               self.rate == other.rate and self.parent_reaction == other.parent_reaction

    def is_equivalent_to(self, other: 'Rule') -> bool:
        if len(self.lhs) != len(other.lhs) or len(self.rhs) != len(other.rhs) or self.rate.name != other.rate.name:
            return False

        if self == other:
            return True

        my_lhs = deepcopy(self.lhs)
        my_rhs = deepcopy(self.rhs)

        other_lhs = deepcopy(other.lhs)
        other_rhs = deepcopy(other.rhs)

        for complex in my_lhs + my_rhs:  # pylint: disable=redefined-builtin
            complex.found = False  # type: ignore

        for complex in other_lhs + other_rhs:
            complex.visited = False  # type: ignore

        for my_complex in my_lhs:
            for other_complex in other_lhs:
                if my_complex.is_equivalent_to(other_complex) and not other_complex.visited:  # type: ignore
                    my_complex.found = True  # type: ignore
                    other_complex.visited = True  # type: ignore
                    break

        for my_complex in my_rhs:
            for other_complex in other_rhs:
                if my_complex.is_equivalent_to(other_complex) and not other_complex.visited:  # type: ignore
                    my_complex.found = True  # type: ignore
                    other_complex.visited = True  # type: ignore
                    break

        all_found = all(complex.found for complex in my_lhs + my_rhs)  # type: ignore
        all_visited = all(complex.visited for complex in other_lhs and other_rhs)  # type: ignore

        return all_found and all_visited


class RuleBasedModel:
    def __init__(self, mol_defs: List[MolDef], initial_conditions: List[InitialCondition], parameters: List[Parameter],
                 observables: List[Observable], rules: List[Rule]) -> None:
        self.mol_defs, self.initial_conditions, self.parameters, self.observables, self.rules = mol_defs, initial_conditions, \
                                                                                                parameters, observables, rules

    @property
    def rate_parameters(self) -> List[Parameter]:
        return sorted(set(rule.rate for rule in self.rules), key=lambda x: x.name)


def calc_state_paths(states: List[State]) -> Dict[State, List[List[State]]]:
    specs = sorted(set(spec.to_component_spec() for state in states for spec in state.specs),
                   key=lambda x: x.struct_index)

    def spec_to_bond_states(spec: Spec) -> List[State]:
        assert spec.is_component_spec and spec.is_structured
        return [state for state in states if spec in (s.to_component_spec() for s in state.specs) if
                len(state.components) == 2]

    def neighbor(spec: Spec, state: State) -> Spec:
        assert spec.is_component_spec and spec.is_structured and len(state.components) == 2
        return [neigh_spec.to_component_spec() for neigh_spec in state.specs if neigh_spec.to_component_spec() != spec][
            0]

    spec_paths = {spec: [] for spec in specs}  # type: Dict[Spec, List[List[State]]]

    for state in states:
        state.visited = []  # type: ignore

    nums = count(1)

    def visit_nodes(current_path: List[State], current_spec: Spec, current_num: int) -> None:
        spec_paths[current_spec].append(current_path)
        bonds_to_visit = [state for state in spec_to_bond_states(current_spec) if
                          current_num not in state.visited]  # type: ignore
        for i, state in enumerate(bonds_to_visit):
            if i == 0:
                next_num = current_num
            else:
                next_num = next(nums)

            for visited_state in current_path + [state]:
                visited_state.visited.append(next_num)  # type: ignore

            visit_nodes(current_path + [state], neighbor(current_spec, state), next_num)

    for spec in specs:
        assert spec.struct_index is not None
        if spec.struct_index > 1:
            break
        visit_nodes([], spec, 0)

    for spec, paths in spec_paths.items():
        if not paths:
            raise AssertionError('Could not find path to {}'.format(str(spec)))

    state_paths = {state: [] for state in states}  # type: Dict[State, List[List[State]]]
    for state in states:
        for component in state.components:
            for spec_path in spec_paths[component]:
                if len(spec_path) > 0 and spec_path[-1] == state:
                    path = spec_path[:-1]
                else:
                    path = spec_path

                if path not in state_paths[state]:
                    state_paths[state].append(path)

    for state in states:
        del state.visited  # type: ignore

    return state_paths


def calc_connected_complexes(states: List[State]) -> List[List[State]]:
    complexes = [[]]  # type: List[List[State]]

    while states:
        state = states.pop()

        if state.is_global:
            continue

        for complex in complexes:  # pylint: disable=redefined-builtin
            if state in complex:
                continue
            elif not any(state.is_mutually_exclusive_with(other) for other in complex):
                complex.append(state)
            else:
                other = next(other for other in complex if state.is_mutually_exclusive_with(other))
                new_complex = deepcopy(complex)
                new_complex.remove(other)
                new_complex.append(state)
                complexes.append(new_complex)

    # Make sure that the complexes are completely connected to the reactants.
    # The function calc_state_paths raises an exception if it's not.
    connected_complexes = []
    for complex in complexes:
        try:
            calc_state_paths(complex)
            connected_complexes.append(complex)
        except AssertionError:
            pass

    return connected_complexes


def with_connectivity_constraints(cont_set: VennSet[State]) -> VennSet:
    complexes = calc_connected_complexes(cont_set.values)
    complex_constraints = []

    for complex in complexes:  # pylint: disable=redefined-builtin
        state_paths = calc_state_paths(complex)
        constraint = UniversalSet()  # type:  VennSet[State]

        for state in complex:
            assert not state.is_global, 'Global state {} appearing in connectivity constraints.'.format(state)

            if any(path == [] for path in state_paths[state]):
                continue

            state_constraints = [Complement(ValueSet(state))]  # type: List[VennSet[State]]
            for path in state_paths[state]:
                state_constraints.append(Intersection(*(ValueSet(x) for x in path)))

            constraint = Intersection(constraint, Union(*state_constraints))  # pylint: disable=redefined-variable-type

        complex_constraints.append(constraint.to_simplified_set())

    if complex_constraints:
        LOGGER.debug('{} : Complex constraints {}'.format(current_function_name(),
                                                          ' XOR '.join(str(x) for x in complex_constraints)))
        return Intersection(cont_set, DisjunctiveUnion(*complex_constraints))
    else:
        return cont_set


class QuantContingencyConfigs(Iterator[VennSet[State]]):  # pylint: disable=too-few-public-methods
    def __init__(self, q_contingencies: List[Contingency]) -> None:
        self.q_contingencies = deepcopy(q_contingencies)

        combis = [[]]  # type: List[List[Contingency]]
        for contingency in self.q_contingencies:
            new_combis = []
            for combi in combis:
                new_combis.append(
                    combi + [Contingency(contingency.reaction, ContingencyType.inhibition, contingency.effector)])
                combi.append(Contingency(contingency.reaction, ContingencyType.requirement, contingency.effector))
            combis.extend(new_combis)

        self.combi_sets = []  # type: List[VennSet[State]]

        if combis == [[]]:
            self.combi_sets = [UniversalSet()]
        else:
            self.combi_sets = [Intersection(*(x.to_venn_set() for x in combi)) for combi in combis]

        self.current_combi_set = -1

    def __iter__(self) -> 'QuantContingencyConfigs':
        return self

    def __next__(self) -> VennSet[State]:
        try:
            self.current_combi_set += 1
            return self.combi_sets[self.current_combi_set]
        except IndexError:
            raise StopIteration


def rule_based_model_from_rxncon(rxncon_sys: RxnConSystem) -> RuleBasedModel:  # pylint: disable=too-many-locals
    def mol_defs_from_rxncon(rxncon_sys: RxnConSystem) -> List[MolDef]:
        mol_defs = {}
        for spec in rxncon_sys.components():
            LOGGER.debug('{} : Creating MolDefBuilder for {}'.format(current_function_name(), str(spec)))
            builder = MolDefBuilder(spec)
            for state in rxncon_sys.states_for_component(spec):
                LOGGER.debug(
                    '{} : Applying State {} of type {}'.format(current_function_name(), str(state), type(state)))
                for func in STATE_TO_MOL_DEF_BUILDER_FN[type(state)]:  # type: ignore
                    func(state, builder)

            mol_defs[spec] = builder.build()

        return list(mol_defs.values())

    def remove_global_states(solutions: List[Dict[State, bool]]) -> List[Dict[State, bool]]:
        filtered_solutions = []  # type: List[Dict[State, bool]]
        for soln in solutions:
            cleaned_solution = {}
            for state, val in soln.items():
                if state.is_global:
                    LOGGER.warning(
                        '{} : REMOVING INPUT STATE {} from contingencies.'.format(current_function_name(), state))
                else:
                    cleaned_solution[state] = val

            if cleaned_solution not in filtered_solutions:
                filtered_solutions.append(cleaned_solution)

        return filtered_solutions

    def is_satisfiable(states: Iterable[State]) -> bool:
        for pair in combinations(states, 2):
            if pair[0].is_mutually_exclusive_with(pair[1]):
                return False

        return True

    def calc_positive_solutions(rxncon_sys: RxnConSystem, solution: Dict[State, bool]) -> List[List[State]]:
        def complementary_state_combos(state: State) -> List[List[State]]:
            combos = product(
                *(rxncon_sys.complement_states_for_component(spec.to_component_spec(), state) for spec in
                  state.specs))
            return [list(combo) for combo in combos if is_satisfiable(combo)]

        def structure_states(states: List[State]) -> List[State]:
            cur_index = max(spec.struct_index for state in states for spec in state.specs if spec.is_structured)
            assert cur_index is not None

            spec_to_index = {}  # type: Dict[Spec, int]
            struct_states = []  # type: List[State]

            for state in states:
                if state.is_structured:
                    struct_states.append(state)
                    continue

                for spec in state.specs:
                    if spec.is_structured:
                        continue

                    try:
                        state = state.to_structured_from_spec(
                            spec.with_struct_index(spec_to_index[spec.to_component_spec()]))
                    except KeyError:
                        cur_index += 1
                        state = state.to_structured_from_spec(spec.with_struct_index(cur_index))
                        spec_to_index[spec.to_component_spec()] = cur_index

                struct_states.append(state)

            return struct_states

        ordered_solution = OrderedDict(sorted(solution.items(), key=lambda x: x[0]))

        trues = [state for state, val in ordered_solution.items() if val]
        falses = [state for state, val in ordered_solution.items() if not val
                  and not any(state.is_mutually_exclusive_with(x) for x in trues)]

        if not falses:
            return [trues] if is_satisfiable(trues) else []

        positivized_falses = [list(chain(*x)) for x in
                              product(*(complementary_state_combos(state) for state in falses))]

        solutions = []

        for positivized_false in positivized_falses:
            possible_solution = []  # type: List[State]
            for soln_state in structure_states(trues + positivized_false):
                if soln_state not in possible_solution:
                    possible_solution.append(soln_state)

            if is_satisfiable(possible_solution):
                solutions.append(possible_solution)

        return solutions

    def calc_rule(reaction: Reaction, cont_soln: List[State]) -> Rule:
        def calc_complexes(terms: List[ReactionTerm], states: List[State]) -> List[Complex]:
            if not all(x.is_structured for x in states):
                unstructs = [x for x in states if not x.is_structured]
                raise AssertionError('Error in building rule for Reaction {}, States {} appear unstructured'
                                     .format(str(reaction), ', '.join(str(x) for x in unstructs)))

            if not is_satisfiable(cont_soln):
                raise AssertionError(
                    'Cannot satisfy contingencies {} simultaneously'.format(' & '.join(str(s) for s in cont_soln)))

            states = copy(states)
            builder = ComplexExprBuilder(reaction=reaction)
            struct_index = 0
            for term in terms:
                struct_states = deepcopy(term.states)
                for spec in term.specs:
                    struct_spec = copy(spec)
                    struct_spec.struct_index = struct_index
                    builder.add_mol(struct_spec, is_reactant=True)
                    struct_states = [state.to_structured_from_spec(struct_spec) for state in struct_states]
                    struct_index += 1

                states += struct_states

            assert all(x.is_structured for x in states)

            for state in states:
                for func in STATE_TO_COMPLEX_BUILDER_FN[type(state)]:  # type: ignore
                    func(state, builder)

            return builder.build()

        lhs = calc_complexes(reaction.terms_lhs, cont_soln)
        rhs = calc_complexes(reaction.terms_rhs, cont_soln)

        rate = Parameter('k', '1.0')

        return Rule(lhs, rhs, rate, parent_reaction=reaction)

    def calc_initial_conditions(mol_defs: List[MolDef]) -> List[InitialCondition]:
        return \
            [InitialCondition(mol_def.create_neutral_complex(),
                              Parameter('Num{}'.format(mol_def.name), str(INITIAL_MOLECULE_COUNT))) for mol_def in
             mol_defs]

    def calc_observables(rxncon_sys: RxnConSystem) -> List[Observable]:
        def observable_complex(states: List[State]) -> Complex:
            builder = ComplexExprBuilder()

            assert all(x.is_structured for x in states), 'Unstructured states appearing in observable ' \
                                                         '{}'.format(states)

            for state in states:
                for func in STATE_TO_COMPLEX_BUILDER_FN[type(state)]:  # type: ignore
                    func(state, builder)

            complexes = builder.build(only_reactants=False)

            assert len(complexes) == 1, 'Multiple complexes appearing in observable {}'.format(states)
            return complexes[0]

        observables = []
        output_rxns = [rxn for rxn in rxncon_sys.reactions if isinstance(rxn, OutputReaction)]
        for rxn in output_rxns:
            LOGGER.debug('{} : calculating observable {}'.format(current_function_name(), str(rxn)))
            solns = Intersection(*(x.to_venn_set() for x
                                   in rxncon_sys.contingencies_for_reaction(rxn))).calc_solutions()
            positive_solns = []  # type: List[List[State]]
            for soln in solns:
                positive_solns += calc_positive_solutions(rxncon_sys, soln)

            for index, positive_soln in enumerate(positive_solns):
                LOGGER.debug('{} : solution {} : {}'.format(current_function_name(), index, positive_soln))
                observables.append(Observable('{}{}'.format(rxn.name, index), observable_complex(positive_soln)))

        return observables

    LOGGER.debug('{} : Entered function'.format(current_function_name()))

    mol_defs = mol_defs_from_rxncon(rxncon_sys)
    LOGGER.debug(
        '{} : Generated MolDefs: {}'.format(current_function_name(), ', '.join(str(mol_def) for mol_def in mol_defs)))

    rules = []  # type: List[Rule]

    for reaction in (x for x in rxncon_sys.reactions if not isinstance(x, OutputReaction)):
        LOGGER.debug('{} : Generating rules for reaction {}'.format(current_function_name(), str(reaction)))
        strict_cont_set = Intersection(*(x.to_venn_set() for x in
                                         rxncon_sys.s_contingencies_for_reaction(reaction)))  # type: VennSet[State]
        quant_contingencies = QuantContingencyConfigs(rxncon_sys.q_contingencies_for_reaction(reaction))
        LOGGER.debug('{} : Strict contingencies {}'.format(current_function_name(), str(strict_cont_set)))

        for quant_contingency_set in quant_contingencies:
            LOGGER.debug(
                '{} : quantitative contingency config: {}'.format(current_function_name(), str(quant_contingency_set)))

            cont_set = Intersection(strict_cont_set, quant_contingency_set)  # type: VennSet[State]
            cont_set = with_connectivity_constraints(cont_set)
            solutions = cont_set.calc_solutions()
            solutions = remove_global_states(solutions)

            LOGGER.debug('{} : contingency solutions {}'.format(current_function_name(), str(solutions)))
            positive_solutions = []  # type: List[List[State]]
            for solution in solutions:
                positive_solutions += calc_positive_solutions(rxncon_sys, solution)

            for positive_solution in positive_solutions:
                LOGGER.debug('{} : positivized contingency solution {}'
                             .format(current_function_name(), ' & '.join(str(x) for x in positive_solution)))
                rule = calc_rule(reaction, positive_solution)
                if not any(rule.is_equivalent_to(existing) for existing in rules):
                    rules.append(rule)

    return RuleBasedModel(mol_defs, calc_initial_conditions(mol_defs), [], calc_observables(rxncon_sys), rules)


def mol_from_str(mol_str: str) -> Mol:
    name, config_str = mol_str.split('(')
    assert config_str[-1] == ')'
    site_strs = config_str[:-1].split(',')

    site_to_mod = {}  # type: Dict[str, str]
    site_to_bond = {}  # type: Dict[str, Optional[int]]

    for site_str in site_strs:
        if not site_str:
            continue
        elif '~' not in site_str and '!' not in site_str:
            site_to_bond[site_str] = None
        elif '~' not in site_str:
            site, bond = site_str.split('!')
            site_to_bond[site] = int(bond)
        elif '!' not in site_str:
            site, mod = site_str.split('~')
            site_to_mod[site] = mod
        else:
            raise AssertionError

    return Mol(name, site_to_mod, site_to_bond, False)


def complex_from_str(complex_str: str) -> Complex:
    return Complex([mol_from_str(x) for x in complex_str.split('.')])


def rule_from_str(rule_str: str) -> Rule:
    lhs_strs = []  # type: List[str]
    rhs_strs = []  # type: List[str]
    rate = None  # type: Optional[Parameter]

    current = lhs_strs
    done = False

    for part in rule_str.split():
        if part == '+':
            done = False
            continue
        elif part == '->':
            done = False
            current = rhs_strs
        elif not done:
            current.append(part)
            done = True
        elif done:
            rate = Parameter(part, None)

    assert rate is not None
    return Rule([complex_from_str(x) for x in lhs_strs], [complex_from_str(x) for x in rhs_strs], rate)


def initial_condition_from_str(ic_str: str) -> InitialCondition:
    complex_str, param = ic_str.split()
    return InitialCondition(complex_from_str(complex_str), Parameter(param, None))
