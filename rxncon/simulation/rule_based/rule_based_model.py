from typing import Dict, List, Optional, Tuple, Iterable
from itertools import combinations, product, chain, permutations, count
from copy import copy, deepcopy
from collections import defaultdict, OrderedDict
from re import match

from rxncon.core.rxncon_system import RxnConSystem
from rxncon.core.reaction import Reaction, ReactionTerm, OutputReaction
from rxncon.core.state import State, StateModifier
from rxncon.core.spec import Spec
from rxncon.core.contingency import Contingency, ContingencyType
from rxncon.core.effector import Effector, AndEffector, OrEffector, NotEffector, StateEffector
from rxncon.venntastic.sets import Set as VennSet, Intersection, Union, Complement, ValueSet, UniversalSet


NEUTRAL_MOD = '0'
INITIAL_MOLECULE_COUNT = 100
SITE_NAME_REGEX = '^[a-zA-Z0-9]+$'

class MolDef:
    def __init__(self, name: str, site_defs: Dict[str, List[str]]):
        self.name, self.site_defs = name, site_defs

    def __str__(self) -> str:
        site_strs = []
        for site_name, site_def in self.site_defs.items():
            site_strs.append('{0}:{1}'.format(site_name, '~'.join(x for x in site_def))) if site_def else site_strs.append(site_name)
        return '{0}({1})'.format(self.name, ','.join(site_strs))

    def __repr__(self) -> str:
        return 'MolDef<{}>'.format(str(self))

    def mods_for_site(self, site: str) -> List[str]:
        return self.site_defs[site]

    @property
    def sites(self) -> List[str]:
        return list(self.site_defs.keys())

    def create_neutral_complex(self) -> 'Complex':
        site_to_mod  = {}
        site_to_bond = {}
        for site, mods in self.site_defs.items():
            if mods:
                site_to_mod[site] = NEUTRAL_MOD
            else:
                site_to_bond[site] = None

        return Complex([Mol(self.name, site_to_mod, site_to_bond, False)])


class MolDefBuilder:
    def __init__(self, spec: Spec):
        self.name      = str(spec.to_non_struct_spec())  # type: str
        self.site_defs = {}                              # type: Dict[str, List[str]]

    def build(self) -> MolDef:
        return MolDef(self.name, self.site_defs)

    def add_site(self, site: Spec):
        if site_name(site) not in self.site_defs:
            self.site_defs[site_name(site)] = []

    def add_mod(self, site: Spec, mod: StateModifier):
        self.site_defs[site_name(site)].append(str(mod.value))


class Mol:
    def __init__(self, name: str, site_to_mod: Dict[str, str], site_to_bond: Dict[str, Optional[int]], is_reactant: bool):
        self.name         = name
        self.site_to_mod  = site_to_mod
        self.site_to_bond = site_to_bond
        self.is_reactant  = is_reactant
        self._assert_valid_site_names()

    def __str__(self) -> str:
        mod_str  = ','.join('{}{}'.format(site, '~' + mod)
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

    def __eq__(self, other: 'Mol') -> bool:
        # This equality skips checking the is_reactant property, since it does not appear
        # in the BNGL representation.
        return self.name  == other.name and self.site_to_mod == other.site_to_mod and self.site_to_bond == other.site_to_bond

    def __lt__(self, other: 'Mol') -> bool:
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
        return [s for s, b in self.site_to_bond.items() if b is not None and bond == b]

    def with_relabeled_bonds(self, bond_to_bond: Dict[int, int]) -> 'Mol':
        new_mol = self.clone()
        for site, old_bond in new_mol.site_to_bond.items():
            if old_bond is not None:
                new_mol.site_to_bond[site] = bond_to_bond[old_bond]

        return new_mol

    def _assert_valid_site_names(self):
        for site in self.sites:
            assert match(SITE_NAME_REGEX, site), 'Invalid site name: {}'.format(site)

def site_name(spec: Spec) -> str:
    bad_chars = ['[', ']', '/', '(', ')']
    spec_str = (spec.locus.domain + 'D' if spec.locus.domain else '') + \
               (spec.locus.subdomain + 'S' if spec.locus.subdomain else '') + \
               (spec.locus.residue + 'R' if spec.locus.residue else '')

    for bad_char in bad_chars:
        spec_str = spec_str.replace(bad_char, '')

    return spec_str


class MolBuilder:
    def __init__(self, spec: Spec, is_reactant: bool=False):
        self.name         = str(spec.to_non_struct_spec())
        self.site_to_mod  = {}  # type: Dict[str, str]
        self.site_to_bond = {}  # type: Dict[str, int]
        self.is_reactant  = is_reactant

    def build(self) -> Mol:
        return Mol(self.name, self.site_to_mod, self.site_to_bond, self.is_reactant)

    def set_bond_index(self, spec: Spec, bond_index: Optional[int]):
        self.site_to_bond[site_name(spec)] = bond_index

    def set_mod(self, spec: Spec, mod: Optional[StateModifier]):
        self.site_to_mod[site_name(spec)] = str(mod.value) if mod else None


class Complex:
    def __init__(self, mols: List[Mol]):
        self.mols = sorted(mols)

        self._assert_bonds_connected()
        self._assert_mols_connected()

    def __str__(self):
        return '.'.join(str(mol) for mol in self.mols)

    def __repr__(self):
        return 'Complex<{}>'.format(str(self))

    def __eq__(self, other: 'Complex') -> bool:
        return self.mols == other.mols

    def __lt__(self, other: 'Complex') -> bool:
        return self.mols < other.mols

    @property
    def bonds(self):
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
    def is_reactant(self):
        return any(mol.is_reactant for mol in self.mols)

    def _assert_mols_connected(self):
        connected = []
        to_visit  = [self.mols[0]]
        while to_visit:
            current = to_visit.pop()
            neighbors = self.neighbors(current)
            for neighbor in neighbors:
                if neighbor not in connected and neighbor not in to_visit:
                    to_visit.append(neighbor)

            connected.append(current)

        assert sorted(connected) == self.mols

    def _assert_bonds_connected(self):
        bond_to_site_count = defaultdict(int)

        for mol in self.mols:
            for bond in mol.bonds:
                bond_to_site_count[bond] += len(mol.sites_by_bond(bond))

        assert all(site_count == 2 for _, site_count in bond_to_site_count.items())


class ComplexExprBuilder:
    def __init__(self):
        self._mol_builders = {}  # type: Dict[Spec, MolBuilder]
        self._current_bond = 0
        self._bonds        = []  # type: List[Tuple[Spec, Spec]]

    def build(self, only_reactants: bool=True) -> List[Complex]:
        complexes = []

        for group in self._grouped_specs():
            possible_complex = Complex([self._mol_builders[spec].build() for spec in group])

            if not only_reactants or (only_reactants and possible_complex.is_reactant):
                complexes.append(possible_complex)

        return complexes

    def add_mol(self, spec: Spec, is_reactant: bool=False):
        if spec.to_component_spec() not in self._mol_builders:
            self._mol_builders[spec.to_component_spec()] = MolBuilder(spec.to_component_spec(), is_reactant)

    def set_bond(self, first: Spec, second: Spec):
        if (first.to_component_spec(), second.to_component_spec()) not in self._bonds and \
           (second.to_component_spec(), first.to_component_spec()) not in self._bonds:
            self._bonds.append((first.to_component_spec(), second.to_component_spec()))
            self._current_bond += 1
            self.set_half_bond(first, self._current_bond)
            self.set_half_bond(second, self._current_bond)

    def set_half_bond(self, spec: Spec, value: Optional[int]):
        self._mol_builders[spec.to_component_spec()].set_bond_index(spec, value)

    def set_mod(self, spec: Spec, mod: Optional[StateModifier]):
        self._mol_builders[spec.to_component_spec()].set_mod(spec, mod)

    def _grouped_specs(self) -> List[List[Spec]]:
        grouped_specs = [[spec] for spec in self._mol_builders.keys()]

        for bond in self._bonds:
            for i, group in enumerate(grouped_specs):
                if bond[0] in group:
                    grouped_specs.pop(i)
                    other_group = next((other_group for other_group in grouped_specs if bond[1] in other_group), [])
                    other_group += group
                    break

        return grouped_specs


STATE_TO_COMPLEX_BUILDER_FN = {
    # Covalent modification state.
    '$x-{$y}': [
        lambda state, builder: builder.add_mol(state['$x']),
        lambda state, builder: builder.set_mod(state['$x'], state['$y'])
    ],
    # Interaction state.
    '$x--$y': [
        lambda state, builder: builder.add_mol(state['$x']),
        lambda state, builder: builder.add_mol(state['$y']),
        lambda state, builder: builder.set_bond(state['$x'], state['$y'])
    ],
    # Self-interaction state.
    '$x--[$y]': [
        lambda state, builder: builder.add_mol(state['$x']),
        lambda state, builder: builder.set_bond(state.specs[0], state.specs[1])
    ],
    # Empty binding state.
    '$x--0': [
        lambda state, builder: builder.add_mol(state['$x']),
        lambda state, builder: builder.set_half_bond(state['$x'], None)
    ]
}


STATE_TO_MOL_DEF_BUILDER_FN = {
    # Covalent modification state.
    '$x-{$y}': [
        lambda state, builder: builder.add_site(state['$x']),
        lambda state, builder: builder.add_mod(state['$x'], state['$y'])
    ],
    # Interaction state.
    '$x--$y': [
        lambda state, builder: builder.add_site(state['$x']) if builder.name == state['$x'].component_name \
                               else builder.add_site(state['$y']),
    ],
    # Self-interaction state.
    '$x--[$y]': [
        lambda state, builder: builder.add_site(state.specs[0]),
        lambda state, builder: builder.add_site(state.specs[1])
    ],
    # Empty binding state.
    '$x--0': [
        lambda state, builder: builder.add_site(state['$x'])
    ]
}


class Parameter:
    def __init__(self, name: str, value: Optional[str]):
        self.name, self.value = name, value

    def __eq__(self, other: 'Parameter'):
        assert isinstance(other, Parameter)
        return self.name == other.name and self.value == other.value

    def __hash__(self):
        return hash(str(self))

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self)


class InitialCondition:
    def __init__(self, complex: Complex, value: Parameter):
        self.complex, self.value = complex, value

    def __str__(self):
        return '{} {}'.format(str(self.complex), str(self.value))

    def __repr__(self):
        return 'InitialCondition<{}>'.format(str(self))

    def is_equivalent_to(self, other: 'InitialCondition') -> bool:
        return self.complex.is_equivalent_to(other.complex) and self.value.name == other.value.name


class Observable:
    def __init__(self, name: str, complex: Complex):
        self.name, self.complex = name, complex

    def __str__(self):
        return '{} {}'.format(self.name, str(self.complex))

    def __repr__(self):
        return 'Observable<{}>'.format(str(self))


class Rule:
    def __init__(self, lhs: List[Complex], rhs: List[Complex], rate: Parameter, parent_reaction: Reaction=None):
        self.lhs, self.rhs, self.rate = sorted(lhs), sorted(rhs), rate
        self.parent_reaction = parent_reaction

    def __str__(self):
        return ' + '.join(str(x) for x in self.lhs) + ' -> ' + ' + '.join(str(x) for x in self.rhs) + \
               ' ' + str(self.rate) + ' ' + str(self.parent_reaction)

    def __repr__(self):
        return str(self)

    def __eq__(self, other: 'Rule') -> bool:
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

        for complex in my_lhs + my_rhs:
            complex.found = False

        for complex in other_lhs + other_rhs:
            complex.visited = False

        for my_complex in my_lhs:
            for other_complex in other_lhs:
                if my_complex.is_equivalent_to(other_complex) and not other_complex.visited:
                    my_complex.found = True
                    other_complex.visited = True
                    break

        for my_complex in my_rhs:
            for other_complex in other_rhs:
                if my_complex.is_equivalent_to(other_complex) and not other_complex.visited:
                    my_complex.found = True
                    other_complex.visited = True
                    break

        eq = all(complex.found for complex in my_lhs + my_rhs) and \
             all(complex.visited for complex in other_lhs and other_rhs)

        return eq


class RuleBasedModel:
    def __init__(self, mol_defs: List[MolDef], initial_conditions: List[InitialCondition], parameters: List[Parameter],
                 observables: List[Observable], rules: List[Rule]):
        self.mol_defs, self.initial_conditions, self.parameters, self.observables, self.rules = mol_defs, initial_conditions, \
            parameters, observables, rules

    @property
    def rate_parameters(self):
        return sorted(set(rule.rate for rule in self.rules), key=lambda x: x.name)


def calc_state_paths(states: List[State]) -> Dict[State, List[List[State]]]:
    specs = sorted(set(spec.to_component_spec() for state in states for spec in state.specs), key=lambda x: x.struct_index)

    def spec_to_bond_states(spec: Spec) -> List[State]:
        assert spec.is_component_spec and spec.is_structured
        return [state for state in states if spec in (s.to_component_spec() for s in state.specs) if len(state.specs) == 2]

    def neighbor(spec: Spec, state: State) -> Spec:
        assert spec.is_component_spec and spec.is_structured and len(state.specs) == 2
        return [neigh_spec.to_component_spec() for neigh_spec in state.specs if neigh_spec.to_component_spec() != spec][0]

    spec_paths = {spec: [] for spec in specs}  # type: Dict[Spec, List[List[State]]]

    for state in states:
        state.visited = []

    nums = count(1)

    def visit_nodes(current_path: List[State], current_spec: Spec, current_num: int):
        spec_paths[current_spec].append(current_path)
        bonds_to_visit = [state for state in spec_to_bond_states(current_spec) if current_num not in state.visited]
        for i, state in enumerate(bonds_to_visit):
            if i == 0:
                next_num = current_num
            else:
                next_num = next(nums)

            for s in current_path + [state]:
                s.visited.append(next_num)

            visit_nodes(current_path + [state], neighbor(current_spec, state), next_num)

    for spec in specs:
        if spec.struct_index > 1:
            break
        visit_nodes([], spec, 0)

    state_paths = {state: [] for state in states}  # type: Dict[State, List[List[State]]]
    for state in states:
        for spec in state.specs:
            for spec_path in spec_paths[spec.to_component_spec()]:
                if len(spec_path) > 0 and spec_path[-1] == state:
                    path = spec_path[:-1]
                else:
                    path = spec_path

                if path not in state_paths[state]:
                    state_paths[state].append(path)

    for state in states:
        del state.visited

    return state_paths


def with_connectivity_constraints(cont_set: VennSet) -> VennSet:
    dnf_terms = cont_set.to_dnf_list()
    connected_dnf_terms = []

    for dnf_term in dnf_terms:
        states = deepcopy([state for state in dnf_term.values if state])
        state_paths = calc_state_paths(states)

        constraint = UniversalSet()
        for state in states:
            if any(path == [] for path in state_paths[state]):
                continue

            state_constraints = [Complement(ValueSet(state))]
            for path in state_paths[state]:
                state_constraints.append(Intersection(*(ValueSet(x) for x in path)))

            constraint = Intersection(constraint, Union(*state_constraints))

        connected_dnf_terms.append(Intersection(dnf_term, constraint))

    return Union(*connected_dnf_terms)


def rule_based_model_from_rxncon(rxncon_sys: RxnConSystem) -> RuleBasedModel:
    def mol_defs_from_rxncon(rxncon_sys: RxnConSystem) -> Dict[Spec, MolDef]:
        mol_defs = {}
        for spec in rxncon_sys.components():
            builder = MolDefBuilder(spec)
            for state in rxncon_sys.states_for_component(spec):
                for func in STATE_TO_MOL_DEF_BUILDER_FN[state.repr_def]:
                    func(state, builder)

            mol_defs[spec] = builder.build()

        return mol_defs

    def venn_from_contingency(contingency: Contingency) -> VennSet:
        def parse_effector(eff: Effector) -> VennSet:
            if isinstance(eff, StateEffector):
                return ValueSet(eff.expr)
            elif isinstance(eff, NotEffector):
                return Complement(parse_effector(eff.expr))
            elif isinstance(eff, OrEffector):
                return Union(*(parse_effector(x) for x in eff.exprs))
            elif isinstance(eff, AndEffector):
                return Intersection(*(parse_effector(x) for x in eff.exprs))
            else:
                raise AssertionError

        if contingency.type in [ContingencyType.requirement]:
            return parse_effector(contingency.effector)
        elif contingency.type in [ContingencyType.inhibition]:
            return Complement(parse_effector(contingency.effector))
        else:
            return UniversalSet()

    def calc_positive_solutions(rxncon_sys: RxnConSystem, solution: Dict[State, bool]) -> List[List[State]]:
        def is_satisfiable(states: Iterable[State]) -> bool:
            for pair in combinations(states, 2):
                if pair[0].is_mutually_exclusive_with(pair[1]):
                    return False

            return True

        def complementary_state_combos(state: State) -> List[List[State]]:
            combos = product(
                *(rxncon_sys.complementary_states_for_component(spec.to_component_spec(), state) for spec in
                  state.specs))
            return [list(combo) for combo in combos if is_satisfiable(combo)]

        def structure_states(states: List[State]) -> List[State]:
            cur_index = max(spec.struct_index for state in states for spec in state.specs if spec.is_structured)

            spec_to_index = {}

            struct_states = []

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
                        state = state.to_structured_from_spec(spec.with_struct_index(cur_index))
                        cur_index += 1
                        spec_to_index[spec.to_component_spec()] = cur_index

                struct_states.append(state)

            return struct_states

        ordered_solution = OrderedDict(sorted(solution.items(), key=lambda x: x[0]))

        trues  = [state for state, val in ordered_solution.items() if val]
        falses = [state for state, val in ordered_solution.items() if not val
                  and not any(state.is_mutually_exclusive_with(x) for x in trues)]

        if not falses:
            return [trues] if is_satisfiable(trues) else []

        positivized_falses = [list(chain(*x)) for x in
                              product(*(complementary_state_combos(state) for state in falses))]

        solutions = []

        for positivized_false in positivized_falses:
            possible_solution = []
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

            states = copy(states)
            builder = ComplexExprBuilder()
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
                for func in STATE_TO_COMPLEX_BUILDER_FN[state.repr_def]:
                    func(state, builder)

            return builder.build()

        lhs = calc_complexes(reaction.terms_lhs, cont_soln)
        rhs = calc_complexes(reaction.terms_rhs, cont_soln)

        rate = Parameter('k', '1.0')

        return Rule(lhs, rhs, rate, parent_reaction=reaction)

    def calc_initial_conditions(mol_defs: List[MolDef]) -> List[InitialCondition]:
        return \
            [InitialCondition(mol_def.create_neutral_complex(),
                              Parameter('Num{}'.format(mol_def.name), str(INITIAL_MOLECULE_COUNT))) for mol_def in mol_defs]

    def calc_observables(rxncon_sys: RxnConSystem) -> List[Observable]:
        def observable_complex(states: List[State]) -> Complex:
            builder = ComplexExprBuilder()

            assert all(x.is_structured for x in states)

            for state in states:
                for func in STATE_TO_COMPLEX_BUILDER_FN[state.repr_def]:
                    func(state, builder)

            complexes = builder.build(only_reactants=False)

            assert len(complexes) == 1
            return complexes[0]

        observables = []
        output_rxns = [rxn for rxn in rxncon_sys.reactions if isinstance(rxn, OutputReaction)]
        for rxn in output_rxns:
            solns = Intersection(*(venn_from_contingency(x) for x
                                   in rxncon_sys.contingencies_for_reaction(rxn))).calc_solutions()
            positive_solns = []
            for soln in solns:
                positive_solns += calc_positive_solutions(rxncon_sys, soln)

            for index, positive_soln in enumerate(positive_solns):
                observables.append(Observable('{}{}'.format(rxn.name, index), observable_complex(positive_soln)))

        return observables

    mol_defs = list(mol_defs_from_rxncon(rxncon_sys).values())

    rules = []

    for reaction in (x for x in rxncon_sys.reactions if not isinstance(x, OutputReaction)):
        cont_set = Intersection(*(venn_from_contingency(x) for x
                                in rxncon_sys.s_contingencies_for_reaction(reaction)))

        cont_set = with_connectivity_constraints(cont_set)
        solutions = cont_set.calc_solutions()

        positive_solutions = []
        for solution in solutions:
            positive_solutions += calc_positive_solutions(rxncon_sys, solution)

        for positive_solution in positive_solutions:
            rule = calc_rule(reaction, positive_solution)
            if not any(rule.is_equivalent_to(existing) for existing in rules):
                rules.append(rule)

    return RuleBasedModel(mol_defs, calc_initial_conditions(mol_defs), [], calc_observables(rxncon_sys), rules)


def mol_from_str(mol_str: str) -> Mol:
    name, config_str = mol_str.split('(')
    assert config_str[-1] == ')'
    site_strs = config_str[:-1].split(',')

    site_to_mod  = {}
    site_to_bond = {}

    for site_str in site_strs:
        if not site_str:
            continue
        elif '~' not in site_str and '!' not in site_str:
            site_to_bond[site_str] = None
        elif '~' not in site_str:
            site, bond = site_str.split('!')
            site_to_bond[site] = bond
        elif '!' not in site_str:
            site, mod = site_str.split('~')
            site_to_mod[site] = mod
        else:
            raise AssertionError

    return Mol(name, site_to_mod, site_to_bond, False)


def complex_from_str(complex_str: str) -> Complex:
    return Complex([mol_from_str(x) for x in complex_str.split('.')])


def rule_from_str(rule_str: str) -> Rule:
    lhs_strs = []
    rhs_strs = []
    rate     = None

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

    return Rule([complex_from_str(x) for x in lhs_strs], [complex_from_str(x) for x in rhs_strs], rate)


def initial_condition_from_str(ic_str: str) -> InitialCondition:
    complex_str, param = ic_str.split()
    return InitialCondition(complex_from_str(complex_str), Parameter(param, None))
