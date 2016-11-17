from typing import Dict, List, Optional, Tuple, Iterable
from itertools import combinations, product, chain
from copy import copy, deepcopy

from rxncon.core.rxncon_system import RxnConSystem
from rxncon.core.reaction import Reaction, ReactionTerm
from rxncon.core.state import State, StateModifier
from rxncon.core.spec import Spec
from rxncon.core.contingency import Contingency, ContingencyType
from rxncon.core.effector import Effector, AndEffector, OrEffector, NotEffector, StateEffector
from rxncon.venntastic.sets import Set as VennSet, Intersection, Union, Complement, ValueSet, UniversalSet


class MolDef:
    def __init__(self, name: str, site_defs: Dict[str, List[str]]):
        self.name, self.site_defs = name, site_defs

    def __str__(self):
        site_strs = []
        for site_name, site_def in self.site_defs.items():
            site_strs.append('{0}:{1}'.format(site_name, '~'.join(x for x in site_def))) if site_def else site_strs.append(site_name)
        return '{0}({1})'.format(self.name, ','.join(site_strs))

    def __repr__(self):
        return str(self)

    def mods_for_site(self, site: str):
        return self.site_defs[site]

    @property
    def sites(self):
        return self.site_defs.keys()


class MolDefBuilder:
    def __init__(self, name: str):
        self.name      = name  # type: str
        self.site_defs = {}    # type: Dict[str, List[str]]

    def build(self):
        return MolDef(self.name, self.site_defs)

    def add_site(self, site: Spec):
        if site_name(site) not in self.site_defs:
            self.site_defs[site_name(site)] = []

    def add_mod(self, site: Spec, mod: StateModifier):
        self.site_defs[site_name(site)].append(str(mod.value))


class Mol:
    def __init__(self, spec: Spec, site_to_mod: Dict[str, Optional[str]], site_to_bond: Dict[str, int], is_reactant: bool):
        assert spec.is_component_spec
        self.spec         = spec
        self.site_to_mod  = site_to_mod
        self.site_to_bond = site_to_bond
        self.is_reactant  = is_reactant

    def __str__(self):
        mod_str  = ','.join('{}{}'.format(site, '~' + mod if mod else '') for site, mod in self.site_to_mod.items())
        bond_str = ','.join('{}{}'.format(site, '!' + str(bond) if bond else '') for site, bond in self.site_to_bond.items())

        strs = []
        if mod_str:
            strs.append(mod_str)
        if bond_str:
            strs.append(bond_str)

        return '{}({})'.format(self.spec.component_name, ','.join(strs))

    def __repr__(self):
        return str(self)

    @property
    def name(self):
        return self.spec.component_name


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
        self.spec         = spec
        self.site_to_mod  = {}  # type: Dict[str, str]
        self.site_to_bond = {}  # type: Dict[str, int]
        self.is_reactant  = is_reactant

    def build(self) -> Mol:
        return Mol(self.spec, self.site_to_mod, self.site_to_bond, self.is_reactant)

    def set_bond_index(self, spec: Spec, bond_index: Optional[int]):
        self.site_to_bond[site_name(spec)] = bond_index

    def set_mod(self, spec: Spec, mod: Optional[StateModifier]):
        self.site_to_mod[site_name(spec)] = str(mod.value) if mod else None


class Complex:
    def __init__(self, mols: List[Mol]):
        self.mols  = mols

    def __str__(self):
        return '.'.join(str(mol) for mol in self.mols)

    @property
    def is_reactant(self):
        return any(mol.is_reactant for mol in self.mols)


class ComplexExprBuilder:
    def __init__(self):
        self._mol_builders = {}  # type: Dict[Spec, MolBuilder]
        self._current_bond = 0
        self._bonds        = []  # type: List[Tuple[Spec, Spec]]

    def build(self) -> List[Complex]:
        complexes = []

        for group in self._grouped_specs():
            possible_complex = Complex([self._mol_builders[spec].build() for spec in group])

            if possible_complex.is_reactant:
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
    def __init__(self, name: Optional[str], value: Optional[str]):
        assert name or value
        self.name, self.value = name, value

    def __str__(self):
        if not self.value:
            return self.name
        elif not self.name:
            return self.value
        else:
            return '{} = {}'.format(self.name, self.value)

    def __repr__(self):
        return str(self)

class InitialCondition:
    def __init__(self, complex: Complex, value: Parameter):
        self.complex, self.value = complex, value


class Observable:
    def __init__(self, name: str, complex: Complex):
        self.name, self.complex = name, complex


class Rule:
    def __init__(self, lhs: List[Complex], rhs: List[Complex], rate: Parameter, parent_reaction: Reaction=None):
        self.lhs, self.rhs, self.rate = lhs, rhs, rate
        self.parent_reaction = parent_reaction

    def __str__(self):
        return ' + '.join(str(x) for x in self.lhs) + ' -> ' + ' + '.join(str(x) for x in self.rhs) + \
               ' ' + str(self.rate) + ' ' + str(self.parent_reaction)

    def __repr__(self):
        return str(self)


class RuleBasedModel:
    def __init__(self, mol_defs: List[MolDef], initial_conditions: List[InitialCondition], parameters: List[Parameter],
                 observables: List[Observable], rules: List[Rule]):
        self.mol_defs, self.initial_conditions, self.parameters, self.observables, self.rules = mol_defs, initial_conditions, \
            parameters, observables, rules


def calc_positive_solutions(rxncon_sys: RxnConSystem, solution: Dict[State, bool]) -> List[List[State]]:
    def is_satisfiable(states: Iterable[State]) -> bool:
        for pair in combinations(states, 2):
            if pair[0].is_mutually_exclusive_with(pair[1]):
                return False

        return True

    def complementary_state_combos(state: State) -> List[List[State]]:
        combos = product(*(rxncon_sys.complementary_states_for_component(spec.to_component_spec(), state) for spec in state.specs))
        return [list(combo) for combo in combos if is_satisfiable(combo)]

    trues  = [state for state, val in solution.items() if val]
    falses = [state for state, val in solution.items() if not val]

    if not falses:
        return [trues]

    positivized_falses = [list(chain(*x)) for x in product(*(complementary_state_combos(state) for state in falses))]

    solutions = []

    for positivized_false in positivized_falses:
        possible_solution = trues + positivized_false
        if is_satisfiable(possible_solution):
            solutions.append(possible_solution)

    return solutions


def rule_based_model_from_rxncon(rxncon_sys: RxnConSystem) -> RuleBasedModel:
    def mol_defs_from_rxncon(rxncon_sys: RxnConSystem) -> Dict[Spec, MolDef]:
        mol_defs = {}
        for spec in rxncon_sys.components():
            builder = MolDefBuilder(spec.component_name)
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
                return Union(parse_effector(eff.left_expr), parse_effector(eff.right_expr))
            elif isinstance(eff, AndEffector):
                return Intersection(parse_effector(eff.left_expr), parse_effector(eff.right_expr))
            else:
                raise AssertionError

        if contingency.type in [ContingencyType.requirement, ContingencyType.positive]:
            return parse_effector(contingency.effector)
        elif contingency.type in [ContingencyType.inhibition, ContingencyType.negative]:
            return Complement(parse_effector(contingency.effector))
        else:
            return UniversalSet()

    def calc_rule(reaction: Reaction, cont_soln: List[State]) -> Rule:
        def calc_complexes(terms: List[ReactionTerm], states: List[State]) -> List[Complex]:
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

            for state in states:
                for func in STATE_TO_COMPLEX_BUILDER_FN[state.repr_def]:
                    func(state, builder)

            return builder.build()

        lhs = calc_complexes(reaction.terms_lhs, cont_soln)
        rhs = calc_complexes(reaction.terms_rhs, cont_soln)

        return Rule(lhs, rhs, Parameter('k', None), parent_reaction=reaction)

    mol_defs = list(mol_defs_from_rxncon(rxncon_sys).values())

    rules = []

    for reaction in rxncon_sys.reactions:
        solutions = Intersection(*[venn_from_contingency(x) for x
                                   in rxncon_sys.contingencies_for_reaction(reaction)]).calc_solutions()

        positive_solutions = []
        for solution in solutions:
            positive_solutions += calc_positive_solutions(rxncon_sys, solution)

        for positive_solution in positive_solutions:
            rules.append(calc_rule(reaction, positive_solution))

    return RuleBasedModel(mol_defs, [], [], [], rules)
