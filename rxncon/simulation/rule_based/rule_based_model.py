from typing import Dict, List, Optional, Tuple

from rxncon.core.rxncon_system import RxnConSystem
from rxncon.core.state import State, StateModifier
from rxncon.core.spec import Spec


class MolDef:
    def __init__(self, name: str, site_defs: Dict[str, List[str]]):
        self.name, self.site_defs = name, site_defs

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

    def add_site(self, site: str):
        if site not in self.site_defs:
            self.site_defs[site] = []

    def add_mod_for_site(self, site: str, mod: str):
        self.site_defs[site] += mod


class Mol:
    def __init__(self, spec: Spec, site_to_mod: Dict[str, str], site_to_bond: Dict[str, int]):
        assert spec.is_component_spec
        self.spec         = spec
        self.site_to_mod  = site_to_mod
        self.site_to_bond = site_to_bond

    @property
    def name(self):
        return self.spec.component_name


def site_name(spec: Spec) -> str:
    bad_chars = ['[', ']', '/', '(', ')']
    spec_str = str(spec.locus)
    for bad_char in bad_chars:
        spec_str = spec_str.replace(bad_char, '')

    return spec_str


class MolBuilder:
    def __init__(self, spec: Spec):
        self.spec         = spec
        self.site_to_mod  = {}  # type: Dict[str, str]
        self.site_to_bond = {}  # type: Dict[str, int]

    def build(self) -> Mol:
        return Mol(self.spec, self.site_to_mod, self.site_to_bond)

    def set_bond_index(self, spec: Spec, bond_index: int):
        self.site_to_bond[site_name(spec)] = bond_index

    def set_mod(self, spec: Spec, mod: Optional[StateModifier]):
        self.site_to_mod[site_name(spec)] = str(mod) if mod else None


class ComplexExprBuilder:
    def __init__(self):
        self._mol_builders = {}  # type: Dict[Spec, MolBuilder]
        self._current_bond = 0
        self._bonds        = []  # type: List[Tuple[Spec, Spec]]

    def build(self) -> List[Complex]:
        complexes = []

        for group in self._grouped_specs():
            complexes.append(Complex([self._mol_builders[spec].build() for spec in group]))

        return complexes

    def add_mol(self, spec: Spec):
        if spec.to_component_spec() not in self._mol_builders:
            self._mol_builders[spec.to_component_spec()] = MolBuilder(spec.to_component_spec())

    def set_bond(self, first: Spec, second: Spec):
        self._bonds.append((first.to_component_spec(), second.to_component_spec()))
        self._current_bond += 1
        self._mol_builders[first.to_component_spec()].set_bond_index(first, self._current_bond)
        self._mol_builders[second.to_component_spec()].set_bond_index(second, self._current_bond)

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
        lambda state, builder: builder.set_mod(state['$x'], None)
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



class Complex:
    def __init__(self, mols: List[Mol]):
        self.mols  = mols
        self._validate()

    def _validate(self):
        pass


class Parameter:
    def __init__(self, name: Optional[str], value: Optional[str]):
        assert name or value
        self.name, self.value = name, value


class InitialCondition:
    def __init__(self, complex: Complex, value: Parameter):
        self.complex, self.value = complex, value


class Observable:
    def __init__(self, name: str, complex: Complex):
        self.name, self.complex = name, complex


class Rule:
    def __init__(self, lhs: List[Complex], rhs: List[Complex], rate: Parameter):
        self.lhs, self.rhs, self.rate = lhs, rhs, rate


class RuleBasedModel:
    def __init__(self, mol_defs: List[MolDef], initial_conditions: List[InitialCondition], parameters: List[Parameter],
                 observables: List[Observable], rules: List[Rule]):
        self.mol_defs, self.initial_conditions, self.parameters, self.observables, self.rules = mol_defs, initial_conditions, \
            parameters, observables, rules


def rule_based_model_from_rxncon(rxncon_sys: RxnConSystem) -> RuleBasedModel:
    def mol_defs_from_rxncon(rxncon_sys: RxnConSystem) -> Dict[Spec, MolDef]:
        mol_defs = {}
        for spec in rxncon_sys.components():
            builder = MolDefBuilder(spec.component_name)
            for state in rxncon_sys.states_for_component(spec):
                STATE_TO_MOL_DEF_BUILDER_FN[state.repr_def](state, builder)

            mol_defs[spec] = builder.build()

        return mol_defs

    mol_defs = mol_defs_from_rxncon(rxncon_sys)



