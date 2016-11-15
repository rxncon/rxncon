from typing import Dict, List, Optional, Tuple

from rxncon.core.state import State
from rxncon.core.spec import Spec


class MolDef:
    def __init__(self, name: str, site_defs: Dict[str, List[str]]):
        self.name, self.site_defs = name, site_defs

    def values_for_site(self, site: str):
        return self.site_defs[site]

    @property
    def sites(self):
        return self.site_defs.keys()


class Mol:
    def __init__(self, mol_def: MolDef, site_to_value: Dict[str, str]):
        self.mol_def = mol_def
        self.sites   = mol_def.site_defs.keys()

        self.site_to_value = {site: None for site in self.sites}
        for site, modifier in site_to_value.items():
            self.site_to_value[site] = modifier

        self._validate()

    @property
    def name(self):
        return self.mol_def.name

    def _validate(self):
        for site, state in self.site_to_value.items():
            assert state in self.mol_def.values_for_site(site)


def str_from_spec(spec: Spec) -> str:
    bad_chars = ['[', ']', '/', '(', ')']
    spec_str = str(spec.to_non_struct_spec())
    for bad_char in bad_chars:
        spec_str = spec_str.replace(bad_char, '')

    return spec_str


class BuildingArea:
    pass


STATE_TO_THING = {
    # Covalent modification state.
    '$x-{$y}': [
        lambda state, building: building.ensure_mol(state, '$x'),
        lambda state, building: building.set_site_to_value(state, '$x', site_value(state))
    ],
    # Interaction state.
    '$x--$y': [
        lambda state, building: ensure_mol(state, building, '$x'),
        lambda state, building: ensure_mol(state, building, '$y'),
        lambda state, building: set_bond(state, building, '$x', '$y')
    ],
    # Self-interaction state.
    '$x--[$y]': [
        lambda state, building: ensure_mol(state, building, '$x'),
        lambda state, building: set_bond(state, building, '$x', '$x.to_component_spec_[$y]')
    ],
    # Empty binding state.
    '$x--0': [
        lambda state, building: ensure_mol(state, building, '$x'),
        lambda state, building: set_site_to_state()
    ]

}


class Complex:
    def __init__(self, mols: List[Mol], bonds: List[Tuple[Spec, Spec]]):
        self.mols  = mols
        self.bonds = bonds
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











