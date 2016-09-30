from typing import Dict, List, Optional

from rxncon.core.state import State
from rxncon.core.spec import MolSpec

SiteName      = str
SiteModifier  = str
MolIndex      = int

class MolDef:
    def __init__(self, name: str, site_defs: Dict[SiteName, List[SiteModifier]]):
        self.name, self.site_defs = name, site_defs

    def valid_modifiers(self, site: SiteName):
        return self.site_defs[site]


class Mol:
    def __init__(self, mol_def: MolDef, site_modifiers: Dict[SiteName, SiteModifier]):
        self.mol_def, self.site_modifiers = mol_def, site_modifiers
        self.site_bonds = {site: None for site in self.site_modifiers.keys()}
        self._validate()

    def set_bond(self, site: SiteName, bond_num: int):
        assert not self.site_bonds[site]
        self.site_bonds[site] = bond_num

    def set_modifier(self, site: SiteName, modifier: SiteModifier):
        assert not self.site_modifiers[site]
        assert modifier in self.mol_def.valid_modifiers(site)
        self.site_modifiers[site] = modifier

    def _validate(self):
        for site, modifier in self.site_modifiers.items():
            assert modifier in self.mol_def.valid_modifiers(site)


class Complex:
    def __init__(self):
        self.mols       = {}  # type: Dict[MolIndex, Mol]
        self.bond_index = 0

    def apply_state(self, state: State):
        assert state.is_elemental


    def set_mol_at_index(self, mol: Mol, index: MolIndex):
        assert not self.mols[index]
        self.mols[index] = mol

    def set_bond(self, first_index: MolIndex, first_site: SiteName, second_index: MolIndex, second_site: SiteName):
        self.bond_index += 1
        self.mols[first_index].set_bond(first_site, self.bond_index)
        self.mols[second_index].set_bond(second_site, self.bond_index)


def str_from_spec(spec: MolSpec) -> str:
    bad_chars = ['[', ']', '/', '(', ')']
    spec_str = str(spec.to_non_struct_spec())
    for bad_char in bad_chars:
        spec_str = spec_str.replace(bad_char, '')

    return spec_str


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











