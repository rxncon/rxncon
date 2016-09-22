from typing import List, Dict, Tuple
from collections import defaultdict
from typecheck import typecheck

from rxncon.core.spec import MolSpec, LocusResolution
from rxncon.core.state import State, StateDef
from rxncon.core.rxncon_system import RxnConSystem
from rxncon.util.utils import all_eq


class MutualExclusivityError(Exception):
    pass


class MolDef:
    @typecheck
    def __init__(self, component: MolSpec, valid_states: Dict[Tuple[MolSpec, StateDef], List[State]]):
        assert component.has_resolution(LocusResolution.component)
        assert not component.struct_index

        self.component    = component
        self.valid_states = valid_states
        self._valid_states_nested = defaultdict(dict)
        for (spec, state_def), states in self.valid_states.items():
            self._valid_states_nested[spec][state_def] = states

        self._validate()

    def __str__(self) -> str:
        return 'MolDef<{0}>, states: <{1}>'.format(str(self.component), self.valid_states)

    def __repr__(self) -> str:
        return str(self)

    @typecheck
    def valid_states_by_spec(self, spec: MolSpec) -> Dict[StateDef, List[State]]:
        return self._valid_states_nested[spec]

    def _validate(self):
        for (spec, state_def), states in self.valid_states.items():
            assert all(spec in state.mol_specs for state in states)
            assert all_eq([state_def] + [state.definition for state in states])
            assert not spec.struct_index
            assert not any(state.is_struct_state for state in states)


class Mol:
    @typecheck
    def __init__(self, mol_def: MolDef, states: Dict[Tuple[MolSpec, StateDef], State]):
        self.mol_def      = mol_def
        self.states       = states
        self.struct_index = None

        self._determine_struct_index()
        self._validate()

    @property
    def component(self) -> MolSpec:
        return self.mol_def.component

    def _determine_struct_index(self):
        possible_indices = {'all'}

        for state in self.states.values():
            indices = {component.struct_index for component in state.components if component.is_non_struct_equiv_to(self.component)
                       and component.struct_index is not None}

            if indices and possible_indices == {'all'}:
                possible_indices = set(indices)
            elif indices:
                possible_indices = possible_indices & indices

        if possible_indices == set():
            raise Exception('Inconsistent structure indices on Mol {}'.format(str(self)))
        elif len(possible_indices) > 1:
            raise Exception('Insufficient structure annotation on Mol {}'.format(str(self)))
        elif len(possible_indices) == 1:
            self.struct_index = list(possible_indices)[0]

    def _validate(self):
        for (spec, state_def), state in self.states:
            assert state.to_non_struct_state() in self.mol_def.valid_states[(spec, state_def)]


@typecheck
def mol_def_for_component(rxncon_sys: RxnConSystem, component: MolSpec) -> MolDef:
    return MolDef(component, rxncon_sys.states_for_component_grouped(component))


# class Complex:
#     def __init__(self, mols: List[Mol], ):
#


