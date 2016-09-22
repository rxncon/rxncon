from typing import List, Dict, Tuple
from collections import defaultdict
from typecheck import typecheck
from copy import deepcopy

from rxncon.core.spec import MolSpec, LocusResolution
from rxncon.core.state import State, StateDef
from rxncon.core.rxncon_system import RxnConSystem
from rxncon.util.utils import all_eq


class MutualExclusivityError(Exception):
    pass


class InconsistentStructureError(Exception):
    pass


class InsufficientStructureError(Exception):
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
        return 'MolDef<{0}>, valid_states: <{1}>'.format(str(self.component), self.valid_states)

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
        self.mol_def   = mol_def
        self.states    = states
        self.component = deepcopy(mol_def.component)
        self._mutation_post_process()

    def __str__(self):
        return 'Mol<{0}>, states: <{1}>'.format(str(self.mol_def.component), str(self.states))

    def __repr__(self):
        return str(self)

    def add_state(self, mol_spec: MolSpec, state: State):
        assert not mol_spec.struct_index

        if (mol_spec, state.definition) in self.states.keys():
            raise MutualExclusivityError

        self.states[(mol_spec, state.definition)] = state
        self._mutation_post_process()

    def _mutation_post_process(self):
        self._determine_struct_index()
        self._validate()

    def _determine_struct_index(self):
        UNIVERSAL_SET = {'all'}

        possible_indices = UNIVERSAL_SET

        for state in self.states.values():
            indices = {component.struct_index for component in state.components if component.is_non_struct_equiv_to(self.component)
                       and component.struct_index is not None}

            if indices and possible_indices == UNIVERSAL_SET:
                possible_indices = set(indices)
            elif indices:
                possible_indices = possible_indices & indices

        if possible_indices == set():
            raise InconsistentStructureError('Inconsistent structure annotation on Mol {}'.format(str(self)))
        elif len(possible_indices) > 1:
            raise InsufficientStructureError('Insufficient structure annotation on Mol {}'.format(str(self)))
        elif possible_indices != UNIVERSAL_SET and len(possible_indices) == 1:
            self.component.struct_index = list(possible_indices)[0]

    def _validate(self):
        for (spec, state_def), state in self.states.items():
            assert state.to_non_struct_state() in self.mol_def.valid_states[(spec, state_def)]


@typecheck
def mol_def_for_component(rxncon_sys: RxnConSystem, component: MolSpec) -> MolDef:
    return MolDef(component, rxncon_sys.states_for_component_grouped(component))


# class Complex:
#     def __init__(self, mols: List[Mol], ):
#


