from typing import List, Dict, Tuple

from rxncon.core.spec import MolSpec, LocusResolution
from rxncon.core.state import State, StateDef
from rxncon.core.rxncon_system import RxnConSystem
from rxncon.util.utils import all_eq


class MutualExclusivityError(Exception):
    pass


class MolDef:
    def __init__(self, spec: MolSpec, grouped_states: Dict[Tuple[MolSpec, StateDef], List[State]]):
        assert spec.has_resolution(LocusResolution.component)
        self.spec = spec
        self.grouped_states = grouped_states
        self._validate()

    def __str__(self) -> str:
        return 'MolDef<{0}>, states: <{1}>'.format(str(self.spec), self.grouped_states)

    def __repr__(self) -> str:
        return str(self)

    def _validate(self):
        for (spec, state_def), states in self.grouped_states.items():
            assert all(spec in state.mol_specs for state in states)
            assert all_eq([state_def] + [state.definition for state in states])



def mol_def_for_component(rxncon_sys: RxnConSystem, component: MolSpec) -> MolDef:
    return MolDef(component, rxncon_sys.states_for_component_grouped(component))
