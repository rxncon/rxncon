from collections import defaultdict

from rxncon.core.specification import Specification, SpecificationResolution
from rxncon.core.state import State


class MutualExclusivityError(Exception):
    pass


class MoleculeDef:
    def __init__(self, spec: Specification):
        assert spec.has_resolution(SpecificationResolution.component)
        self.spec = spec
        self.props_def = defaultdict(set)

    def add_state(self, state: State):
        for spec in state.specs:
            if spec.to_component_specification() == self.spec:
                self.props_def[spec].add(state)

    def states_for_spec(self, spec: Specification):
        return self.props_def[spec]

    @property
    def states(self):
        res = set()
        for spec, states in self.props_def.items():
            res = res.union(states)

        return res

    @property
    def specs(self):
        return self.props_def.keys()


class Molecule:
    def __init__(self, molecule_def: MoleculeDef, structure_index: int):
        self.molecule_def = molecule_def
        self.structure_index = structure_index
        self.props = {}

    def __eq__(self, other):
        return self.molecule_def == other.molecule_def and self.structure_index == other.structure_index and \
            self.props == other.props

    @property
    def bonds(self):
        return [Bond(*state.specs) for state in self.states if state.is_bond]

    @property
    def states(self):
        return self.props.values()

    def component_matches(self, other: 'Molecule'):
        return self.molecule_def == other.molecule_def and self.structure_index == other.structure_index

    def merge_with(self, other: 'Molecule'):
        assert self.component_matches(other)
        for state in other.states:
            self.set_state(state)

    def set_state(self, state: State):
        for spec in state.specs:
            if spec.to_component_specification() == self.molecule_def.spec:
                if spec in self.props and self.props[spec] != state:
                    raise MutualExclusivityError
                self.props[spec] = state


class Bond:
    def __init__(self, left_spec: Specification, right_spec: Specification):
        self.left_spec, self.right_spec = sorted([left_spec, right_spec])

    def __eq__(self, other: 'Bond'):
        return self.left_spec == other.left_spec and self.right_spec == other.right_spec





