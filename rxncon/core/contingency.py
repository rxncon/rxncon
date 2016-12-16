from enum import unique
from copy import deepcopy

from rxncon.core.effector import Effector, StructEquivalences, QualSpec, StateEffector
from rxncon.core.reaction import Reaction
from rxncon.util.utils import OrderedEnum


@unique
class ContingencyType(OrderedEnum):
    undefined   = None
    requirement = '!'
    inhibition  = 'x'
    positive    = 'k+'
    negative    = 'k-'
    no_effect   = '0'
    unknown     = '?'


class Contingency:
    def __init__(self, target: Reaction, type: ContingencyType, effector: Effector):
        self.target, self.type, self.effector = target, type, effector

    def __eq__(self, other: 'Contingency') -> bool:
        return self.target == other.target and self.type == other.type and self.effector == other.effector

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'Contingency({0}, {1}, {2}'.format(str(self.target), str(self.type), str(self.effector))

    def clone(self) -> 'Contingency':
        return deepcopy(self)

    def to_structured(self):
        if isinstance(self.effector, StateEffector) and not self.effector.is_structured:
            sc = self.clone()

            struct_components = {spec: spec.with_struct_index(index) for index, spec in enumerate(self.target.components_lhs)}
            updates = {}
            for spec in self.effector.expr.specs:
                try:
                    updates[spec] = spec.with_struct_from_spec(struct_components[spec.to_component_spec()])
                except KeyError:
                    pass

            sc.effector.expr.update_specs(updates)
            return sc
        elif isinstance(self.effector, StateEffector) and self.effector.is_structured:
            return self



        sc = self.clone()
        sc.effector = sc.effector.to_merged_struct_effector()
        return sc

