from enum import Enum, unique
import typecheck as tc

import rxncon.core.effector as eff
import rxncon.core.reaction as rxn


@unique
class ContingencyType(Enum):
    undefined   = None
    requirement = '!'
    inhibition  = 'x'
    positive    = 'k+'
    negative    = 'k-'
    no_effect   = '0'
    unknown     = '?'


class Contingency:
    @tc.typecheck
    def __init__(self, target: rxn.Reaction, type: ContingencyType, effector: eff.Effector):
        self.target = target
        self.type = type
        self.effector = effector

    @tc.typecheck
    def __eq__(self, other: 'Contingency') -> bool:
        return self.target == other.target and self.type == other.type and self.effector == other.effector

    def __str__(self) -> str:
        return 'Contingency(target={0}, type={1}, effector={2}'.format(self.target, self.type, self.effector)

