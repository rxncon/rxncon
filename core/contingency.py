from enum import Enum, unique

import core.effector as eff
import core.reaction as rxn


@unique
class ContingencyType(Enum):
    undefined   = None
    requirement = '!'
    inhibition  = 'x'
    positive    = 'K+'
    negative    = 'K-'
    no_effect   = '0'
    unknown     = '?'


class Contingency:
    def __init__(self, target: rxn.Reaction, type: ContingencyType, effector: eff.Effector):
        self.target = target
        self.type = type
        self.effector = effector

    def __eq__(self, other: 'Contingency') -> bool:
        assert isinstance(other, Contingency)
        return self.target == other.target and self.type == other.type and self.effector == other.effector

    def __str__(self) -> str:
        return 'Contingency(target={0}, type={1}, effector={2}'.format(self.target, self.type, self.effector)

