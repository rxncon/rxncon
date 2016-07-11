from enum import unique
from typecheck import typecheck

from rxncon.core.effector import Effector
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
    @typecheck
    def __init__(self, target: Reaction, type: ContingencyType, effector: Effector):
        self.target, self.type, self.effector = target, type, effector

    @typecheck
    def __eq__(self, other: 'Contingency') -> bool:
        return self.target == other.target and self.type == other.type and self.effector == other.effector

    @typecheck
    def __repr__(self) -> str:
        return str(self)

    @typecheck
    def __str__(self) -> str:
        return 'Contingency(target={0}, type={1}, effector={2}'.format(str(self.target),
                                                                       str(self.type),
                                                                       str(self.effector))

