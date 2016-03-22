import typing as tg
import rxncon.venntastic.sets as venn
import rxncon.core.reaction as rxn
import rxncon.core.state as sta
import rxncon.core.specification as spec

class Bipartite_Boolean_Model:
    def __init__(self, rules: tg.List["Rule"], init_conditions: tg.List['InitConditions']):
        self.rules = rules
        self.init_conditions = init_conditions
        self._validate()

    def _validate(self):
        pass

class InitConditions:
    def __init__(self, target: 'Node', value: tg.Optional[tg.Union[bool, 'Node']]):
        self.target = target
        self.value = value

    def __eq__(self, other: 'InitConditions'):
        return self.target == other.target and self.value == self.value

class Node:
    def __init__(self, value: tg.Union[rxn.Reaction, sta.State]):
        self.value = value

    def __eq__(self, other: 'Node'):
        if isinstance(self.value, rxn.Reaction) and isinstance(other.value, rxn.Reaction) and self.value == other.value:
            return True
        elif isinstance(self.value, sta.State) and isinstance(other.value, sta.State) and self.value == other.value:
            return True
        elif isinstance(self.value, spec.Specification) and isinstance(other.value, spec.Specification) and self.value == other.value:
            return True
        else:
            return False

    def __hash__(self):
        return hash(str(self))

    def __repr__(self):
        return str(self)

    def __str__(self):
        return str(self.value)


class Rule:
    def __init__(self, target: Node, factor: 'Factor'):  # should factor expect class 'Factor' ?
        self.target = target
        self.factor = factor
        self._validate()

    def _validate(self):
        pass


class Factor(venn.Set):
    def __init__(self, factor: venn.Set):
        # venn.set e.g.:Union(Intersection(A--B, A_{P}), Intersection(A--C, C-{P})) --> (A--B & A_{P}) |(A--C & C_{P})
        self.value = factor
        self._validate()

    def _validate(self):
        pass