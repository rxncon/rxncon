from typing import List, Union
import rxncon.venntastic.sets as venn
import rxncon.core.reaction as rxn
import rxncon.core.state as sta


class Bipartite_Boolean_Model:
    def __init__(self, rules: List["Rule"]):
        self.rules = rules
        self._validate()

    def _validate(self):
        pass


class Node:
    def __init__(self, value: Union[rxn.Reaction, sta.State]):
        self.value = value

    def __eq__(self, other: 'Node'):
        if isinstance(self.value, rxn.Reaction) and isinstance(other.value, rxn.Reaction) and self.value == other.value:
                return True
        elif isinstance(self.value, sta.State) and isinstance(other.value, sta.State) and self.value == other.value:
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
    def __init__(self, target: Node, factor: "Factor"):
        self.target = target
        self.factor = factor
        self._validate()

    def _validate(self):
        pass


class Factor:
    def __init__(self, factor: venn.Set ):
        # venn.set e.g.:Union(Intersection(A--B, A_{P}), Intersection(A--C, C-{P})) --> (A--B & A_{P}) |(A--C & C_{P})
        self.factor = factor
        self._validate()

    def _validate(self):
        pass