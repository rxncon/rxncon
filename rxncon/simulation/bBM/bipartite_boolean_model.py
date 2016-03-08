from typing import List
import rxncon.venntastic.sets as venn

class Bipartite_Boolean_Model:
    def __init__(self, rules: List["Rule"]):
        self.rules = rules
        self._validate()

    def _validate(self):
        pass

class Rule:
    def __init__(self, target: "Target", factor: "Factor"):
        self.target = target
        self.factor = factor
        self._validate()

    def _validate(self):
        pass

class Target:
    def __init__(self, name: str ):
        self.name=name
        self._validate()

    def _validate(self):
        pass

class Factor:
    def __init__(self, factor: venn.Set ):
        # venn.set e.g.:Union(Intersection(A--B, A_{P}), Intersection(A--C, C-{P})) --> (A--B & A_{P}) |(A--C & C_{P})
        self.factor = factor
        self.validate()


    def _validate(self):
        pass