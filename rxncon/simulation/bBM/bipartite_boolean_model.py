import rxncon.venntastic.sets as venn
import rxncon.core.reaction as rxn
import rxncon.core.state as sta
import rxncon.core.specification as spec
import typecheck as tc
import typing as tg

class BipartiteBooleanModel:
    @tc.typecheck
    def __init__(self, rules: tg.List["Rule"], init_conditions: tg.List['InitCondition']):
        self.rules = rules
        self.init_conditions = init_conditions
        self._validate()

    def _validate(self):
        # check for all rules if all factors also have rules
        rule_targets = [rule.target for rule in self.rules]
        init_targets = [init_condition.target for init_condition in self.init_conditions]
        targets = rule_targets + init_targets
        for rule in self.rules:
            assert rule.factor.value is not None
            for list in rule.factor.value.to_nested_list_form():
                for property in list:
                    if hasattr(property, "expr"):
                        assert (property.expr.value in targets) is True
                    else:
                        assert (property.value in targets) is True


class InitCondition:
    @tc.typecheck
    def __init__(self, target: 'Node', value: tg.Optional[tg.Union[bool, float]]):
        self.target = target
        self.value = value

    @tc.typecheck
    def __eq__(self, other: 'InitCondition'):
        return self.target == other.target and self.value == self.value

    def __repr__(self):
        return str(self)
    def __str__(self):
        return "target: {0}, value: {1}".format(self.target, self.value)

class Node:
    @tc.typecheck
    def __init__(self, value: tg.Union[rxn.Reaction, sta.State]):
        self.value = value

    @tc.typecheck
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
    @tc.typecheck
    def __init__(self, target: Node, factor: 'Factor'):  # should factor expect class 'Factor' ?
        self.target = target
        self.factor = factor
        self._validate()

    def _validate(self):
        assert self.target.value is not None
        assert isinstance(self.target, Node)
        assert self.factor.value is not None
        assert isinstance(self.factor, Factor)
        #assert self.factor.value != venn.UniversalSet()
        #assert self.factor.value != venn.EmptySet

    def __str__(self):
        return "target: {0}, factors: {1}".format(self.target, self.factor)


class Factor(venn.Set):
    @tc.typecheck
    def __init__(self, factor: venn.Set):
        # venn.set e.g.:Union(Intersection(A--B, A_{P}), Intersection(A--C, C-{P})) --> (A--B & A_{P}) |(A--C & C_{P})
        self.value = factor
        self._validate()

    def _validate(self):
        #assert self.value != venn.UniversalSet()
        #assert self.value != venn.EmptySet
        pass