import rxncon.simulation.bBM.bipartite_boolean_model as bBm
from typing import List
import rxncon.venntastic.sets as venn



class bBm_System:
    def __init__(self, bipartite_boolean_model: bBm.Bipartite_Boolean_Model):
        self.bipartite_boolean_model = bipartite_boolean_model

    def to_string(self) -> str:
        bipartite_boolean_model_strings = [self._header_string(), self._build_rules()]
        return "\n".join(bipartite_boolean_model_strings)

    def to_file(self, path: str):
        pass

    def _header_string(self):
        return "target, factors"

    def _build_rules(self):
        rules=[self._rule_to_string(rule) for rule in self.bipartite_boolean_model.rules]
        return "\n".join(rules)

    def _rule_to_string(self, rule: bBm.Rule):
        return "{0}, {1}".format(self._target_to_string(rule.target), self._factors_to_string(rule.factor))

    def _target_to_string(self, target: bBm.Target):
        return target.name

    def _factors_to_string(self, factor: bBm.Factor):
        nested_factor_list=factor.to_union_list_form()
        return self._nested_factor_list_to_string(nested_factor_list)

    def _nested_factor_list_to_string(self, nested_factor_list: List[List[venn.Set]]) -> str:
        result=[]
        for bool_and in nested_factor_list:
            result += self._bool_and_to_string(bool_and)

        return " | ".join(result)

    def _bool_and_to_string(self, bool_and: List[venn.Set]):
        elements=[element.name for element in bool_and]

        if len(elements) == 1:
            return elements[0]
        elif len(elements) > 1:
            return "({0})".format(" & ".join(elements))
        else:
            raise NotImplementedError
