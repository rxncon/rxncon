import rxncon.simulation.bBM.bipartite_boolean_model as bBm
from typing import List, Union
import re
import rxncon.venntastic.sets as venn
import rxncon.core.reaction as rxn
import rxncon.core.state as sta


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

    def _target_to_string(self, target: bBm.Node):
        return self._generate_name(target.value)

    def _factors_to_string(self, factor: bBm.Factor):
        nested_factor_list=factor.to_union_list_form()
        return self._nested_factor_list_to_string(nested_factor_list)

    def _nested_factor_list_to_string(self, nested_factor_list: List[List[venn.Set]]) -> str:
        result=[]
        for bool_and in nested_factor_list:
            result += self._bool_and_to_string(bool_and)

        return " | ".join(result)

    def _bool_and_to_string(self, bool_and: List[venn.Set]):
        elements=[self._generate_name(element) for element in bool_and]

        if len(elements) == 1:
            return elements[0]
        elif len(elements) > 1:
            return "({0})".format(" & ".join(elements))
        else:
            raise NotImplementedError


    def _generate_name(self, property_set: venn.PropertySet):
        if isinstance(property_set.value, rxn.Reaction):
            return string_from_reaction(property_set.value)

        elif isinstance(property_set.value, sta.InterProteinInteractionState):
            return string_from_inter_protein_interaction_state(property_set.value)

        elif isinstance(property_set.value, sta.IntraProteinInteractionState):
            return string_from_intra_protein_interaction_state(property_set.value)

        elif isinstance(property_set.value, sta.CovalentModificationState):
            return string_from_covalent_modification_state(property_set.value)

        elif isinstance(property_set.value, sta.TranslocationState):
            return string_from_translocation_state(property_set.value)

        elif isinstance(property_set.value, sta.SynthesisDegradationState):
            return string_from_synthesis_degradation_state(property_set.value)



        #re,sub("\W", "", self.name)  # replace everything not in [a-zA-Z0-9_] with nothing

def string_from_reaction(reaction_node: bBm.Node) -> str:
    def generate_reaction_verb_name(verb: rxn.Verb):
        if re.search("-", verb.value):
            return re.sub("-", "minus", verb.value)
        elif re.search("\+", verb.value):
            return re.sub("\+", "plus", verb.value)

    verb = generate_reaction_verb_name(reaction_node.value.verb)
    return "{0}_{1}_{2}".format(reaction_node.value.subject, generate_reaction_verb_name(reaction_node.value.verb), reaction_node.value.object)

def replace_breakets(value):
    value = re.sub('[\[{(]',"OpenBracket_", value )
    value = re.sub('[\]})]', "_CloseBracket", value )
    return value

def string_from_inter_protein_interaction_state(state) -> str:
    result =  '{0}__{1}'.format(state.first_component, state.second_component)
    return replace_breakets(result)


def string_from_intra_protein_interaction_state(state) -> str:
    # A_[m]--[n]
    result = '{0}__OpenBracket_{1}_CloseBracket'.format(state.first_component, state.second_component.domain)
    return replace_breakets(result)

def string_from_covalent_modification_state(state) -> str:
    # A-{P}
    result = '{0}_{1}'.format(state.substrate, state.modifier.value)
    return replace_breakets(result)


def string_from_translocation_state(state) -> str:
    # A-{cyto}
    result = '{0}_{1}'.format(state.substrate, state.compartment.value)
    return replace_breakets(result)


def string_from_synthesis_degradation_state(state) -> str:
    result = '{}'.format(state.component)
    return replace_breakets(result)
