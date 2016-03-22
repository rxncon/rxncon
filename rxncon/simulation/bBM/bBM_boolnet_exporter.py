import rxncon.simulation.bBM.bipartite_boolean_model as bbm
import rxncon.simulation.bBM.bipartite_boolean_model as bbm
from typing import List, Union
import re
import rxncon.venntastic.sets as venn
import rxncon.core.reaction as rxn
import rxncon.core.state as sta
import rxncon.core.specification as spec


class BoolNet_System:
    def __init__(self, bipartite_boolean_model: bbm.Bipartite_Boolean_Model, ):
        self.bipartite_boolean_model = bipartite_boolean_model

    def to_string(self) -> str:
        bipartite_boolean_model_strings = [self._header_string(), self._build_init_conditions(), self._build_rules()]
        return "\n".join(bipartite_boolean_model_strings)

    def to_file(self, path: str):
        pass

    def _header_string(self):
        return "target, factors"

    def _build_rules(self):
        rules=[self._string_from_rule(rule) for rule in self.bipartite_boolean_model.rules]
        return "\n".join(rules)

    def _build_init_conditions(self):
        init_conditions = [self._string_from_init_cont(init_cont) for init_cont in self.bipartite_boolean_model.init_conditions]
        return "\n".join(init_conditions)

    def _string_from_init_cont(self, init_cont: bbm.InitConditions):
        if init_cont.value is None:
            return "{0}, {0}".format(self._generate_name(init_cont.target))
        else:
            raise NotImplementedError

    def _string_from_rule(self, rule: bbm.Rule):
        return "{0}, {1}".format(self._target_to_string(rule.target), self._factor_to_string(rule.factor.simplified_form().value))

    def _target_to_string(self, target: bbm.Node):
        return self._generate_name(target)

    def _factor_to_string(self, factor: bbm.Factor):
        #if factor == venn.PropertySet(venn.EmptySet()) or \
        #                factor == venn.PropertySet(venn.UniversalSet()):
        #    return ''
        if isinstance(factor, venn.PropertySet):
            return self._generate_name(factor.value)
        elif isinstance(factor, venn.Complement):
            return "! {0}".format(self._factor_to_string(factor.expr))
        elif isinstance(factor, venn.Intersection):
            return '({0} & {1})'.format(self._factor_to_string(factor.left_expr), self._factor_to_string(factor.right_expr))
        elif isinstance(factor, venn.Union):
            return '({0} | {1})'.format(self._factor_to_string(factor.left_expr), self._factor_to_string(factor.right_expr))
        else:
            raise NotImplementedError

    def _generate_name(self, node: bbm.Node):
        if isinstance(node.value, rxn.Reaction):
            return string_from_reaction(node.value)

        elif isinstance(node.value, sta.InterProteinInteractionState):
            return string_from_inter_protein_interaction_state(node.value)

        elif isinstance(node.value, sta.IntraProteinInteractionState):
            return string_from_intra_protein_interaction_state(node.value)

        elif isinstance(node.value, sta.CovalentModificationState):
            return string_from_covalent_modification_state(node.value)

        elif isinstance(node.value, sta.TranslocationState):
            return string_from_translocation_state(node.value)

        elif isinstance(node.value, sta.SynthesisDegradationState):
            return string_from_synthesis_degradation_state(node.value)

        elif isinstance(node.value, spec.Specification):
            return node.value


def string_from_specification(specification: spec.Specification):
    spec_str = str(specification)
    return replace_not_valid_signs(spec_str)


def string_from_reaction(reaction: rxn.Reaction) -> str:
    # todo: hieraus und aus folgenden eine tabelle fuer marcus machen
    def generate_reaction_verb_name(verb: rxn.Verb):
        if re.search("-", verb.value):
            return re.sub("-", "minus", verb.value)
        elif re.search("\+", verb.value):
            return re.sub("\+", "plus", verb.value)
        else:
            return verb.value
    reaction_str = "{0}_{1}_{2}".format(reaction.subject, generate_reaction_verb_name(reaction.verb), reaction.object)
    return replace_not_valid_signs(reaction_str)


def replace_not_valid_signs(value):
    value = re.sub("-", ".", value)
    value = re.sub('[\[{(]',"_", value )
    value = re.sub('[\]})]', "_", value )
    value = re.sub('\/', "_", value )
    return value


def string_from_inter_protein_interaction_state(state: rxn.sta.InterProteinInteractionState) -> str:
    return replace_not_valid_signs(str(state))


def string_from_intra_protein_interaction_state(state: rxn.sta.IntraProteinInteractionState) -> str:
    # A_[m]--[n]
    return replace_not_valid_signs(str(state))


def string_from_covalent_modification_state(state: rxn.sta.CovalentModificationState) -> str:
    # A-{P}
    return replace_not_valid_signs(str(state))



def string_from_translocation_state(state) -> str:
    # A-{cyto}
    return replace_not_valid_signs(str(state))



def string_from_synthesis_degradation_state(state) -> str:
    return replace_not_valid_signs(str(state))

