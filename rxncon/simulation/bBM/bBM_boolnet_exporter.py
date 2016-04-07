import re
import os

import rxncon.simulation.bBM.bipartite_boolean_model as bbm
import rxncon.venntastic.sets as venn
import rxncon.core.reaction as rxn
import rxncon.core.state as sta
import rxncon.core.specification as spec


class BoolNetSystem:
    def __init__(self, bipartite_boolean_model: bbm.BipartiteBooleanModel):
        self.bipartite_boolean_model = bipartite_boolean_model

    def to_string(self) -> str:
        bipartite_boolean_model_strings = [self._header_string(), self._build_init_conditions(), self._build_rules()]
        return "\n".join(bipartite_boolean_model_strings)

    def to_file(self, file_path: str):
        path, file = os.path.split(file_path)
        if path and os.path.exists(path):
            if not os.path.isfile(file_path):
                self._write_to_file(file_path)
            else:
                raise FileExistsError("{0} exists! remove file and run again".format(file_path))
        elif not path:
            if not os.path.isfile(file):
                self._write_to_file(file_path)
            else:
                print(os.path.dirname(file_path))
                raise FileExistsError("{0} exists! remove file and run again".format(file_path))
        elif path and not os.path.exists(path):
            raise NotADirectoryError("Path {0} does not exists.".format(path))

    def _write_to_file(self, file_path: str):
        with open(file_path, mode='w') as writehandle:
            writehandle.write(self.to_string())

    def _header_string(self):
        return "target, factors"

    def _build_rules(self):
        rules = [self._string_from_rule(rule) for rule in self.bipartite_boolean_model.rules]
        return "\n".join(rules)

    def _build_init_conditions(self):
        init_conditions = [self._string_from_init_conditions(init_condition)
                           for init_condition in self.bipartite_boolean_model.init_conditions]
        return "\n".join(init_conditions)

    def _string_from_init_conditions(self, init_condition: bbm.InitCondition):
        if init_condition.value is None:
            return "{0}, {0}".format(self._generate_name(init_condition.target))
        else:
            raise NotImplementedError

    def _string_from_rule(self, rule: bbm.Rule):
        return "{0}, {1}".format(self._target_to_string(rule.target),
                                 self._factor_to_string(rule.factor.simplified_form().value))

    def _target_to_string(self, target: bbm.Node):
        return self._generate_name(target)

    def _factor_to_string(self, factor: bbm.Factor):
        if isinstance(factor, venn.PropertySet):
            return self._generate_name(factor.value)
        elif isinstance(factor, venn.Complement):
            return "! {0}".format(self._factor_to_string(factor.expr))
        elif isinstance(factor, venn.Intersection):
            return '({0} & {1})'.format(self._factor_to_string(factor.left_expr),
                                        self._factor_to_string(factor.right_expr))
        elif isinstance(factor, venn.Union):
            return '({0} | {1})'.format(self._factor_to_string(factor.left_expr),
                                        self._factor_to_string(factor.right_expr))
        else:
            raise NotImplementedError

    def _generate_name(self, node: bbm.Node):
        if isinstance(node.value, rxn.Reaction):
            return string_from_reaction(node.value)

        elif isinstance(node.value, sta.State):
            return replace_invalid_chars(str(node.value))

        else:
            raise NotImplementedError


def string_from_reaction(reaction: rxn.Reaction) -> str:
    # todo: make a table out of this for marcus
    def generate_reaction_verb_name(verb: rxn.Verb):
        if re.search("-", verb.value):
            return re.sub("-", "minus", verb.value)
        elif re.search("\+", verb.value):
            return re.sub("\+", "plus", verb.value)
        else:
            return verb.value
    reaction_str = "{0}_{1}_{2}".format(reaction.subject, generate_reaction_verb_name(reaction.verb), reaction.object)
    return replace_invalid_chars(reaction_str)


def replace_invalid_chars(value):
    value = re.sub("-", "_", value)
    value = re.sub('[\[{(]', ".", value)
    value = re.sub('[\]})]', ".", value)
    value = re.sub('/', ".", value)
    return value
