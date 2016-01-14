

class RuleBasedModel:
    def export_to_bngl(self):
        bngl = bngl_export.bngl_file_from_rule_based_model(self)
        bngl.export()



class Rule:
    def __init__(self, dsfafda):
        self.first_reactant = None
        self.second_reactant = None
        self.reaction_constant = 1.0



class Reactant:
    def __init__(self):
        self.states = []




