from typing import List, Optional
import rxncon.simulation.rule_based.bngl_export as bngle

class RuleBasedModel:
    def __init__(self):
        self.rules = []  # type: List[Rule]

    def export_to_bngl(self):
        bngl_export = bngle.BNGLSystem(self.rules)
        bngl_export.bngl_file_from_rule_based_model()
        bngl_export.write_to_file()


class Rule:
    def __init__(self, reaction: object):
        self.name = None  # name of the reaction

        # if we have a unified complex we just have a left_reactant
        self.left_reactant = None  # type: Reactant
        self.right_reactant = None  # type: Reactant
        self.rate = 1.0     # type: float
        self.forward_rate = None
        self.reverse_rate = None

        self.reaction_class = reaction.reaction_class
        self.directionality = reaction.directionality
        # influence: positive or negative
        self.influence = reaction.influence
        self.isomerism = reaction.isomerism
        self.modifier = reaction.modifier

        self.source = reaction.source
        self.product = reaction.product

    def __eq__(self, other) -> bool:
        assert isinstance(other, Rule)

        if self.name == other.name and \
                        self.left_reactant == other.left_reactant and \
                        self.right_reactant == other.right_reactant:
            return True
        else:
            return False

    @property
    def rate(self):
        if self.directionality == 1:
            self.forward_rate = True
            self.reverse_rate = True
        elif self.directionality == 1:
            self.forward_rate = True
            self.reverse_rate = False

class Reactant:
    def __init__(self, name: str):
        self.states = []  # type: List[State]
        self.name = None

    def add_state(self, state: object, influence: bool):
        assert isinstance(state, State)
        if isinstance(influence, bool) or influence is None:
            state = State(state, influence)
            self.states.append(state)
        else:
            raise AttributeError()

    def remove_state(self, state: object, influence: bool):
        assert isinstance(state, State)
        state = State(state, influence)
        self.states.remove(state)

    def __eq__(self, other: Reactant) -> bool:
        assert isinstance(other, Reactant)

        if self.name == other.name and self.states == other.states:
            return True
        else:
            return False

    def __str__(self) -> str:
        return "Reactant({1},{2})".format(self.name, str(self.states))


class State:
    def __init__(self, state: object, influence: Optional[bool]):
        self.state = state
        #True = positive influence = !
        #False = negative influence = x
        self.influence = influence  # type: bool

    def __eq__(self, other) -> bool:
        assert isinstance(other, State)

        return self.state == other.state

    def __str__(self) -> str:
        return self.state.__str__()





