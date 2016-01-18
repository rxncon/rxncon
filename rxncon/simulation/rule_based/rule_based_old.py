from typing import List, Optional
import rxncon.simulation.rule_based.bngl_export as bngle

class RuleBasedModel:
    def __init__(self):
        self.rules = []  # type: List[Rule]

    def export_to_bngl(self):
        # @basti My fault for putting it like this: BNGL_export should depend on rule_based instead of the other way around.
        #        So this method should go. I was trying to illustrate that the data in this class should be sufficient to
        #        create a full BNGL file.
        #        Think of it like this: the Rule based model is the "real thing", the BNGL file is merely a "representation"
        #        of it. The representation should depend on the real thing, not vice-versa.
        bngl_export = bngle.BNGLSystem(self.rules)
        bngl_export.bngl_file_from_rule_based_model()
        bngl_export.write_to_file()


class Rule:
    # @basti Isn't the 'rule' sort of like a 'reaction'? I think therefore all the stuff you are setting to 'none' here,
    #        such as left_reactant, rate, etc. should be in the function header of the __init__ method. And the 'reaction'
    #        in the header should then probably go.
    #
    #        In any case, this should _not_ yet depend on the rxncon Reaction object (I dunno if you meant that)
    #        We will have a separate translation step from rxncon -> rule_based that will create one or more Rule objects
    #        from each Reaction object.
    def __init__(self, reaction: object):
        self.name = None  # name of the reaction

        # if we have a unified complex we just have a left_reactant
        self.left_reactant = None  # type: Reactant
        self.right_reactant = None  # type: Reactant
        self.rate = 1.0     # type: float
        self.forward_rate = None
        self.reverse_rate = None

        # @basti I think this is not applicable for rule based models right? Haven't they just got reactants, a rate
        #        and product states?
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
        # @basti If you want to do something like this, make 'directionality' into an Enum. That way you can name the
        #        different options.
        if self.directionality == 1:
            self.forward_rate = True
            self.reverse_rate = True
        elif self.directionality == 1:
            self.forward_rate = True
            self.reverse_rate = False

class Reactant:
    # @basti I think you should ditch the add_state and remove_state methods. Keep it simpler and make sure that once a
    #        Reactant object is initialized, it is valid. This would mean that the constructor would look something like:
    #
    # def __init__(self, name:str, states: List[State])
    #     bla

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
        # @basti Why does the state need to know the influence? This is part of the Rule I guess?
        self.state = state
        #True = positive influence = !
        #False = negative influence = x
        self.influence = influence  # type: bool

    def __eq__(self, other) -> bool:
        assert isinstance(other, State)

        return self.state == other.state

    def __str__(self) -> str:
        return self.state.__str__()





