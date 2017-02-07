from typing import Tuple, Dict

from rxncon.simulation.boolean.boolean_model import BooleanModel, ReactionTarget, KnockoutTarget, OverexpressionTarget, \
    StateTarget, UpdateRule
from rxncon.venntastic.sets import Set as VennSet, ValueSet, Complement, Intersection, Union, EmptySet, UniversalSet


def boolnet_from_boolean_model(boolean_model: BooleanModel) -> Tuple[str, Dict[str, str], Dict[str, bool]]:
    """
    Translates the boolean model into the BoolNet syntax.

    Note:
        BoolNet is an R package that provides tools for assembling, analyzing and visualizing Boolean networks.

    Args:
        boolean_model: The boolean model.

    Returns:
        1. The first return value is the boolean model in BooleNet syntax.
        2. The second return value is an abbreviation, target mapping.
        3. The third return value is a mapping of the initial condition.

    """
    def initialize_boolnet_names() -> None:
        nonlocal boolnet_names
        nonlocal reaction_index
        nonlocal state_index
        nonlocal knockout_index
        nonlocal overexpression_index
        for update_rule in boolean_model.update_rules:
            target = update_rule.target
            if target in boolnet_names.keys():
                raise AssertionError
            elif isinstance(target, ReactionTarget):
                name = 'R{}'.format(reaction_index)
                boolnet_names[target] = name
                reaction_index += 1
            elif isinstance(target, KnockoutTarget):
                name = 'K{}'.format(knockout_index)
                boolnet_names[target] = name
                knockout_index += 1
            elif isinstance(target, OverexpressionTarget):
                name = 'O{}'.format(overexpression_index)
                boolnet_names[target] = name
                overexpression_index += 1
            elif isinstance(target, StateTarget):
                name = 'S{}'.format(state_index)
                boolnet_names[target] = name
                state_index += 1
            else:
                raise AssertionError

    def str_from_factor(factor: VennSet) -> str:
        """
        Translates a factor into a string.

        Note:
            During this process the names of the targets are replaced by abbreviations. Reaction targets are replaced
            by R{} and states are replaced by S{} where {} is a continuous numerating for state and reaction targets
            respectively.

        Args:
            factor:

        Returns:
            The string of the factor.

        Raises:
            AssertionError: If the factor is not a valid VennSet object an error is raised.

        """
        if isinstance(factor, ValueSet):
            return boolnet_names[factor.value]
        elif isinstance(factor, Complement):
            return '!({})'.format(str_from_factor(factor.expr))
        elif isinstance(factor, Intersection):
            return '({})'.format(' & '.join(str_from_factor(x) for x in factor.exprs))
        elif isinstance(factor, Union):
            return '({})'.format(' | '.join(str_from_factor(x) for x in factor.exprs))
        elif isinstance(factor, EmptySet):
            return '0'
        elif isinstance(factor, UniversalSet):
            return '1'
        else:
            raise AssertionError('Could not parse factor {}'.format(factor))

    def str_from_update_rule(update_rule: UpdateRule) -> str:
        """
        Creates a string from an update rule.

        Args:
            update_rule: A target and its factor.

        Returns:
            The string of the update rule.

        """
        return '{0}, {1}'.format(boolnet_names[update_rule.target],
                                 str_from_factor(update_rule.factor))

    boolnet_names  = {}  # type: Dict[Target, str]

    reaction_index       = 0
    state_index          = 0
    knockout_index       = 0
    overexpression_index = 0

    initialize_boolnet_names()

    def sort_key(rule_str: str) -> Tuple[str, int]:
        """
        Function for sorting.

        Args:
            rule_str: A string representing and update rule of the boolean system.

        Returns:
            A tuple of string and integer.

        """
        target = rule_str.split(',')[0].strip()
        return target[0], int(target[1:])

    rule_strs = sorted([str_from_update_rule(x) for x in boolean_model.update_rules], key=sort_key)

    return 'targets, factors\n' + '\n'.join(rule for rule in rule_strs) + '\n', \
           {name: str(target) for target, name in boolnet_names.items()}, \
           {boolnet_names[target]: value for target, value in boolean_model.initial_conditions.target_to_value.items()}