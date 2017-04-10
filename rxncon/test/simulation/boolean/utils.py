from rxncon.core.reaction import reaction_from_str
from rxncon.core.spec import spec_from_str
from rxncon.core.state import state_from_str
from rxncon.simulation.boolean.boolean_model import Target, StateTarget, ReactionTarget, ComponentStateTarget


def target_from_str(target_str: str) -> Target:
    """
    Generates a target from string input.

    Args:
        target_str: The string representation of a StateTarget, ReactionTarget or ComponentStateTarget.

    Returns:
        A Target object e.g. StateTarget, ReactionTarget or ComponentStateTarget

    Raises:
        SyntaxError: If the string does not correspond to a predefined Target object an error is raised.

    """
    try:
        return StateTarget(state_from_str(target_str))
    except SyntaxError:
        pass

    try:
        if '#' in target_str:
            rxn_str, index_strs = target_str.split('#')
            target = ReactionTarget(reaction_from_str(rxn_str))
            contingency_variant_index, interaction_variant_index = None, None

            for variant_str in index_strs.split('/'):
                if variant_str[0] == 'c':
                    contingency_variant_index = int(variant_str[1:])
                elif variant_str[0] == 'i':
                    interaction_variant_index = int(variant_str[1:])

            target.contingency_variant_index = contingency_variant_index
            target.interaction_variant_index = interaction_variant_index
            return target
        else:
            return ReactionTarget(reaction_from_str(target_str))
    except SyntaxError:
        pass

    try:
        return ComponentStateTarget(spec_from_str(target_str))
    except SyntaxError:
        raise SyntaxError('Could not parse target str {}'.format(target_str))