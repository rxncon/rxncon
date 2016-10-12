from typing import List

from rxncon.core.reaction import BIDIRECTIONAL_VERBS

def preprocessed_reaction_strs(rxn_str: str) -> List[str]:
    for verb in BIDIRECTIONAL_VERBS:
        if '_{}_'.format(verb) in rxn_str:
            return [
                rxn_str.replace('_{}_'.format(verb), '_{}+_'.format(verb)),
                rxn_str.replace('_{}_'.format(verb), '_{}-_'.format(verb)),
            ]

    return [rxn_str]