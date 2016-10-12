from typing import List
import re

from rxncon.core.reaction import BIDIRECTIONAL_VERBS

def split_bidirectional_reaction_str(rxn_str: str) -> List[str]:
    for verb in BIDIRECTIONAL_VERBS:
        if '_{}_'.format(verb).lower() in rxn_str.lower():
            verb = re.findall('(?i){}'.format(verb), rxn_str)[0]

            return [
                rxn_str.replace('_{}_'.format(verb), '_{}+_'.format(verb)),
                rxn_str.replace('_{}_'.format(verb), '_{}-_'.format(verb)),
            ]

    return [rxn_str]
