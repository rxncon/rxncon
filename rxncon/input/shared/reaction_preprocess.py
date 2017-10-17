from typing import List
import re
import logging

from rxncon.core.reaction import BIDIRECTIONAL_REACTIONS

logger = logging.getLogger(__name__)


def split_bidirectional_reaction_str(rxn_str: str) -> List[str]:
    for verb in BIDIRECTIONAL_REACTIONS:
        if '_{}_'.format(verb).lower() in rxn_str.lower():
            verb = re.findall('(?i)_{}_'.format(verb), rxn_str)[0][1:-1]

            logger.info('split_bidirectional_reaction_str: {}'.format(rxn_str))

            return [
                rxn_str.replace('_{}_'.format(verb), '_{}+_'.format(verb)),
                rxn_str.replace('_{}_'.format(verb), '_{}-_'.format(verb)),
            ]

    return [rxn_str]
