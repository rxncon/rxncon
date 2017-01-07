from enum import unique
from copy import deepcopy
import logging
from rxncon.util.utils import current_function_name

from rxncon.core.effector import Effector, StructEquivalences, QualSpec, StateEffector, TrivialStructEquivalences
from rxncon.core.reaction import Reaction
from rxncon.util.utils import OrderedEnum


logger = logging.getLogger(__name__)


@unique
class ContingencyType(OrderedEnum):
    undefined   = None
    requirement = '!'
    inhibition  = 'x'
    positive    = 'k+'
    negative    = 'k-'
    no_effect   = '0'
    unknown     = '?'


class Contingency:
    def __init__(self, target: Reaction, type: ContingencyType, effector: Effector):
        self.target, self.type, self.effector = target, type, effector

    def __eq__(self, other: 'Contingency') -> bool:
        return self.target == other.target and self.type == other.type and self.effector == other.effector

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'Contingency({0}, {1}, {2}'.format(str(self.target), str(self.type), str(self.effector))

    def clone(self) -> 'Contingency':
        return deepcopy(self)

    def to_structured(self):
        logger.debug('{}: {}'.format(current_function_name(), str(self)))

        if isinstance(self.effector, StateEffector) and self.effector.is_structured:
            # A fully structured StateEffector is fine.
            return self
        elif isinstance(self.effector, StateEffector) and not self.effector.is_structured:
            # For a non-structured StateEffector, assume the Specs appearing in the Effector
            # match those appearing in the Reaction.
            equivs = StructEquivalences()
            struct_components = {spec.to_non_struct_spec(): spec for spec in self.target.components_lhs_structured}
            for spec in self.effector.expr.specs:
                try:
                    equivs.add_equivalence(QualSpec([], struct_components[spec.to_component_spec()]),
                                           QualSpec([str(self.target)], spec.to_component_spec()))
                except KeyError:
                    pass
            sc = self.clone()
            sc.effector = sc.effector.to_merged_struct_effector(equivs, None, [str(self.target)])
            logger.info('{}: {} :: {} -> {}'.format(current_function_name(), str(self.target),
                                                    str(self.effector.expr), str(sc.effector.expr)))
            return sc
        elif self.effector.is_structured:
            # A fully structured Boolean Effector needs to have its structure indices merged.
            sc = self.clone()
            sc.effector = sc.effector.to_merged_struct_effector()
            logger.info('{}: {} :: {} -> {}'.format(current_function_name(), str(self.target),
                                                    str(self.effector), str(sc.effector)))
            return sc
        else:
            # For a non-structured Boolean Effector, assume all Specs that could match, actually do match.
            sc = self.clone()
            struct_components = {spec.to_non_struct_spec(): spec for spec in self.target.components_lhs_structured}
            sc.effector = sc.effector.to_merged_struct_effector(TrivialStructEquivalences(struct_components), None, None)
            logger.info('{}: {} :: {} -> {}'.format(current_function_name(), str(self.target),
                                                    str(self.effector), str(sc.effector)))
            return sc


