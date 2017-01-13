from enum import unique, Enum
from copy import deepcopy
from collections import defaultdict
import logging
from rxncon.util.utils import current_function_name

from rxncon.core.effector import Effector, StructEquivalences, QualSpec, StateEffector, TrivialStructEquivalences, NotEffector, \
    AndEffector, OrEffector
from rxncon.core.reaction import Reaction
from rxncon.venntastic.sets import Set, ValueSet, Intersection, Complement, Union, UniversalSet


logger = logging.getLogger(__name__)


@unique
class ContingencyType(Enum):
    requirement = '!'
    inhibition  = 'x'
    positive    = 'k+'
    negative    = 'k-'
    no_effect   = '0'
    unknown     = '?'


class Contingency:
    def __init__(self, target: Reaction, type: ContingencyType, effector: Effector) -> None:
        self.target, self.type, self.effector = target, type, effector

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Contingency):
            return NotImplemented
        return self.target == other.target and self.type == other.type and self.effector == other.effector

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'Contingency({0}, {1}, {2}'.format(str(self.target), str(self.type), str(self.effector))

    def clone(self) -> 'Contingency':
        return deepcopy(self)

    def to_structured(self) -> 'Contingency':
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
            assert isinstance(sc.effector, StateEffector)
            logger.info('{}: {} :: {} -> {}'.format(current_function_name(), str(self.target),
                                                    str(self.effector.expr), str(sc.effector.expr)))
            sc._validate_structure_indices()
            return sc
        elif self.effector.is_structured:
            # A fully structured Boolean Effector needs to have its structure indices merged.
            sc = self.clone()
            sc.effector = sc.effector.to_merged_struct_effector()
            logger.info('{}: {} :: {} -> {}'.format(current_function_name(), str(self.target),
                                                    str(self.effector), str(sc.effector)))
            sc._validate_structure_indices()
            return sc
        else:
            # For a non-structured Boolean Effector, assume all Specs that could match, actually do match.
            sc = self.clone()
            struct_components = {spec.to_non_struct_spec(): spec for spec in self.target.components_lhs_structured}
            sc.effector = sc.effector.to_merged_struct_effector(TrivialStructEquivalences(struct_components), None, None)
            logger.info('{}: {} :: {} -> {}'.format(current_function_name(), str(self.target),
                                                    str(self.effector), str(sc.effector)))
            sc._validate_structure_indices()
            return sc

    def to_venn_set(self, k_plus_strict=False, k_minus_strict=False, structured=True, state_wrapper=lambda x: x) -> Set:
        def parse_effector(eff: Effector) -> Set:
            if isinstance(eff, StateEffector):
                if structured:
                    return ValueSet(state_wrapper(eff.expr))
                else:
                    return ValueSet(state_wrapper(eff.expr.to_non_structured()))
            elif isinstance(eff, NotEffector):
                return Complement(parse_effector(eff.expr))
            elif isinstance(eff, OrEffector):
                return Union(*(parse_effector(x) for x in eff.exprs))
            elif isinstance(eff, AndEffector):
                return Intersection(*(parse_effector(x) for x in eff.exprs))
            else:
                raise AssertionError

        if k_plus_strict:
            positive = (ContingencyType.requirement, ContingencyType.positive)
        else:
            positive = (ContingencyType.requirement,)

        if k_minus_strict:
            negative = (ContingencyType.inhibition, ContingencyType.negative)
        else:
            negative = (ContingencyType.inhibition,)

        if self.type in positive:
            return parse_effector(self.effector)
        elif self.type in negative:
            return Complement(parse_effector(self.effector))
        else:
            return UniversalSet()

    def _validate_structure_indices(self):
        # Assert that every index is only used once.
        specs = [spec for state in self.effector.states for spec in state.specs]
        index_to_specs = defaultdict(set)

        for spec in specs:
            assert spec.is_structured
            index_to_specs[spec.struct_index].add(spec.to_component_spec())

        assert all(len(x) == 1 for _, x in index_to_specs.items())



