from enum import unique, Enum
from copy import deepcopy
from collections import defaultdict
import logging
from typing import Any, Callable, Set as TgSet, Dict, List, Optional  # pylint: disable=unused-import

from rxncon.util.utils import current_function_name
from rxncon.core.effector import Effector, StructEquivalences, QualSpec, StateEffector, TrivialStructEquivalences, NotEffector, \
    AndEffector, OrEffector, StructCounter
from rxncon.core.reaction import Reaction
from rxncon.venntastic.sets import Set as VennSet, ValueSet, Intersection, Complement, Union, UniversalSet
from rxncon.core.state import State


LOGGER = logging.getLogger(__name__)


@unique
class ContingencyType(Enum):
    requirement = '!'
    inhibition  = 'x'
    positive    = 'k+'
    negative    = 'k-'
    no_effect   = '0'
    unknown     = '?'


class Contingency:
    def __init__(self, reaction: Reaction, contingency_type: ContingencyType, effector: Effector) -> None:
        self.reaction, self.contingency_type, self.effector = reaction, contingency_type, effector

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Contingency):
            return NotImplemented
        return self.reaction == other.reaction and self.contingency_type == other.contingency_type and self.effector == other.effector

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return 'Contingency({0}, {1}, {2}'.format(str(self.reaction), str(self.contingency_type), str(self.effector))

    def clone(self) -> 'Contingency':
        return deepcopy(self)

    def with_merged_struct_effector(self, equivs: Optional[StructEquivalences]=None, counter: Optional[StructCounter]=None,
                                    namespace: Optional[List[str]]=None) -> 'Contingency':
        structured = self.clone()
        structured.effector = structured.effector.to_merged_struct_effector(equivs, counter, namespace)
        structured.validate_struct_indices()
        return structured

    def to_structured(self) -> 'Contingency':
        LOGGER.debug('{}: {}'.format(current_function_name(), str(self)))

        if isinstance(self.effector, StateEffector) and self.effector.is_structured:
            # A fully structured StateEffector is fine.
            return self
        elif isinstance(self.effector, StateEffector) and not self.effector.is_structured:
            # For a non-structured StateEffector, assume the Specs appearing in the Effector
            # match those appearing in the Reaction.
            equivs = StructEquivalences()
            struct_components = {spec.to_non_struct_spec(): spec for spec in self.reaction.components_lhs_structured}
            for spec in self.effector.expr.specs:
                try:
                    equivs.add_equivalence(QualSpec([], struct_components[spec.to_component_spec()]),
                                           QualSpec([str(self.reaction)], spec.to_component_spec()))
                except KeyError:
                    pass

            return self.with_merged_struct_effector(equivs, None, [str(self.reaction)])
        elif self.effector.is_structured:
            # A fully structured Boolean Effector needs to have its structure indices merged.
            return self.with_merged_struct_effector()
        else:
            # For a non-structured Boolean Effector, assume all Specs that could match, actually do match.
            struct_components = {spec.to_non_struct_spec(): spec for spec in self.reaction.components_lhs_structured}
            equivs = TrivialStructEquivalences(struct_components)  # pylint: disable=redefined-variable-type
            return self.with_merged_struct_effector(equivs)

    def to_venn_set(self, k_plus_strict: bool=False, k_minus_strict: bool=False, structured: bool=True,
                    state_wrapper: Callable[[State], Any]=lambda x: x) -> VennSet[Any]:
        def parse_effector(eff: Effector) -> VennSet:
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
            positive = (ContingencyType.requirement,)  # type: ignore

        if k_minus_strict:
            negative = (ContingencyType.inhibition, ContingencyType.negative)
        else:
            negative = (ContingencyType.inhibition,)  # type: ignore

        if self.contingency_type in positive:
            return parse_effector(self.effector)
        elif self.contingency_type in negative:
            return Complement(parse_effector(self.effector))
        else:
            return UniversalSet()

    def validate_struct_indices(self) -> None:
        # Assert that every index is only used once.
        specs = [spec for state in self.effector.states for spec in state.specs]
        index_to_specs = defaultdict(set)  # type: Dict[int, TgSet]

        for spec in specs:
            assert spec.struct_index is not None
            index_to_specs[spec.struct_index].add(spec.to_component_spec())

        assert all(len(x) == 1 for _, x in index_to_specs.items())



